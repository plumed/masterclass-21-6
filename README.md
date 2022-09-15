# PLUMED Masterclass 21.6: Dimensionality reduction

## Aims

The primary aim of this Masterclass is to show you how you might use PLUMED in your work.
We will see how to call PLUMED from a python notebook and discuss some strategies for selecting
collective variables.

## Objectives

Once this Masterclass is completed, users will be able to:

- Calculate CVs and import them into a python notebook.
- Generate visualizations of data using chemiscope
- Use the projection of a vector as a CV.
- Use path collective variables.
- Run dimensionality reduction algorithms with PLUMED
- Use Gibbs dividing surfaces to determine the extent of a phase. 

## Acknowledgements

Throughout this exercise, we use the [atomistic simulation environment](https://wiki.fysik.dtu.dk/ase/) 
and [chemiscope](https://chemiscope.org/).  Please look at the information at the links I have provided here 
for more information about these codes.

## Setting up the software 

For this Masterclass, you will need to set up the plumed and gromacs as you did for [this masterclass](https://www.plumed.org/doc-v2.8/user-doc/html/masterclass-21-1.html).  You thus install plumed and gromacs using:

````
conda install --strict-channel-priority -c plumed/label/masterclass-mpi -c conda-forge plumed
conda install --strict-channel-priority -c plumed/label/masterclass-mpi -c conda-forge gromacs
````

You can install the other software you need using:

````
conda install -c conda-forge py-plumed Numpy pandas matplotlib notebook mdtraj mdanalysis git ase
````

Notice that you need a package called ase (the atomic simulation environment) for these exercises and the other packages you have been using.
You also need to install chemiscope, which you can do by using the following command:

````
pip install chemiscope
````

## Resources

The data needed to complete this Masterclass can be found on [GitHub](https://github.com/plumed/masterclass-21-6).
You can clone this repository locally on your machine using the following command:

````
git clone https://github.com/plumed/masterclass-21-6.git
````

I recommend that you run each exercise in a separate sub-directory (i.e. Exercise-1, Exercise-2, ...), which you can create inside the root directory `masterclass-21-6`. Organizing your data this way will help you to keep things clean.

_All the exercises have been tested with PLUMED version 2.7.0._

## Exercises

### Exercise 1: Running PLUMED from a python notebook

I like working with Python notebooks.  Using notebooks allows me to keep the codes that I am working on close to the explanations of how these codes work.
I can also make all my figures within a notebook.  By using notebooks, I can thus have a single file that contains:

- Discussion of the work done.
- The analysis code.
- The final figures that were generated.

In previous masterclasses we have run plumed driver through bash from within a notebook by using commands like the one shown below:

```python
!cd ../Exercises/Exercise_1 && plumed driver --noatoms > /dev/null
# Read in colvar file produced
data = np.loadtxt("../Exercises/Exercise_1/colvar")
```

We have then read in the colvar file produced and analyzed it from within the notebook.  We can avoid using plumed driver and can call 
plumed directly from python using the python interface.  The code below, for instance, generates the plot underneath the 
code.  The code reads in the trajectory from the file traj.pdb that you obtained from the GitHub repository.  The plot then shows how 
the $\phi$ angle on the second residue of the protein changes with time.     
 
```python
import matplotlib.pyplot as plt
import numpy as np
import plumed
import ase 
import ase.io 

# Read in trajectory using ase
traj = ase.io.read('../data/traj.pdb',':')

# Setup plumed object to do calculation
p = plumed.Plumed()
p.cmd("setMDEngine","python")
# Read PDB so need to multiply by 0.1 to convert to nm
p.cmd("setMDLengthUnits", 0.1)
p.cmd("setTimestep", 1.)
p.cmd("setKbT", 1.)
natoms = len(traj[0].positions)
p.cmd("setNatoms",natoms)
p.cmd("setLogFile","test.log")
p.cmd("init")

# Read plumed input 
p.cmd("readInputLine","MOLINFO STRUCTURE=../data/bhp.pdb")
# If you are doing many variables I would represent putting these 
# next three PLUMED commands into a function
p.cmd("readInputLine", "t1: TORSION ATOMS=@phi-2" )
# Now setup some memory to hold the variable that is shared 
# between plumed and the underlying code
shape = np.zeros( 1, dtype=np.int_ )
p.cmd("getDataRank t1 ", shape )
t1 = np.zeros((1))
p.cmd("setMemoryForData t1", t1)

# Loop over trajectory and get data from plumed
nfram, tt, v1, box = 0, [], [], np.array([[100.,0,0],[0,100.,0],[0,0,100]])
charges, forces, virial = np.zeros(natoms,dtype=np.float64), np.zeros([natoms,3]), np.zeros((3,3),dtype=np.float64)
for ts in traj : 
    # Set all the input variables to PLUMED
    p.cmd("setStep",nfram)
    p.cmd("setBox",box )
    p.cmd("setMasses", ts.get_masses() )
    p.cmd("setCharges", charges )
    pos = np.array(ts.get_positions(), dtype=np.float64 )
    p.cmd("setPositions", pos )
    p.cmd("setForces", forces )
    p.cmd("setVirial", virial )
    # Run the plumed calculation
    p.cmd("calc")
    tt.append(nfram)
    # We can now extract the value of the torsion by accessing the shared memory we set up earlier
    v1.append(t1[0])
    nfram = nfram + 1

# Plot teh graph of the torsional angle as a function of time
plt.plot( tt, v1, 'ko')
plt.xlabel("Simulation step")
plt.ylabel("Torsion angle / radians")
plt.show()
```

![The value of the phi angle as a function of time.](figures/masterclass-21-6-phi-time.png)

__Your task in this first exercise is to modify the code above and to produce a figure similar to the one shown below.__  This figure 
shows all the values of the \f$\phi\f$ and \f$\psi\f$ angles in the second residue of the protein during the simulation.

![A plot showing phi against psi for the second residue of the protein.](figures/masterclass-21-6-phi-psi.png)

### Exercise 2: Generating a chemiscope representation

Plots showing the trajectory in CV space similar to those you generated at the end of the previous exercise are helpful.  What would be more 
useful, however, is some way of understanding the relationship between the positions in CV space and the structures of the various atoms.
In other words, what we would like is something like this:

![A representation of the time against phi plot from the previous exercise that was generated using chemiscope.](masterclass-21-6-chemiscope.png)

You can see the frame in the trajectory that each point in the plot corresponds to from the above figure.  The snapshots on the right correspond to the structures the system
had at the points highlighted in red, yellow, green and blue respectively in the plot on the left.  

The figure above was generated using chemiscope.  
This server allows you to generate and interact with plots like the one shown above.  __Your task in this exercise is to generate your own chemiscope representation
of the data in traj.pdb.__  To create a chemiscope representation of the $\phi$ angles that we generated using the first python script from the previous exercise, you would 
add the following python code:

```python
from ase.data import atomic_masses
from chemiscope import write_input

# This ensures that the atomic masses are used in place of the symbols
# when constructing the atomic configurations' chemiscope representations.
# Using the symbols will not work because ase is written by chemists and not 
# biologists.  For a chemist, HG1 is mercury as opposed to the first hydrogen
# on a guanine residue.  
for frame in traj:
    frame.numbers = np.array(
        [
            np.argmin(np.subtract(atomic_masses, float(am)) ** 2)
            for am in frame.arrays["occupancy"]
        ]
    )

# This constructs the dictionary of properties for chemiscope
properties = {
    "time": {
        "target": "structure",
        "values": tt,
        "description": "Simulation step number",
    },
    "t2": {
        "target": "structure",
        "values": v1,
        "description": "Phi angle of second residue",
    },
}

# This generates our chemiscope output
write_input("torsion_chemiscope.json.gz", frames=traj, properties=properties )
```

You would then upload the `torsion_chemiscope.json.gz` file that is generated by this script at [https://chemiscope.org](https://chemiscope.org).

__See if you can generate your own chemiscope representation of the data in traj.pdb.__ I would recommend calculating and uploading a chemiscope representation of all the protein's torsional angles. At the very least, you need to do at least two backbone torsional angles.  However, if you do more than two torsions, you can generate a plot like the one shown below.

!["A representation of the trajectory that was generated using chemiscope.  Here psi_2 is plotted on the x-axis, phi_2 is plotted on the y-axis, and phi_3 is plotted on the z-axis.  The points are coloured following the value of psi_3, and the points' sizes are proportional to phi_4."](figures/masterclass-21-6-chemiscope2.png)

### Exercise 3: Dimensionality reduction

Chemiscope comes into its own when you are working with a machine learning algorithm.  These algorithms can (in theory) learn the collective variables you need to use from the trajectory data.
To make sense of the coordinates that have been learned, you have to carefully visualize where structures are projected in the low dimensional space.  You can use chemiscope to complete this 
process of visualizing the representation the computer has found for you.  In this next set of exercises, we will apply various dimensionality reduction algorithms to the data 
contained in the file traj.pdb.  If you visualize the contents of that file using VMD, you will see that this file contains a short protein trajectory.  You are perhaps unsure what 
CV to analyze this data and thus want to see if you can shed any light on the contents of the trajectory by using machine learning.

Typically, PLUMED analyses one set of atomic coordinates at a time.  To run a machine learning algorithm, however, you need to gather information on multiple configurations.
Therefore, the first thing you need to learn to use these algorithms is how to store configurations for later analysis with a machine learning algorithm.  The following input
illustrates how to complete this task using PLUMED.
        
```plumed
#SOLUTIONFILE=work/plumed_ex1.dat
# This reads in the template pdb file and thus allows us to use the @nonhydrogens
# special group later in the input
MOLINFO STRUCTURE=__FILL__ MOLTYPE=protein

# This stores the positions of all the nonhydrogen atoms for later analysis
cc: COLLECT_FRAMES __FILL__=@nonhydrogens

# This should output the atomic positions for the frames that were collected to a pdb file called traj.pdb
OUTPUT_ANALYSIS_DATA_TO_PDB USE_OUTPUT_DATA_FROM=__FILL__ FILE=traj.pdb
```

__Copy the input above into a plumed file and fill in the blanks.__  You should then be able to run the command using:

Then, once all the blanks are filled in, run the command using:

````
plumed driver --mf_pdb traj.pdb
````

You can also store the values of collective variables for later analysis with these algorithms.  __Modify the input above so that all
Thirty backbone dihedral angles in the protein are stored, and output using [OUTPUT_ANALYSIS_DATA_TO_COLVAR](https://www.plumed.org/doc-master/user-doc/html/_o_u_t_p_u_t__a_n_a_l_y_s_i_s__d_a_t_a__t_o__c_o_l_v_a_r.html) and rerun the calculation.__

You can find more information on the dimensionality reduction algorithms that we are using in this section in [this paper](https://arxiv.org/abs/1907.04170).

#### PCA

Having learned how to store data for later analysis with a dimensionality reduction algorithm, lets now apply principal component analysis (PCA) upon
our stored data.  In principle component analysis, low dimensional projections for our trajectory are constructed by:

- Computing a covariance matrix from the trajectory data
- Diagonalizing the covariance matrix.
- Calculating the projection of each trajectory frame on a subset of the eigenvectors of the covariance matrix.

To perform PCA using PLUMED, we are going to use the following input with the blanks filled in: 

```plumed
#SOLUTIONFILE=work/plumed_ex2.dat
# This reads in the template pdb file and thus allows us to use the @nonhydrogens
# special group later in the input
MOLINFO STRUCTURE=__FILL__ MOLTYPE=protein

# This stores the positions of all the nonhydrogen atoms for later analysis
cc: COLLECT_FRAMES __FILL__=@nonhydrogens
# This diagonalizes the covariance matrix
pca: PCA USE_OUTPUT_DATA_FROM=__FILL__ METRIC=OPTIMAL NLOW_DIM=2
# This projects each of the trajectory frames onto the low dimensional space that was
# identified by the PCA command
dat: PROJECT_ALL_ANALYSIS_DATA USE_OUTPUT_DATA_FROM=__FILL__ PROJECTION=__FILL__

# This should output the PCA projections of all the coordinates
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=__FILL__ ARG=dat.* FILE=pca_data

# These next three commands calculate the secondary structure variables.  These
# variables measure how much of the structure resembles an alpha helix, an antiparallel beta-sheet
# and a parallel beta-sheet.  Configurations that have different secondary structures should be projected
# in different parts of the low dimensional space.
alpha: ALPHARMSD RESIDUES=all
abeta: ANTIBETARMSD RESIDUES=all STRANDS_CUTOFF=1.0
pbeta: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1.0

# These commands collect and output the secondary structure variables so that we can use this information to
# determine how good our projection of the trajectory data is.
cc2: COLLECT_FRAMES ARG=alpha,abeta,pbeta
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=cc2 ARG=cc2.* FILE=secondary_structure_data
```

To generate the projection, you run the command:

````
plumed driver --mf_pdb traj.pdb
````

You can generate a projection of the data above using chemiscope by using the following script:

```python
# This ase command should read in the traj.pdb file that was analyzed.  Notice that the analysis
# actions below ignore the first frame in this trajectory, so we need to take that into account
# when we generated the chemiscope
traj = ase.io.read('../data/traj.pdb',':')
# This reads in the PCA projection that are output by by the OUTPUT_ANALYSIS_DATA_TO_COLVAR command
# above
projection = np.loadtxt("pca_data")
# We also read in the secondary structure data by colouring points following the secondary
# structure. We can get a sense of how good our projection is.
structure = np.loadtxt("secondary_structure_data")

# This ensures that the atomic masses are used in place of the symbols
# when constructing the atomic configurations' chemiscope representations.
# Using the symbols will not work because ase is written by chemists and not
# biologists.  For a chemist, HG1 is mercury as opposed to the first hydrogen
# on a guanine residue.
for frame in traj:
    frame.numbers = np.array(
        [
            np.argmin(np.subtract(atomic_masses, float(am)) ** 2)
            for am in frame.arrays["occupancy"]
        ]
    )

# This constructs the dictionary of properties for chemiscope
properties = {
    "pca1": {
        "target": "structure",
        "values": projection[:,0],
        "description": "First principle component",
    },
    "pca2": {
        "target": "structure",
        "values": projection[:,1],
        "description": "Second principle component",
    },
    "alpha": {
        "target": "structure",
        "values": structure[:,0],
        "description": "Alpha helical content",
    },
    "antibeta": {
        "target": "structure",
        "values": structure[:,1],
        "description": "Anti parallel beta sheet content",
    },
    "parabeta": {
        "target": "structure",
        "values": structure[:,2],
        "description": "Parallel beta sheet content",
    },
}

# This generates our chemiscope output
write_input("pca_chemiscope.json.gz", frames=traj[1:], properties=properties )
```

When the output from this set of commands is loaded into chemiscope, we can construct figures like the one shown below.  On the axes here, we have plotted the PCA coordinates.  The
points are then coloured according to the alpha-helical content.

![A representation of the PCA projection that was generated using chemiscope.](figures/masterclass-21-6-chemiscope-pca.png)

__See if you can use PLUMED and chemiscope to generate a figure similar to the one above.__  Try to experiment with the way the points are coloured.  Look at the beta-sheet content as well.

#### MDS

In the previous section, we performed PCA on the atomic positions directly.  In the section before last, however, we also saw how we could store high-dimensional vectors of collective variables and then
use them as input to a dimensionality reduction algorithm.  Therefore, we might legitimately ask if we can do PCA using these high-dimensional vectors as input rather than atomic positions.
The answer to this question is yes as long as the CV is not periodic.  If any of our CVs are not periodic, we cannot analyze them using the PCA action.  We can, however, formulate the PCA algorithm
differently.  In this alternative formulation, which is known as classical multidimensional scaling (MDS), we do the following:

- We calculate the matrix of distances between configurations
- We perform an operation known as centring the matrix.
- We diagonalize the centred matrix
The eigenvectors multiplied by the corresponding eigenvalue's square root can then be used as a set of projections for our input points.

This method is used less often than PCA as the matrix that we have to diagonalize here in the third step can be considerably larger than the matrix that we have to diagonalize when we perform PCA.  
To avoid this expensive diagonalization step, we often select a subset of so-called landmark points on which to run the algorithm directly.  Projections for the remaining points are then found
by using a so-called out-of-sample procedure.  This is what has been done in the following input:

```plumed
#SOLUTIONFILE=work/plumed_ex3.dat
# This command reads in the template pdb file and thus allows us to use the @nonhydrogens
# group later in the input
MOLINFO STRUCTURE=__FILL__ MOLTYPE=protein

# This stores the positions of all the nonhydrogen atoms for later analysis
cc: COLLECT_FRAMES ATOMS=@nonhydrogens

# The following commands compute all the Ramachandran angles of the protein for you
r2-phi: TORSION ATOMS=@phi-2
r2-psi: TORSION ATOMS=@psi-2
r3-phi: TORSION ATOMS=@phi-3
r3-psi: TORSION ATOMS=@psi-3
r4-phi: TORSION __FILL__
r4-psi: TORSION __FILL__
r5-phi: TORSION __FILL__
r5-psi: TORSION __FILL__
r6-phi: TORSION __FILL__
r6-psi: TORSION __FILL__
r7-phi: TORSION __FILL__
r7-psi: TORSION __FILL__
r8-phi: TORSION __FILL__
r8-psi: TORSION __FILL__
r9-phi: TORSION __FILL__
r9-psi: TORSION __FILL__
r10-phi: TORSION __FILL__
r10-psi: TORSION __FILL__
r11-phi: TORSION __FILL__
r11-psi: TORSION __FILL__
r12-phi: TORSION __FILL__
r12-psi: TORSION __FILL__
r13-phi: TORSION __FILL__
r13-psi: TORSION __FILL__
r14-phi: TORSION __FILL__
r14-psi: TORSION __FILL__
r15-phi: TORSION __FILL__
r15-psi: TORSION __FILL__
r16-phi: TORSION __FILL__
r16-psi: TORSION __FILL__

# This command stores all the Ramachandran angles that were computed
angles: COLLECT_FRAMES __FILL__=r2-phi,r2-psi,r3-phi,r3-psi,r4-phi,r4-psi,r5-phi,r5-psi,r6-phi,r6-psi,r7-phi,r7-psi,r8-phi,r8-psi,r9-phi,r9-psi,r10-phi,r10-psi,r11-phi,r11-psi,r12-phi,r12-psi,r13-phi,r13-psi,r14-phi,r14-psi,r15-phi,r15-psi,r16-phi,r16-psi
# Lets now compute the matrix of distances between the frames in the space of the Ramachandran angles
distmat: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=__FILL__ METRIC=EUCLIDEAN
# Now select 500 landmark points to analyze
fps: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=__FILL__ NLANDMARKS=500
# Run MDS on the landmarks
mds: CLASSICAL_MDS __FILL__=fps NLOW_DIM=2
# Project the remaining trajectory data
osample: PROJECT_ALL_ANALYSIS_DATA USE_OUTPUT_DATA_FROM=__FILL__ PROJECTION=__FILL__

# This command outputs all the projections of all the points in the low dimensional space
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=__FILL__ ARG=osample.* FILE=mds_data

# These next three commands calculate the secondary structure variables.  These
# variables measure how much of the structure resembles an alpha helix, an antiparallel beta-sheet
# and a parallel beta-sheet.  Configurations that have different secondary structures should be projected
# in different parts of the low dimensional space.
alpha: ALPHARMSD RESIDUES=all
abeta: ANTIBETARMSD RESIDUES=all STRANDS_CUTOFF=1.0
pbeta: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1.0

# These commands collect and output the secondary structure variables so that we can use this information to
# determine how good our projection of the trajectory data is.
cc2: COLLECT_FRAMES ARG=alpha,abeta,pbeta
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=cc2 ARG=cc2.* FILE=secondary_structure_data
```

This input collects all the torsional angles for the configurations in the trajectory.  Then, at the end of the calculation, the matrix of distances between these points is computed, and a set of landmark points
is selected using a method known as farthest point sampling.  A matrix that contains only those distances between the landmarks is then constructed and diagonalized by the CLASSICAL_MDS action so that
projections of the landmarks can be built.  The final step is then to project the remainder of the trajectory using the PROJECT_ALL_ANALYSIS_DATA action.  Try to fill in the blanks in the input above
and run this calculation now using the command:

````
plumed driver --mf_pdb traj.pdb
````

Once the calculation has completed, you can, once again, visualize the data using chemiscope by using a suitably modified version of the script from the previous exercise.  The image below shows the
MDS coordinates coloured according to the alpha-helical content.

![A representation of the MDS projection that was generated using chemiscope.](figures/masterclass-21-6-chemiscope-mds.png)

__Try to generate an image that looks like this one yourself by completing the input above and then using what you learned about generating chemiscope representations in the previous exercise.__

#### Sketch-map

The two algorithms (PCA and MDS) that we have looked at thus far are linear dimensionality reduction algorithms.  In addition to these, there are a whole class of non-linear dimensionality reduction
reduction algorithms which work by transforming the matrix of dissimilarities between configurations, calculating geodesic rather than Euclidean distances between configurations or by changing the form of the
loss function that is optimized.  In this final exercise, we will use an algorithm that uses the last of these three strategies to construct a non-linear projection.  The algorithm is known as sketch-map
and input for sketch-map is provided below:

```plumed
#SOLUTIONFILE=work/plumed_ex4.dat
# This reads in the template pdb file and thus allows us to use the @nonhydrogens
# special group later in the input
MOLINFO STRUCTURE=__FILL__ MOLTYPE=protein

# This stores the positions of all the nonhydrogen atoms for later analysis
cc: COLLECT_FRAMES __FILL__=@nonhydrogens
# This should output the atomic positions for the frames that were collected and analyzed using MDS
OUTPUT_ANALYSIS_DATA_TO_PDB USE_OUTPUT_DATA_FROM=__FILL__ FILE=traj.pdb

# The following commands compute all the Ramachandran angles of the protein for you
r2-phi: TORSION ATOMS=@phi-2
r2-psi: TORSION ATOMS=@psi-2
r3-phi: TORSION ATOMS=@phi-3
r3-psi: TORSION ATOMS=@psi-3
r4-phi: TORSION __FILL__
r4-psi: TORSION __FILL__
r5-phi: TORSION __FILL__
r5-psi: TORSION __FILL__
r6-phi: TORSION __FILL__
r6-psi: TORSION __FILL__
r7-phi: TORSION __FILL__
r7-psi: TORSION __FILL__
r8-phi: TORSION __FILL__
r8-psi: TORSION __FILL__
r9-phi: TORSION __FILL__
r9-psi: TORSION __FILL__
r10-phi: TORSION __FILL__
r10-psi: TORSION __FILL__
r11-phi: TORSION __FILL__
r11-psi: TORSION __FILL__
r12-phi: TORSION __FILL__
r12-psi: TORSION __FILL__
r13-phi: TORSION __FILL__
r13-psi: TORSION __FILL__
r14-phi: TORSION __FILL__
r14-psi: TORSION __FILL__
r15-phi: TORSION __FILL__
r15-psi: TORSION __FILL__
r16-phi: TORSION __FILL__
r16-psi: TORSION __FILL__

# This command stores all the Ramachandran angles that were computed
angles: COLLECT_FRAMES __FILL__=r2-phi,r2-psi,r3-phi,r3-psi,r4-phi,r4-psi,r5-phi,r5-psi,r6-phi,r6-psi,r7-phi,r7-psi,r8-phi,r8-psi,r9-phi,r9-psi,r10-phi,r10-psi,r11-phi,r11-psi,r12-phi,r12-psi,r13-phi,r13-psi,r14-phi,r14-psi,r15-phi,r15-psi,r16-phi,r16-psi
# Lets now compute the matrix of distances between the frames in the space of the Ramachandran angles
distmat: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=__FILL__ METRIC=EUCLIDEAN
# Now select 500 landmark points to analyze
fps: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=__FILL__ NLANDMARKS=500
# Run sketch-map on the landmarks
smap: SKETCH_MAP __FILL__=fps NLOW_DIM=2 HIGH_DIM_FUNCTION={SMAP R_0=6 A=8 B=2} LOW_DIM_FUNCTION={SMAP R_0=6 A=2 B=2} CGTOL=1E-3 CGRID_SIZE=20 FGRID_SIZE=200 ANNEAL_STEPS=0
# Project the remaining trajectory data
osample: PROJECT_ALL_ANALYSIS_DATA USE_OUTPUT_DATA_FROM=__FILL__ PROJECTION=__FILL__

# This command outputs all the projections of all the points in the low dimensional space
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=__FILL__ ARG=osample.* FILE=smap_data

# These next three commands calculate the secondary structure variables.  These
# variables measure how much of the structure resembles an alpha helix, an antiparallel beta-sheet
# and a parallel beta-sheet.  Configurations that have different secondary structures should be projected
# in different parts of the low dimensional space.
alpha: ALPHARMSD RESIDUES=all
abeta: ANTIBETARMSD RESIDUES=all STRANDS_CUTOFF=1.0
pbeta: PARABETARMSD RESIDUES=all STRANDS_CUTOFF=1.0

# These commands collect and output the secondary structure variables so that we can use this information to
# determine how good our projection of the trajectory data is.
cc2: COLLECT_FRAMES ARG=alpha,abeta,pbeta
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=cc2 ARG=cc2.* FILE=secondary_structure_data
```

This input collects all the torsional angles for the configurations in the trajectory.  Then, at the end of the calculation, the matrix of distances between these points is computed, and a set of landmark points
is selected using a method known as farthest point sampling.  A matrix that contains only those distances between the landmarks is then constructed and diagonalized by the CLASSICAL_MDS action, and this
set of projections is used as the initial configuration for the various minimization algorithms that are then used to optimize the sketch-map stress function.  As in the previous exercise, once the projections of
the landmarks are found, the projections for the remainder of the points in the trajectory are found by using the PROJECT_ALL_ANALYSIS_DATA action.  Try to fill in the blanks in the input above
and run this calculation now using the command:

````
plumed driver --mf_pdb traj.pdb
````

Once the calculation has completed, you can, once again, visualize the data generated using chemiscope by using the scripts from earlier exercises.  My projection is shown in the figure below.  Points are, once again,
coloured following the alpha-helical content of the corresponding structure.

![A representation of the sketch-map projection that was generated using chemiscope.](figures/masterclass-21-6-chemiscope-smap.png)

__Try to see if you can reproduce an image that looks like the one above__

#### Summary

This exercise has shown you that running dimensionality reduction algorithms using PLUMED involves the following stages:

- Data is collected from the trajectory using [COLLECT_FRAMES](https://www.plumed.org/doc-master/user-doc/html/_c_o_l_l_e_c_t__f_r_a_m_e_s.html).
- Landmark points are selected using a [landmarks algorithm](https://www.plumed.org/doc-master/user-doc/html/_analysis.html#landmarks)
- The distances between the trajectory frames are computed using [EUCLIDEAN_DISSIMILARITIES](https://www.plumed.org/doc-master/user-doc/html/_e_u_c_l_i_d_e_a_n__d_i_s_s_i_m_i_l_a_r_i_t_i_e_s.html)
- A loss function is optimized to generate projections of the landmarks.
- Projections of the non-landmark points are found using [PROJECT_ALL_ANALYSIS_DATA](https://www.plumed.org/doc-master/user-doc/html/_p_r_o_j_e_c_t__a_l_l__a_n_a_l_y_s_i_s__d_a_t_a.html).

There are multiple choices to be made in each of the various stages described above.  For example, you can change the particular sort of data this is collected from the
trajectory. There are multiple different ways to select landmarks. You can use the distances directly, or you can transform them. You can use various loss functions and
optimize the loss function using a variety of different algorithms.  When you tackle your own research problems using these methods, you can experiment with the various choices that can
be made.  

__Generating a low-dimensional representation using these algorithms is not enough to publish a paper.__  We want to get some physical understanding from our simulation. Gaining this understanding is the hard part.

### Exercise 4: Path collective variables

In many papers in this area, you will hear people talk about the distinction between collective variables and the reaction coordinate.  The distinction these authors are trying emphasize when they use this language 
is between a descriptor that takes different values for the various important basins in the energy landscape (the collective variable) and the actual pathway the reaction proceeds along (the reaction coordinate).
Furthermore, the authors of such papers will often argue that the reaction coordinate is simply some linear/non-linear combination of collective variables.
In this exercise, we will study alanine dipeptide with path collective variables to show you one way of thinking about this distinction between collective variables and reaction coordinates.  
By studying a system that is this simple, you will also hopefully come to understand how we can interpret the coordinates that we extract using the dimensionality reduction algorithms that were discussed in the previous exercise.

To remind you, the free energy surface as a function of the $\phi$ and $\psi$ angles for alanine dipeptide is shown below:

![The Free energy landscape of alanine dipeptide in Ramachandran angles in the CHARMM27 force field](belfast-2-rama.png)

In other masterclasses, we have discussed how there are two critical states in the above energy landscape.  These states are labelled C7eq and C7ax above.  
The two Ramachandran angles plotted on the x and y axes of the free energy surface
above are examples of what we have called collective variables.  Both of these angles can be used to distinguish between C7eq and C7ax configurations.
The reaction coordinate is the path that the system actually takes as it moves from the C7ax to the C7eq configuration.  Based on the shape of the free
energy surface, we might suppose that this reaction  coordinate looks something like the black line shown below: 

![Paths that connect the C7eq and C7ax configuration](belfast-2-good-bad-path.png) 

The file called alanine-transformation.pdb that you can find in the data directory of the GitHub repository contains a series of configurations that lie close to the transition path that is illustrated in black in the figure above.  Below are plots that show how 
$\phi$ and $\psi$ change as the system moves along this path.  __Try to see if you can use what you have learned in previous masterclasses to reproduce the figure above before continuing.__  

![Changes in the Ramachandran angles as the protein transitions from the C7ax to C7eq state](figures/masterclass-21-6-rama-transition.png)

#### RMSD distances

We know what structures correspond to the various stable states of our system.  We might, therefore, be tempted to ask ourselves if we can not just use the RMSD distance from a structure as the reaction coordinate.
This approach will work if there is a single important configuration in the energy landscape.  We could use the RMSD distance from this lowest energy structure as a CV to extract this configuration's free energy
relative to everything else.  How well does this work if we have two distinct states, though, as we have for alanine dipeptide?  To investigate this question, fill in the PLUMED input below that calculates the RMSD distances 
from the C7ax and C7eq configurations:

```plumed
#SOLUTIONFILE=work/plumed_ex5.dat
c7ax: RMSD __FILL__=../data/C7ax.pdb TYPE=__FILL__ 
c7eq: RMSD __FILL__=../data/C7eq.pdb TYPE=__FILL__
PRINT ARG=c7ax,c7eq FILE=colvar STRIDE=100
```

You can run an MD simulation to monitor these distances using the python script below.

```python
def generate_gromacs_input( directory, plumed_input, nsteps, temp ) :
    # Setup the initial configuration by picking a frame at random from the startpoints file
    allc = io.read('../data/startpoints.pdb',':')
    nframes = len(allc)
    io.write( directory + '/conf.pdb', allc[int(np.random.randint(0,nframes))] )
    # Copy the topology to the appropriate directory
    shutil.copyfile("../data/topol.top", directory + "/topol.top")
    # Setup the mdp file
    mdp = open("../data/md.mdp","r")
    contents = mdp.read()
    mdp.close()
    new_content = contents.replace("SEED", str(np.random.randint(0,1000000))).replace("NSTEPS",str(nsteps)).replace("TEMP",str(temp))
    mdpout = open(directory + "/md.mdp", "w")
    mdpout.write(new_content)
    mdpout.close()
    # And write the plumed input
    pout = open(directory + "/plumed.dat", "w")
    pout.write( plumed_input )
    pout.close()
    return

plm = '''
INSERT PLUMED INNPUT HERE
'''
# Create a directory to run the calculation in
!rm -rf ../Unbiased_MD && mkdir ../Unbiased_MD
# Generate gromacs input for MD simulation at 1000 K
generate_gromacs_input( '../Unbiased_MD', plm, 500000, 1000 )
# Now run gromacs
!cd ../Unbiased_MD/ && gmx_mpi grompp -f md.mdp -c conf.pdb -p topol.top -maxwarn 2 &> /dev/null
!cd ../Unbiased_MD/ && gmx_mpi mdrun --plumed plumed.dat &> /dev/null
# And read in the colvar files
unbiased_data = np.loadtxt('../Unbiased_MD/colvar')
```

I ran simulations at 300K and 1000K using the above script.  When the simulations completed, I was able to construct the figures below:

![The black points are the configurations visited in MD simulations at 300K (left) and 1000K (right).  The top panel shows these trajectories in the Ramachandran plane.  In the lower panel, I have projected the trajectories on the RMSD distances from the two reference states.](figures/masterclass-21-6-rmsd-distances.png)

The black points above are the RMSD values for the trajectories.  I have also shown the RMSD values for the frames in
alanine-transformation in red.  Notice that all the configurations explored in the MD are very distant from the C7ax and C7eq states.
__Try now to reproduce the above figure yourself.__

You can see that at 300K, the simulation we ran did not visit the C7eq state.  Furthermore, you can also clearly see this from the 
projections of the trajectories on the RMSD distances from these two configurations.  Notice in these figures, however, that the distances
from these two reference configurations were often considerably larger than the distance between these reference configurations and the configurations on the reference path in alanine-transformation.pdb.

#### Using dot products to calculate CVs

Instead of calculating two distances, we might ask ourselves if the linear combination of $\phi$ and $\psi$ that is illustrated in the figure below:

![An illustration showing how PCAVARS coordinates work.  The vector connecting some reference state to any state the system is in can be in (purple and orange points) can be projected onto the vector connecting the two states of interest (black arrow) by using the dot product of the vectors shown here.](marvel-2-pca-coordinates.png)

gives a better description for the transition.  We can define this CV as follows:

$$
s = \frac{(\phi_2 - \phi_1).(\phi_3 - \phi_1) + (\psi_2 - \psi_1).(\psi_3 - \psi_1)}{ \sqrt{(\phi_2 - \phi_1)^2 + (\psi_2 - \psi_1)^2} }
$$

In this expression $(\phi_1,\psi_1)$ are the Ramachandran angles in the $C_7eq$ configuration and 
$(\phi_2,\psi_2)$ are the Ramachandran angles in the $C_7ax$.  $(\phi_3,\psi_3)$ is the 
instantaneous configuration.  You should thus be able to see that we have arrived at the above expression
by using our knowledge of the dot product between two vectors.

__See if you can write an input to re-analyse the data in alanine-transformation.pdb and the MD simulations from the previous section using this CV.__  
You should be able to get plots of the value of this CV as a function of step number that looks like the ones shown below:

![The black points are the configurations visited in MD simulations at 300K (left) and 1000K (right).  The top panel shows these trajectories in the Ramachandran plane.  In the lower panel, I have projected the trajectories on vector in the Ramachandran plane shown above.](figures/masterclass-21-6-pcavars-transition.png)

I implemented this CV using a combination of TORSION and COMBINE.  I also ignored the fact that the torsions are periodic variables when calculating the linear combination.  
You can't use the sort of linear algebra above with periodic variables, but it is OK for these illustrative purposes.  

Notice that the first frame in alanine-transformation.pdb has the molecule in the $C_7ax$ configuration.  The last frame has the molecule in the $C_7eq$ state.  

#### PCA Collective Variables

In a previous section we saw how we can construct a new collective variables by taking a linear combination of two other 
variables.  This idea can be extended to higher dimensions, however.  As long as we can find the vector that connectes the $C_7eq$ and 
$C_7ax$ states we can project our current coordinates on that particular vector.  We can even define this vector in the space of the 
coordinates of the atoms.  In other words, if the 3$N$ coordinate of atomic positions is $\mathbf{x}^{(1)}$ for the $C_7eq$ 
configuration and $\mathbf{x}^{(2)}$ for the $C_7ax$ configuration and if the instantaneous configuration of the atoms is $\mathbf{x}^{(3)}$ we 
can use the following as a CV:

$$
s = \frac{\sum_{i=1}^{3N} (x^{(2)}_i - x^{(1)}_i ) (x^{(3)}_i - x^{(1)}_i )}{ \sqrt{\sum_{i=1}^{3N} (x^{(2)}_i - x^{(1)}_i )^2} } 
$$

as long as rotational and translational movement is removed before calculating the two displacement vectors above. 

We can also get some sense of how far the system has moved from this vector by computing:

$$
z = \sqrt{ \sum_{i=1}^{3N} (x^{(3)}_i - x^{(1)}_i)^2 - s^2 }
$$

which follows by applying Pythagoras theorem.

You can calculate this collective variable using PLUMED by using the input below:

```plumed
#SOLUTIONFILE=work/plumed_ex6.dat
p: PCAVARS __FILL__=pca-reference.pdb __FILL__=OPTIMAL
PRINT ARG=p.* FILE=colvar
```

__Use this input to re-analyse the data in `alanine-transformation.pdb` and your MD simulations before continuing.__  The figure below shows the results that I obtained by analyzing the 
data using the PCAVARS command above.

![The black points are the configurations visited in MD simulations at 300K (left) and 1000K (right).  The top panel shows these trajectories in the Ramachandran plane.  The projection on the vector connecting the C7ax and C7eq states (left) and the perpendicular distance from this vector (right) and shown in the bottom panels.](figures/masterclass-21-6-pcavars.png)

The figure above shows that this coordinate is good.  We move quite far from this initial vector when we move to the C7ax state, however.

#### Two PCA Collective variables

Instead of using a single PCA variable and calculating the projection on the vector connecting these two points, we can instead use the projection of the instantaneous coordinates 
in the plane that contains three reference configurations.  You can calculate these collective variables using PLUMED by using the input below:

```plumed
#SOLUTIONFILE=work/plumed_ex7.dat
p: PCAVARS __FILL__=pca2-reference.pdb __FILL__=OPTIMAL
PRINT ARG=p.* FILE=colvar
```

Notice that there are now three structures within `pca2-reference.pdb`, and thus two vectors are defined.  __Use this input to re-analyse the data in `alanine-transformation.pdb` and your MD simulations before continuing.__  
The figure below shows the results that I obtained by analyzing the
data using the PCAVARS command above.

![The black points are the configurations visited in MD simulations at 300K (left) and 1000K (right).  The top panel shows these trajectories in the Ramachandran plane.   The projection on the vectors defined using the configurations in pca2-reference.pdb is shown in the bottom panel.  The green points in the top panel give the Ramachandran angles for the configurations used to define the reference configurations.](figures/masterclass-21-6-pcavars2.png)

#### Path collective variables

Instead of using the projection on more than one vector, we can extend the ideas in the previous section and use curvilinear paths rather than linear paths.  This trick is what is done with the PATH CVs that are illustrated in the figure below:

![The S variable can be thought as the length of the red segment, while the Z variable is the length of the green one](belfast-2-ab-sz.png)

As you can see, there are two path collective variables, $S(X)$ measures your position on a path and is calculated using:

$$
S(X)=\frac{\sum_{i=1}^{N} i\ \exp^{-\lambda \vert X-X_i \vert }}{ \sum_{i=1}^{N} \exp^{-\lambda \vert X-X_i \vert } }
$$

$Z(x)$, meanwhile, measures the distance from the path using:

$$
Z(X) = - \frac{1}{\lambda} \ln\left[ \sum_{i=1}^{N} \exp^{-\lambda \vert X-X_i \vert } \right]
$$

We can calculate these CVs using a filled-in version of the input that is shown below:

```plumed
#SOLUTIONFILE=work/plumed_ex8.dat
path: PATH __FILL__=../data/alanine-path.pdb TYPE=__FILL__ LAMBDA=15100.
PRINT ARG=* FILE=colvar
```

The figure below shows the $\phi$ and $\psi$ angles for the configurations that define the path as red dots.  Meanwhile, the black dots are the 
$\phi$ and $\psi$ angles that are visited during a high temperature, unbiased MD trajectory.  The green dots are then the configurations in alanine-transformation.pdb 
You can clearly see that the red dots lie along the transition path that connects the C7eq and C7ax configurations and that the green dots follow this pathway: 

![The path that we have defined for our transition in the Ramachandran plane.  Black dots are the configurations visited at 300 K (left) and 1000 K (right)](figures/masterclass-21-6-path-angles.png)

When we calculate the path collective variables for the black and red data points in the above figure, we get the result shown on the right. The panel on the left shows the output from 
a 300 K simulation during which the C7eq configuration was not visited:

![An analysis of the path CV values visited during trajectories at 300 K and 1000 K.](figures/masterclass-21-6-path-path.png)

You can see that the majority of configurations have very small z values.  It would thus seem that we have identified our reaction coordinate.  We can use the S-path coordinate
to tell us where we are on the transition pathway between the two states.  Furthermore, in defining this coordinate, we have used the positions of all the heavy atoms in the protein.  

#### The isocommittor

Notice that when we talk about the reaction coordinate, we are not talking about a single set of configurations that connect the two states.
The previous sections have shown that there are always multiple pathways that connect the two states.  The task of identifying a single reaction coordinate thus appears
impossible.  How, after all, can there be a single reaction path if we know that there are multiple pathways that connect the two states?  

One way of answering this question is to look at how far the transitions you observe deviate from the reference transition path using the Z-coordinate.  An
alternative answer can be found by considering so-called isocommittor surfaces.  When using this method, you suppose there is a
saddle point between the two states of interest (the $C_7ax$ and $C_7eq$ configurations in our alanine dipeptide example).  

Let's suppose that we now start a simulation with the system balanced precariously on this saddle point.  The system will, very-rapidly, fall off the maximum and move
towards either the left or the right basin.  If we are exactly at the saddle, 50% of the trajectories will fall to the right, and 50% will fall to the left. 

We can use the idea of the isocommittor to identify whether any collective variable is giving a good description of the transition state ensemble's location.
The idea in such simulations is to define regions of space that correspond to the various stable states in the free energy landscape.  We then launch lots of simulations
from configurations that are not within these basins and identify the basin that these calculations visit first.  We can then make graphs like those shown below. These graphs 
show the fraction of times a particular CV value was visited in a trajectory that visited the c7eq basin before it visited the c7ax basin. 

![Probabilities of visiting the C7eq basin before the C7ax basin if you are differents points on each of the CV ranges](figures/masterclass-21-6-committor.png)

The PLUMED input that I used when making the figure above is a filled in version of the following:

```plumed 
#SOLUTIONFILE=work/plumed_ex9.dat
phi: TORSION __FILL__=5,7,9,15
psi: TORSION __FILL__=7,9,15,17
pca: PCAVARS __FILL__=../data/pca-reference.pdb TYPE=__FILL__
PRINT __FILL__=phi,psi,pca.eig-1,path.spath FILE=colvar STRIDE=1
__FILL__ ARG=phi STRIDE=1 BASIN_LL1=-3 BASIN_UL1=-1 BASIN_LL2=1 BASIN_UL2=2 FILE=basin
```

I deployed the input above within the following python script that launches a large number of gromacs simulations:

```python
def gen_trajectories( ntraj, plm ) : 
    nc7eq, c7eq_colv, all_data = 0, np.zeros([0,5]), np.zeros([0,5])
    for i in range(ntraj) :
        !rm -rf ../Test && mkdir ../Test
        generate_gromacs_input( '../Test', plm, 10000000, 300 )
        !cd ../Test/ && gmx_mpi grompp -f md.mdp -c conf.pdb -p topol.top -maxwarn 2 &> /dev/null 
        !cd ../Test/ && gmx_mpi mdrun --plumed plumed.dat &> /dev/null
        bfile = open('../Test/basin','r')
        with open( '../Test/basin', "r" ) as myfile :
            for line in myfile :
                if line.startswith("#! SET COMMITTED TO BASIN") : basin = line.split()[5]
        colv_data = np.loadtxt("../Test/colvar")
        if len(colv_data.shape)==1 : colv_data = np.reshape( colv_data, (1,5) )
        all_data = np.concatenate( (all_data, colv_data), axis=0 )
        if basin=="1" : 
            c7eq_colv = np.concatenate( (c7eq_colv, colv_data), axis=0 )
            nc7eq = nc7eq + 1
    print( "NUMBER c7ax", ntraj-nc7eq, "c7eq", nc7eq )
    return c7eq_colv, all_data

p='''
YOUR PLUMED INPUT GOES HERE
'''
bas_data, all_data = gen_trajectories( ntraj, p )
```

Notice how the above script stores the CV values from trajectories that finished in the C7eq basin in the NumPy array called `bas_data`. 
We can use the data in this file to construct the (unnormalized) histograms of CV value that finished in C7eq and the total histogram. 
The script I used to construct these histograms is shown below.  This script also shows how we can arrive at the final conditional probability distributions from 
the figure above by dividing by the total number of configurations that moves to C7eq first by the total number of configuration that visited 
a point.

```python
def histo( data, nbins, xmin, xmax ) :
    delr = ( xmax - xmin ) / nbins 
    hist = np.zeros(nbins)
    for d in data : 
        xbin = int( np.floor( (d-xmin) / delr ) )
        hist[xbin] = hist[xbin] + 1
    return hist / delr 

def get_isocommitor( bas_data, full_data, nbins, xmin, xmax ) :
    bas_histo = histo( bas_data, nbins, xmin, xmax ) 
    full_histo = histo( full_data, nbins, xmin, xmax )
    for i in range(nbins) : 
        if np.abs(full_histo[i])<1E-4 : bas_histo[i] = 0 
        else : bas_histo[i] = bas_histo[i] / full_histo[i]
    return bas_histo 

# This makes the isocommittor as a function of \f$\phi\f$:
commit = get_isocommitor( bas_data[:,1],  all_data[:,1], 30, -np.pi, np.pi )
```

Notice that the script above does not compute the error bars I included in the figure.  __Your task in this exercise is to reproduce this 
figure with the error bars.__   I generated the error bars in my figure by running 10 sets of 500 simulations.  I constructed separate histograms
from each of these batches of simulations and calculated the mean and variance from these ten sets of samples. 

#### Conclusions

The exercises in this section aimed to discuss the difference between collective variables and the reaction coordinate.  I have argued that collective variables are developed by researchers thinking about coordinates that can distinguish between the important states.  Reaction coordinates, meanwhile, are the actual pathways
that reactions proceed along.  These reaction coordinates are found by probing the data we get from simulations.    

I am not convinced that the distinction between reaction coordinates and collective variables is very important.  Every coordinate we have tried in the previous sections has allowed us to determine whether we are in the C7eq or the C7ax basin.   My view is that when we do science, we impose structure on our results by asking questions.  Our question has been
to determine the relative free energies of the C7ax and C7eq basins in this exercise.  To do this, we have to define what configurations are in the C7ax basin and what configurations are in the C7eq basin.  The distinction between these two structures is an arbitrary figment of our imagination.  If we use different variables, we might be defining these states differently, and we might get slightly different answers. 
We should not be surprised by the fact that we get different answers if we ask the question differently.  The important thing is to develop an interesting question and to tell a good story about whatever you find out. 

### Exercise 4: Dealing with indistinguishability

Many researchers have used these rare event methods to study nucleation and crystal growth.  Studying such problems introduces an additional challenge when designing collective variables as all the
atoms are indistinguishable.  You thus know before even running the simulation that there are multiple paths that connect the reactant (liquid) state and the product (solid) state.  The nucleation, 
after all, can start with any atom in the simulation box.  The first step in developing CVs for studying nucleation processes is thus to identify some order parameter that allows us to distinguish atoms that 
are in the solid state from atoms that are in the liquid state.  The following input, once it is filled in will calculate some suitable order parameters for all the atoms in the system: 

```plumed
#SOLUTIONFILE=work/plumed_ex10.dat
UNITS NATURAL
coord: COORDINATIONNUMBER __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} 
cub: FCCUBIC __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
fcub: MTRANSFORM_MORE DATA=__FILL__ SWITCH={SMAP R_0=0.45 D_0=0.0 A=8 B=8}
DUMPMULTICOLVAR __FILL__=coord STRIDE=1000 __FILL__=coord.xyz
DUMPMULTICOLVAR __FILL__=cub STRIDE=1000 __FILL__=cub.xyz 
DUMPMULTICOLVAR __FILL__=fcub STRIDE=1000 __FILL__=fcub.xyz
```

You can run a short MD simulation on Lennard-Jonesium to compute those order parameter by using the following input (in) file for simplemd:

````
inputfile interface.xyz
outputfile output.xyz
temperature 0.6
tstep 0.005
friction 1
forcecutoff 2.5
listcutoff  3.0
nstep 10000
nconfig 100 trajectory.xyz
nstat   100 energies.dat
````

You can run this from a python notebook and read in the values of the order parameters by using the following script: 

```python
p = '''
INSERT YOUR PLUMED INPUT HERE
'''

inp = '''
INSERT YOUR SIMPLEMD INPUT HERE
'''

# Make a directory to run in 
!rm -rf ../LJ-trajectories/First-run && mkdir ../LJ-trajectories/First-run
# Copy the initial configuration
!cp ../data/interface.xyz ../LJ-trajectories/First-run
# Output the plumed file
f = open("../LJ-trajectories/First-run/plumed.dat", 'w')
f.write(p)
f.close()
# Output the input to simplemd
f = open("../LJ-trajectories/First-run/in", 'w')
f.write(inp)
f.close()

# Now run PLUMED
!cd ../LJ-trajectories/First-run/ && plumed simplemd < in &> /dev/null

# Read in the various order parameters (N.B. you can almost certainly write better python here)
!grep X ../LJ-trajectories/First-run/coord.xyz | awk '{print $5}' > ../LJ-trajectories/First-run/coord.dat
coord = np.loadtxt("../LJ-trajectories/First-run/coord.dat")
!grep X ../LJ-trajectories/First-run/cub.xyz | awk '{print $5}' > ../LJ-trajectories/First-run/cub.dat
cub = np.loadtxt("../LJ-trajectories/First-run/cub.dat")
!grep X ../LJ-trajectories/First-run/fcub.xyz | awk '{print $6}' > ../LJ-trajectories/First-run/fcub.dat
fcub = np.loadtxt("../LJ-trajectories/First-run/fcub.dat")
```

A visualization of the order parameters for all the atoms can then be produced using chemiscope, which allows you to see the value of all the individual atom order parameters.  

```python
import ase
import ase.io
from chemiscope import write_input

# Read in the trajectory
traj = ase.io.read('../LJ-trajectories/First-run/trajectory.xyz',':')

# This constructs the dictionary of properties for chemiscope
properties = {
    "coord": {
        "target": "atom",
        "values": coord,
        "description": "Coordination number of atom",
    },
    "cub": {
        "target": "atom",
        "values": cub,
        "description": "FCCCUBIC order parameter",
    },
    "fcub": {
        "target": "atom",
        "values": fcub,
        "description": "Transformed FCCUBIC order parameter",
    },
}

# This generates our chemiscope output
write_input("fccubic_chemiscope.json.gz", cutoff=1.5, frames=traj, properties=properties )
```

__Put all the above together and see if you can produce a chemiscope representation like the one below__

![A representation of a solid-liquid system that was generated using chemiscope](figures/masterclass-21-lj-chemiscope.png)

Chemiscope is invaluable for determining whether our order parameter is good at distinguishing the atoms within the solid from those within the liquid.   If you visualize the FCCCUBIC or the transformed values of these parameters in the example above, you see that roughly half of the atoms are solid-like, and nearly half of the atoms are liquid-like.  Notice also that we can use the dimensionality reduction that were discussed in earlier parts of this masterclass
and clustering algorithms to develop these atomic order parameters.  Generating the data to analyze using these algorithms is easier because we can get multiple sets of coordinates to analyze from each trajectory frame.  The atoms are, 
after all, indistinguishable.

Once we have identified an order parameter we can analyse the distribution of values that it takes by using an input like the one shown below:

```plumed
#SOLUTIONFILE=work/plumed_ex11.dat
UNITS NATURAL
coord: COORDINATIONNUMBER __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} 
cub: FCCUBIC __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
fcub: MTRANSFORM_MORE DATA=__FILL__ SWITCH={SMAP R_0=0.45 D_0=0.0 A=8 B=8}
coord_histo: HISTOGRAM __FILL__=coord STRIDE=10 __FILL__=2 GRID_MAX=15 __FILL__=100 __FILL__=DISCRETE
cub_histo: HISTOGRAM __FILL__=cub STRIDE=10 GRID_MIN=-1 __FILL__=1 __FILL__=100 __FILL__=DISCRETE
fcub_histo: HISTOGRAM __FILL__=fcub STRIDE=10 __FILL__=0 GRID_MAX=1 __FILL__=100 __FILL__=DISCRETE
DUMPGRID __FILL__=coord_histo FILE=coord_histo.dat
DUMPGRID __FILL__=cub_histo FILE=cub_histo.dat
DUMPGRID __FILL__=fcub_histo FILE=fcub_histo.dat
```

The resulting distributions of order parameter values should look like this:

![The distributions of the various order parameters in a simulation of the solid-liquid system.](figures/masterclass-21-6-order-params.png)

__Try to use the information above to reproduce this plot.__ Notice that the distribution of order parameters for the FCC Cubic order parameters is bimodal.  There is a peak corresponding to the value that this quantity takes in the box's liquid part. The second peak then corresponds to the value that the quantity takes in the box's solid part.  The same cannot be said for the coordination number by contrast.  Notice, last of all, that by transforming the values above
we have an order parameter that is one for atoms that are within the solid and zero otherwise. 

It is tempting to argue that the number of solid particles is equal to the number of atoms that have an order parameter that is greater than some threshold.  A better way to calculate the number of solid
particles, however, is to proceed as follows:

- Run simulations of the solid and liquid under the same thermodynamic conditions.
- Calculate the average values per atom value of the order parameter in these simulations of the solid $\phi_s$ and liquid $\phi_l$.
- Calculate the number of solid \f$n_l\f$ and liquid atoms \f$n_l\f$ by solving the following pair of simultaneous equations $N=n_s + n_l$ and $\Phi = n_s \phi_s + n_l \phi_l$, where $N$ is the total number of atoms. $\Phi$, meanwhile, is given by:

$$
\Phi = \sum_{i=1}^N \phi_i
$$

The sum runs over all the atoms and where $\phi_i$ is the order parameter for the $i$th atom in the system.

This procedure works because $\Phi$ is an extensive quantity and is thus additive.  If we have a system that contains a aolid liquid interface we can express the value of any extensive quantity, $\Psi$ using:

$$
\Psi = n_s \psi_s + n_l \psi_l + \psi_i
$$

where $\psi_s$ and $\psi_l$ is the average value per atom value of the quantity $\psi$ in the solid and liquid phases.  $\psi_i$, meanwhile, is the so-called surface excess term which measures the contribution that the presence of the interface makes to the value of $\Psi$.  When we find $n_s$ by solving the simultaneous equations above, we assume that this excess term is zero.

__To complete this exercise, you should run simulations of the solid, liquid and interface.  You should use the ensemble average of the order parameter for your solid and liquid simulations to make a graph like the one below that shows the number of atoms of solid as a function of time for the system that contains the interface__.  The number of solid atoms has been determined by solving the simultaneous equations using the theory above.  Notice that there are input configurations for the solid, liquid and interface systems in the data directory of the GitHub repository.  The solid configuration is in `solid.xyz`. The liquid is in `liquid.xyz`, and the interface is in `interface.xyz`.  The final result you get will look something like this:

![The number of solid atoms calculated by assuming a zero excess term for the sum of coordination numbers (black), the sum of fcc order parameters (red) and the sum of transformed fcc order parameters (green)](figures/masterclass-21-6-solid-liquid.png)

Notice that the green line is the least noisy of these curves.  The reason the line is less noisy is connected to the fact that you have the clearest delineation between solid and liquid atoms when you use this order parameter.

NB. We are running at constant volume because we are keeping things simple and using simplemd.  If you are studying nucleation using these techniques, you should run at constant pressure. 

### Applying these ideas in your research

The exercises in this tutorial have introduced you to the following ideas:

- Using chemiscope to visualize trajectories and collective variables.
- Using dimensionality reduction algorithms
- Using Path collective variables
- Clustering trajectories
- Dealing with indistinguishable atoms

What I would like to do in the follow-up session is to think about how you can use these ideas as you do your research.  The first step in using any of these questions is asking a 
research question that is amenable to answering using one of these methods.  In this final exercise, I would like you to think about using these methods in your research.
The first step in doing so is asking yourself the question __what am I trying to find out?__  The clearer you make this question, the easier you will find the subsequent research project.
With this in mind, here are some examples of terrible research questions:

- Understanding protein-ligand interactions using machine learning algorithms.
- Understanding the mechanism via which limescale forms in pipes
- Understanding protein folding using machine learning.
- Studying transport through membrane proteins.
- Understanding the catalytic mechanism of iron in ammonia formation.

These research questions are terrible because the final result that you will obtain has not been specified.  It is unclear when the project will be finished, or what success in this project looks like.  With this in mind, consider the following alternative questions, which each have a clear goal.  Each of the goals below maps on to one of the terrible research questions above:

- Calculate the free energy change when ligand A binds to protein B.
- Calculate the free energy change when single calcium and carbonate ions bind to copper surfaces.
- Calculate the free energy change when protein A transitions to its folded state.
- Calculate the free energy change as a sodium ion moves through a sodium channel protein.
- Calculate the free energy change when nitrogen, hydrogen and ammonia gas bind to iron surfaces.

For each of the questions above, it is clear what success will look like -- you will have a converged estimate of the free energy difference.  We have shown you how such estimates can be extracted in previous masterclasses and what it looks like when such calculations go wrong.  Answering these research questions is thus simply a matter of applying what you have learned. 
Notice, also, how for each of the questions above, it is clear what you can do in the first step:

- Run a parallel tempering metadynamics simulation with the distance between the centres of mass of the protein and the ligand as a CV.
- Run a metadynamics simulation of a slab of copper in contact with an aqueous solution containing a single calcium ion.  Use the z component of the vector connecting the calcium ion and the centre of mass of the copper slab as the CV.
Run a metadynamics simulation using the distance from the protein's folded state as the collective variable.
- Run a metadynamics simulation of the sodium in the channel protein.  Use the z position of the sodium relative to the centre of the channel protein as a CVs. 
- Run a metadynamics simulation of a slab of iron in contact with a simulation box that is empty except for one nitrogen molecule.  Use the z component of the vector connecting the centre of mass of the nitrogen and the centre of mass of the iron slab as the CV.

These simulations may or may not give you the desired free energy difference.  You will at least get some trajectory data from them, however.  You can then analyze the trajectories using a machine learning
algorithm or similar to refine your approach.  This refining stage is where the experience of your supervisor should be invaluable.

With all this in mind, think about your research project.  Try to develop a research question to work on and an approach that you can use to tackle this question.  Your question should not be some grand challenge.  
It should be something simple that you know how to do.  You do not need to run these simulations.  I want you to think about the research question to discuss them at the follow-up session.  I want you to think about how you 
can use these methods in your research as if you are not doing so, attending these masterclasses is pointless as you will never really apply these methods.  
