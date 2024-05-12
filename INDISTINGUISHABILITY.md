# PLUMED Masterclass 21.6: Dimensionality reduction

## Aims

The primary aim of this Masterclass is to show you how you might use PLUMED in your work.
We will see how to call PLUMED from a python notebook and discuss some strategies for selecting
collective variables.

## Objectives

Once this series of exercise is completed, users will be able to:

- Use Gibbs dividing surfaces to determine the extent of a phase. 

## Resources

The data needed to complete this Masterclass can be found on [GitHub](https://github.com/plumed/masterclass-21-6).
You can clone this repository locally on your machine using the following command:

````
git clone https://github.com/plumed/masterclass-21-6.git
````

I recommend that you run each exercise in a separate sub-directory (i.e. Exercise-1, Exercise-2, ...), which you can create inside the root directory `masterclass-21-6`. Organizing your data this way will help you to keep things clean.

_All the exercises have been tested with PLUMED version 2.7.0._

## Acknowledgements

Throughout this exercise, we use the [atomistic simulation environment](https://wiki.fysik.dtu.dk/ase/)
and [chemiscope](https://chemiscope.org/).  Please look at the information at the links I have provided here
for more information about these codes.

## Exercises

Many researchers have used these rare event methods to study nucleation and crystal growth.  Studying such problems introduces an additional challenge when designing collective variables as all the
atoms are indistinguishable.  You thus know before even running the simulation that there are multiple paths that connect the reactant (liquid) state and the product (solid) state.  The nucleation, 
after all, can start with any atom in the simulation box.  The first step in developing CVs for studying nucleation processes is thus to identify some order parameter that allows us to distinguish atoms that 
are in the solid state from atoms that are in the liquid state.  The following input, once it is filled in will calculate some suitable order parameters for all the atoms in the system: 

```plumed
#SOLUTIONFILE=work/plumed_ex10.dat
UNITS NATURAL
coord: COORDINATIONNUMBER __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} 
cub: FCCUBIC __FILL__=1-5184 __FILL__={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
fcub: MORE_THAN ARG=__FILL__ SWITCH={SMAP R_0=0.45 D_0=0.0 A=8 B=8}
fcub_grp: GROUP ATOMS=cub
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
fcub: MORE_THAN ARG=__FILL__ SWITCH={SMAP R_0=0.45 D_0=0.0 A=8 B=8}
fcub_grp: GROUP ATOMS=cub
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
- Calculate the number of solid $n_l$ and liquid atoms $n_l$ by solving the following pair of simultaneous equations $N=n_s + n_l$ and $\Phi = n_s \phi_s + n_l \phi_l$, where $N$ is the total number of atoms. $\Phi$, meanwhile, is given by:

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
