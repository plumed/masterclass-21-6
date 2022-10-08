# PLUMED Masterclass 21.6: Dimensionality reduction

## Aims

The primary aim of this Masterclass is to show you how you might use PLUMED in your work.
We will see how to call PLUMED from a python notebook and discuss some strategies for selecting
collective variables.

### Applying these ideas in your research

The exercises in this lesson have introduced you to the following ideas:

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

With all this in mind, think about your research project.  Try to develop a research question to work on and an approach that you can use to tackle this question.  Your question should not be some grand challenge.  It should be something simple that you know how to do.  You do not need to run these simulations.  I want you to think about the research question to discuss them at the follow-up session.  I want you to think about how you can use these methods in your research as if you are not doing so, attending these masterclasses is pointless as you will never really apply these methods.  
