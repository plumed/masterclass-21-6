#  PLUMED Masterclass 21.6: Dimensionality reduction

This lesson was given as part of the PLUMED masterclass series in 2021.  It includes:

* A video that explain the theory covered and a second video which shows you how the exercises should be completed.
* A series of exercises that you should try to complete yourself.
* Some supplementary python notebooks that provide further background information on the exercise.

The flow chart shown below indicates the order in which you should consult the resources.  You can click on the nodes to access the various resources.  Follow the thick black lines for the best results.  The resources that are connected by dashed lines are supplmentary resources that you may find useful when completing the exercise.  There is a lot of material in this lesson and the flow chart shows you multiple paths that you can follow through this material.  The best path to follow will depend on your interests. 

This lesson was the sixth masterclass in the 2021 series.  You will likely be able to complete the exercise without completing all the five masterclasses before this one in the series.  However, the first masterclass contains instructions for installing PLUMED using conda that you may need to consult if you have not installed it already.  

```mermaid
flowchart LR;
  A[PLUMED intro] -.-> B[Lecture I] 
  B ==> C[Getting started]
  subgraph one [Introduction];
  C -.-> D[Solutions]
  end
  one --> E[Dimensionality reduction]
  one --> F[Path CVs]
  one --> G[Indistinguishability]
  H ==> I[Lecture II]
  subgraph two [Dimensionality reduction]
  E -.-> J[Solutions]
  end
  subgraph three [Path CVs]
  F -.-> K[Solutions]
  F -.-> M[Generating startpoints]
  M -.-> F
  end
  subgraph four [Indistinguishability]
  G -.-> L[Solutions]
  end
  two --> H[Your research] 
  three --> H
  four --> H
  click A "ref1" "This elementary lesson shows you how to install PLUMED and some basic and more advanced syntax. You will only need to complete up to exercise 1 of this earlier lesson to complete these exercises on dimensionality reduction"
  click B "video1" "A lecture that was given on April 12th 2021 as part of the plumed masterclass series that introduces you to the exercises in this lesson"
  click C "INSTRUCTIONS1.md" "Installation instructions followed by a discussion of how to run PLUMED from python direction and how to use chemiscope to visualize trajectories"
  click D "notebooks/plumed-and-python-exercises-1-and-2.ipynb" "This python notebook contains a full set of solutions for the exercises that show how you can call PLUMED from a jupyter notebook"
  click E "DIMENSIONALITY.md" "An introduction to running various different types of dimensionality reduction algorithms on trajectories"
  click F "PATHS.md" "An introduction to path collective variables and the isocommittor"
  click G "INDISTINGUISHABILITY.md" "An introduction to methods for dealing with systems of indistinguishable particles"
  click H "WRAPUP.md" "Some advice on using the methods discussed in this tutorial in your own research"
  click I "video2" "A lecture that was given on April 19th 2021 as part of the plumed masterclass series that goes through the solutions to the exercises in the lesson"
  click J "notebooks/dimensionality-reduction-exercise-3.ipynb" "This python notebook contains a full set of solutions for the exercises on dimensionality reduction techniques"
  click K "notebooks/path-cvs-exercise-4.ipynb" "This python notebook contains a full set of solutions for the exercises on path collective variables"
  click L "notebooks/indistinguishability-ex-5.ipynb" "This python notebook contains a full set of solutions for the exercises on dealing with systems of indistiguishable particles"
  click M "notebooks/create_startpoints.ipynb" "This python notebook shows how I generated a series of random startpoints for the isocommittor simulations that were used in the exercise on path CVs"
```
