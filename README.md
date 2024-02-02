# ST2 (Théorie des Jeux) - EI Algorithmique Génétique (Group 1)

![CentraleSupelec Logo](https://www.centralesupelec.fr/sites/all/themes/cs_theme/medias/common/images/intro/logo_nouveau.jpg)

## Context
All the cells that make up life on Earth contain one or more DNA molecules that carry genetic information. These molecules, which vary in length, are made up of a succession of nucleotides (or bases: A, C, G and T) which interact with numerous cellular elements, and whose position in space plays an important role in the cell's adaptation to its environment (heat, famine, stress...). While DNA sequences are now widely studied through their textual sequence (succession of A, C, G and T), it is highly instructive to study them through their three-dimensional trajectory. In 1993, biophysicists established a 3D conformation model that transforms a sequence of nucleotides (in the form of letters) into a three-dimensional trajectory. It is now possible to represent any textual DNA sequence as a 3D trajectory.

<img src="documents/RotTable.png" alt="Rotation Table" width="45%"/><img src="documents/Traj3D.png" alt="3D Trajectory Building" width="55%"/>

## Problem
The model was developed for short, naked DNA sequences, it does not take into account all the characteristics of a long chain within the cell (supercoils, nucleosomes, long-distance interactions, etc.). For example, if we observe a bacterial chromosome (a long DNA sequence making up a bacterium) or a plasmid (a small sequence present within bacteria), we will see that this chromosome or plasmid is circular, i.e. the two ends have been "glued" together. The above-mentioned model does not account for this phenomenon when representing the 3D trajectory of a bacterial chromosome or plasmid.

## Objective Function
We choose as our objective function the distance between the start and end point. The Problem is equivalent to the optimization problem of minimizing the distance as a function of the table parameters.

## Code
In order to solve this problem we applied to different algorithms: Recuit simulé and a Genetic Algorithm, its implementation can be visited in the respective python files genetic_algorithm.py, recuit_simule.py. Some minor changes were done in RotTable and Traj3D, we added new methods to facilitate the calculation of our objective function (The distance between start and end-point), and the access and exportation of the tables that we worked on.

### Initialization and Changing the parameters of the Algorithms:

The recuit simulé algorithm takes is implemented as a function <code>recuit_simule()</code>, which takes the following variables in order: dna sequence, a trajectory, initial temperature and max running time (the 2 last ones greater than 0), and a coefficient coeff (see report).  The defaults parameters can be changed on the lines 20 and 21 in __main__, in the section Parameters for Simulated Annealing.

On the other hand the genetic algorithm was implemented as a class <code> GeneticAlgorithm() </code> , it takes the following parameters in order: population size, table of rotations, mutation rate (between 0 and 1), a genetic sequence and a trajectory. The default parameters can be changed on the lines 14 and 15 in __main__, in the section Paramters for Genetic Algorithm.

### Exécution
<code>python -m dna < algorithm_choice > < filepath_to_data > </code>

Where the data can be plasmid_8k.fasta or plasmid_180k.fasta. And the algorithm_choice can be <code> 1 </code> for Recuit simulé or <code> 2 </code> for the genetic algorithm.
