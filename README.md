# Evolutionary Computing

##### By Mirjam Bruinsma, Guy Frankel, Thomas van der Veen and Walter Vianen

This repository is used during the course Evolutionary Computing at the University of Amsterdam in September and October 2018.

We devoloped an evolutionary algorithm that belongs to the class of the island models. This is a model in which the population is divided over multiple interacting islands.
Each island contains their own EA, that is looking through the solution space. After a number of generations they migrate some individuals, to ensure a greater variation in each island. 
The aim of this model is to research the effect of this structural island model on diversity and quality within the subpopulations and the population as a whole.

## File structure
The main model can be found in the 'Island' directory. This folder contains the three algorithms of which the model consists and the way they interact.
The analysis we performed on the three algorithms separately, can be found in the 'Analysis' folder. The analysis on the island model can be found in the 'Island' folder.
The 'Benchmark folder contains the algorithms optimised for winning a competition, related to the assignment of this course. 

## Reproduce figures
To reproduce the figures used in the report, one must run the 'multiple.sh' bash files available in each algorithm folder. This stores the raw data in text files. Consecutively these files can be analysed in python with 'analysis.py'. Plots of the Island model can be made by first running the 'Saving_results_Schaffers' jupyter notebook and then the 'Making_plots_Schaffers' jupyter notebook. 
Be aware that the raw data takes up a minimum of 50mb per seed.
