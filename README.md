# Examining the impacts of a carbon dioxide barrier in outdoor ponds to the behavioral responses of bighead and silver carp

## Purpose of code

This code develops new methods for analyzing bighead and silver carp movement data, which was previously published by [Cupp et al. 2016](https://doi.org/10.1139/cjfas-2015-0472).
Specifically the code includes develops a Markov chain to describe the movement of carp and then uses the Markov chain with [Moran's I test](https://en.wikipedia.org/wiki/Moran%27s_I) to examine spatial correlation. 
The ultimate goal of the project seeks to examine if a carbon dioxide barrier impacts Asian carp behavior. 

## Code authorship

- [Lauren K. Borland](https://github.com/lkborland/) primarily wrote the code with assistance from Collin Mulcahy. 
- [Richard A. Erickson](https://github.com/rerickson-usgs) edited the repository and adapted the code to create publication quality figures. 


## Description of code files

This repository contains the following files:

- [`organizeData.R`](https://github.com/lkborland/REU2018/blob/master/organizeData.R), which merges together data files and calculates distances traveled. 
- [`visualize_and_regression.R`](https://github.com/lkborland/REU2018/blob/master/visualize_and_regression.R), which plots and runs regressions on the 75th quantiles.
- [`MarkovModel.R`](https://github.com/lkborland/REU2018/blob/master/MarkovModel.R), which builds a first order Markov Chain model, run the model, and plots the results. 
- [`markovFunctions.R`](https://github.com/lkborland/REU2018/blob/master/MoranIfunctions.R), which contains helper files for the Markov Model script. 
- [`MoranItest.R`](https://github.com/lkborland/REU2018/blob/master/MoranItest.R), uses Moran's I test to examine if a spatial correlation exists among cells.
- [`MoranIfunctions.R`](https://github.com/lkborland/REU2018/blob/master/MoranIfunctions.R), which contains helper files for the Moran's I test script. 


The files are listed in rough order of analyzes and increasing analyzes complexity. 
