# Examining the impacts of a carbon dioxide barrier in outdoor ponds to the behavioral responses of bighead and silver carp

## Purpose of code

This code develops new methods for analyzing bighead and silver carp movement data, which was previously published by [Cupp et al. 2016](https://doi.org/10.1139/cjfas-2015-0472).

## Code authorship

Lauren K. Borland primarily wrote the code with assistance from Collin Mulcahy. 
Richard A. Erickson edited the repository and adapted the code to create publication quality figures. 


## Overview of the code

Lauren wrote the code that analyses and data exploration in R, including:
* calling the adehabitat package in R to create an object of class ltraj to get movement descriptors for our telemetry data
* removing instances where we had more than one location for one fish at one time, and where time gaps were over 5 seconds
* creating transition matrices for each period of the trial and then each fish in these trials
* creating functions for first- and second-order Markov chains to simulate fish movement in study pools
* conducting Moran's I Monte Carlo analysis for spatial autocorrelation of averages of movement metrics in study tanks
* creating a metric for analyzing if there was a difference between pre-co2 and during-co2 transition matrices by
	* calculating row-wise total variance (variance of movement choice probabilities based on each starting point) for each fish
	* weighting total variance by proportion of time spent in each specific grid cell / starting point
	* summing total variances to get a single metric to analyze the signal (changes from pre-co2 to during-co2) over the noise (changes between fish in the periods)

## Files 

This repository contains the following files:

- https://github.com/lkborland/REU2018/blob/master/CarpMarkov.R 
- (Will add more, with description as RAE reviews them)
