# REU2018

Hello! This is the repository for Lauren K Borland and Collin Mulcahy's Ecological Data Science REU project at the University of Wisconsin at La Crosse. I (Lauren) did a lot of analysis and data exploration in R for this project, including:
* calling the adehabitat package in R to create an object of class ltraj to get movement descriptors for our telemetry data
* removing instances where we had more than one location for one fish at one time, and where time gaps were over 5 seconds
* creating transition matrices for each period of the trial and then each fish in these trials
* creating functions for first- and second-order Markov chains to simulate fish movement in study pools
* conducting Moran's I Monte Carlo analysis for spatial autocorrelation of averages of movement metrics in study tanks
* creating a metric for analyzing if there was a difference between pre-co2 and during-co2 transition matrices by
	* calculating row-wise total variance (variance of movement choice probabilities based on each starting point) for each fish
	* weighting total variance by proportion of time spent in each specific grid cell / starting point
	* summing total variances to get a single metric to analyze the signal (changes from pre-co2 to during-co2) over the noise (changes between fish in the periods)

https://github.com/lkborland/REU2018/blob/master/CarpMarkov.R is where most of my work is currently.
:bulb:
