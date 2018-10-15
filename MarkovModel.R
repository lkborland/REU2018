## Load necessary packages
library(data.table)
library(ggplot2)
library(lmerTest)
library(foreach)
library(doParallel)


source("./markovFunctions.R")

## Load data
trialLocation <- fread("trialLocation.csv")


## Set factor orders
trialLocation[ , Trial := factor(Trial)]
trialLocation[ , Period := factor(Period, levels = unique(Period))]
trialLocation[ , Period2 := factor(Period2, levels = unique(Period2))]

## add in barrier details and zones of CO2 
barrier <- c(2.25, 3.1)
## trialLocation[ , zone := ifelse(x < barrier[1] & y > barrier[2], "High", "Low")]


#making boxes in pool
easting.width <- 4.9
northing.width <- 10
easting.zones <- 5
northing.zones <- 10

easting.boundaries  <- seq(0, easting.width,  length.out = easting.zones)
northing.boundaries <- seq(0, northing.width, length.out = northing.zones)

## Assign fish locations to grids 
trialLocation
easting.boundaries
trialLocation[ , x.zone  := findInterval( x, easting.boundaries)]
trialLocation[ , y.zone  := findInterval( y, northing.boundaries)]
trialLocation[ , grid := paste( x.zone, y.zone, sep = "-")]


## create transition table

positions <- trialLocation[ , .(first.positions = grid,
                                second.positions = lead(grid, 1)
                                ),
                           by = .( TagCodeTrial, Period, Period2, Species, Trial)]
positions[ , firstTransition  := paste( first.positions, second.positions, sep = ":")]


positions[ !is.na(second.positions),
          NfirstOrder := .N,
          by = .(TagCodeTrial, Period, Period2, Species, Trial)]



positions

## create vectorized number of transitions 
firstOrder <- positions[ !is.na(second.positions),
                        .N,
                        by = .(TagCodeTrial, Period, Period2, Species, Trial,
                               first.positions, second.positions, firstTransition)] 

## Use above results to create transition probabilities (a vectorized Markov chain)
firstOrder[ , prob :=  N / sum(N),
           by = .(TagCodeTrial, Period, Period2, Species, Trial,
                  first.positions)]

## Sum to make sure all probs are 1, which is a Markovian property 
firstOrder[ , sum(prob),
           by = .(TagCodeTrial, Period, Period2, Species, Trial,
                  first.positions)]

## Run Markov chain 1000 times on reach fish for each period
registerDoParallel(cores = 4)

tagTrialList <- firstOrder[ , .N, by = .(TagCodeTrial, Period2, Species)]
tagTrialList[ , index := 1:nrow(tagTrialList)]


firstOrderOut <- foreach(i = 1:nrow(tagTrialList),
                         .combine = 'rbind') %dopar% {  
                             firstOrderMC(
                                 tag = tagTrialList[i, TagCodeTrial],
                                 period = tagTrialList[i, Period2],
                                 n.iter = 1000,
                                 firstOrderDT = firstOrder,
                                 species = tagTrialList[i, Species]
                             )   
                         }


tagSummary <- firstOrderOut[ , .(N = sum(N)), by = .( period, fish.cell, x, y, species)]
tagSummary[ , sum(N), by = period]
tagSummary[ , Nperiod := sum(N), by = period]
tagSummary[ , prob :=  N/Nperiod, by = period]
tagSummary[ , sum(prob), by = period]

## create barrier data.frame for plotting
barrierLocation <- data.frame(x = c(0.5, 2.6),
                              y = c(4.5, 4.5))
wallLocation <- data.frame(x = c(2.5, 2.5, 0.5),
                          y = c(9.5, 4.5, 4.5 ))

## Plot Markov output 
firstOrderMarkovOut <- ggplot( ) +
    geom_raster(data = tagSummary, aes(x = x, y = y, fill = prob)) +
    facet_grid( species ~ period) +
    xlab("x cell") +
    ylab("y cell") +
    theme_bw() +
    theme(strip.background = element_blank()) + 
    scale_fill_continuous("Probability") +
    geom_path( data = wallLocation, aes(x = x, y = y), size = 3, color = "white") +
    geom_line( data = barrierLocation, aes(x = x, y = y), size = 5, color = "orange") 
print(firstOrderMarkovOut)
ggsave(plot = firstOrderMarkovOut, file = "firstOrderMarkov.pdf", width = 9, height = 6)
ggsave(plot = firstOrderMarkovOut, file = "firstOrderMarkov.jpg", width = 9, height = 6)


## Examine mean prob as fish stays in the same cell
## Note this would be the equivalent of examining the mean of the diag, if in vector form 
diagMean <- firstOrder[ first.positions == second.positions, .(meanNoMove = mean(prob)),
                      by = .(TagCodeTrial, Period2, Species, Trial)] 

ggSameCell <- ggplot(diagMean, aes(x = Period2,  y = meanNoMove)) +
    ## geom_boxplot() +
    geom_point() + 
    theme_bw() +
    ylab("Probability of staying in same cell")  +
    xlab("Time period") +
    ylim(c(0,1)) +
    geom_violin(draw_quantiles = 0.5, fill = NA)

print(ggSameCell)
ggsave("sameCell.pdf", ggSameCell, width = 6, height = 4)
ggsave("sameCell.jpg", ggSameCell, width = 6, height = 4)

diaMeanLmer <- diagMean[ , lmer( meanNoMove ~ Period2 + Species + (1 | Trial)) ] 
summary(diaMeanLmer)


