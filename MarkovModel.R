## Load necessary packages
library(dplyr)
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

trialLocation[ , .N, by = .(TagCodeTrial, grid)][ , quantile(N, 0.01)]
trialLocation[ , .N, by = .(TagCodeTrial, grid)][ , mean(N)]
trialLocation[ , .N, by = .(TagCodeTrial, grid)][ , median(N)]
trialLocation[ , .N, by = .(TagCodeTrial, grid)][ , min(N)]

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
tagSummary[ , sum(N), by = .(period, species)]
tagSummary[ , Nperiod := sum(N), by = .(period, species)]
tagSummary[ , prob :=  N/Nperiod, by = .(period, species)]
tagSummary[ , sum(prob), by = .(period, species)]

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


trialLocation[ , PeriodSppTotal := .N, by = .(Period2, Species)]
empericalFish <-    trialLocation[ , .N,
                                  by = .(Period2, grid, x.zone, y.zone, Species)]
empericalFish[ , Nperiod := sum(N), by = .(Period2, Species)]
empericalFish[ , prob := N/sum(N), by = .(Period2, Species)]

setnames(empericalFish, c( "period", "fish.cell", "x", "y", "species", "N", "Nperiod", "prob"))

tagSummary[ , Type := "Simulated"]
empericalFish[ , Type := "Observed"]

tagSummary[ , sum(prob), by = .(period, species)]
empericalFish[ , sum(prob), by = .(period, species)]

plotAll <- rbind(tagSummary,
                 empericalFish)

plotAll[ , species_plot := factor(species,
                                  levels = c("BHC", "SVC"),
                                  labels = c("Bighead carp",
                                             "Silver carp")
                                  )]
                                  
firstOrderMarkovWithData <- ggplot( ) +
    geom_raster(data = plotAll, aes(x = x, y = y, fill = prob)) +
    facet_grid( species_plot + Type ~ period) +
    xlab("x cell") +
    ylab("y cell") +
    theme_bw() +
    theme(strip.background = element_blank()) + 
    scale_fill_continuous("Occurrence\nProbability") +
    geom_path( data = wallLocation, aes(x = x, y = y), size = 3, color = "white") +
    geom_line( data = barrierLocation, aes(x = x, y = y), size = 5, color = "orange") 

print(firstOrderMarkovWithData)

ggsave(plot = firstOrderMarkovWithData, file = "firstOrderMarkovWithData.pdf", width = 9, height = 6)
ggsave(plot = firstOrderMarkovWithData, file = "firstOrderMarkovWithData.jpg", width = 9, height = 6)


## Examine mean prob as fish stays in the same cell
## Note this would be the equivalent of examining the mean of the diag, if in vector form 
firstOrder[ , Ntct := sum(N), by = .(TagCodeTrial, Period2, Species, Trial)]
firstOrder[ , Ptct := N/Ntct]
firstOrder[ , sum(Ptct), by = .(TagCodeTrial, Period2, Species, Trial)]
diagMean <- firstOrder[ first.positions == second.positions,
                       .(meanWeigthed = mean(prob * Ptct),
                         meanNoMove = mean(prob)),
                      by = .(TagCodeTrial, Period2, Species, Trial)] 

diagMean[ , mean(meanNoMove), by = Period2]

broom::tidy(lm(meanNoMove ~  Period2 - 1, data = diagMean), conf.int = TRUE)

diagMean[ , species_plot := factor(Species,
                                   levels = c("BHC", "SVC"),
                                   labels = c("Bighead carp",
                                              "Silver carp")
                                   )]


ggSameCell <-
    ggplot(diagMean, aes(x = Period2,  y = meanWeigthed)) +
    geom_point() + 
    theme_bw() +
    theme(strip.background = element_blank()) + 
    facet_grid( . ~ species_plot) +
    ylab("Probability of staying in same cell")  +
    xlab("Time period") +
    ## ylim(c(0,1)) +
    geom_violin(draw_quantiles = 0.5, fill = NA)

print(ggSameCell)
ggsave("sameCell.pdf", ggSameCell, width = 10, height = 4)
ggsave("sameCell.jpg", ggSameCell, width = 10, height = 4)

diaMeanLmer <- diagMean[ , lmer( meanWeigthed ~ Period2 + Species + (1 | Trial)) ] 
summary(diaMeanLmer)

diaMeanLmerCI <- data.frame(cbind(fixef(diaMeanLmer), confint(diaMeanLmer)[ -c(1:2),]))
diaMeanLmerCI$Parameter <- gsub("Period2|\\(|\\)", "", rownames(diaMeanLmerCI))
colnames(diaMeanLmerCI)[1:3] <- c("Coefficient", "L95", "U95")

diaMeanLmerCI$Parameter <- factor(diaMeanLmerCI$Parameter,
                                  levels =
                                      rev(c("Intercept",
                                            "SpeciesSVC",
                                            "IncreasingCO2",
                                            "DuringCO2",
                                            "DecreasingCO2",
                                            "PostCO2")),
                                  labels =
                                      rev(c("Intercept",
                                            "Silver carp",
                                            "Increasing CO2",
                                            "During CO2",
                                            "Decreasing CO2",
                                            "Post CO2")),
                                  )



ggDiaMean <- ggplot(diaMeanLmerCI, aes(x = Parameter, y = Coefficient, ymin = L95,
                          ymax = U95)) +
    geom_point() +
    geom_linerange() +
    theme_bw() +
    coord_flip()+
    geom_hline(yintercept = 0, color = "red")
print(ggDiaMean)


ggsave("ggDiaMean.pdf", ggDiaMean, width = 4, height = 4)
ggsave("ggDiaMean.jpg", ggDiaMean, width = 4, height = 4)
