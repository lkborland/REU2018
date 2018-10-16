## Load necessary packages
library(ggplot2)
library(data.table)
library(lubridate)
library(adehabitatHR)
library(foreach)
library(doParallel)

## Load data
trial1 <- fread("trial1data.csv")
trial2 <- fread("trial2data.csv")
trial3 <- fread("trial3data.csv")

trial <- rbind(trial1, trial2, trial3)

## Create uniuqe tag code for each fish

trial[ , TagCodeTrial := paste(TagCode, Trial, sep = "-")]

## Order Period (for use later in plots and analysis)
trial[ , Period := factor(Period, levels = unique(Period))]
trial[ , unique(Period)]

trial[ , Period2 := factor(Period2, levels = unique(Period2))]
trial[ , unique(Period2)]


## process fish one at a time
fishies <- trial[ , unique(TagCodeTrial), ]

## write function to use in parallel function
## Duplicate funciton is inside of here because it is quicker to
## run the function on one fish at a time
## output is movement data for each fish 
convertToLTRAJ <- function(trialIn = trial, fishiesIn = fishies, index = 3){
    oneFish <- trialIn[ TagCodeTrial == fishiesIn[index], ]
    oneFish <- oneFish[ !duplicated(DateTime, TagCodeTrial),]
    oneFish[ , DateTime := ymd_hms(DateTime)]
    fishOut <- data.table(ld(as.ltraj(xy = oneFish[ , .(Easting, Northing)],
                                      id = oneFish[ , TagCodeTrial],
                                      date = oneFish[ , DateTime],
                                      typeII = TRUE,
                                      infolocs = oneFish[ , 7:11] )))
    return(fishOut)
}

## Run function in parallel 
registerDoParallel(cores = 4)
trialLocation <- foreach(i = 1:length(fishies),
                          .combine = 'rbind') %dopar% {  
    convertToLTRAJ(trial, fishies, i)
}

## examine outputs to make sure all tags present
trialLocation[ , unique(TagCodeTrial)]


## create absolute relative angle 
trialLocation[ , abs.rel.angle := abs(rel.angle)]

## calculate acceleration
head(trialLocation)

acceleration <- trialLocation[ , .(acc = dist/(dt^2)),
                               by = .(TagCodeTrial, date, Trial, Species,
                                      Period, Period2, x, y) ]
acceleration

## file to save
fwrite(x = trial, file = "trial.csv")
fwrite(x = trialLocation, file = "trialLocation.csv")
fwrite(x = acceleration, file = "acceleration.csv")
