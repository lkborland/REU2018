## Load source file to load required packages and
## format data
source("MarkovModel.R")
library(ape)


## Analyze acceleration by cell
acceleration <- fread("acceleration.csv")
acceleration[ , xCell := findInterval(x, easting.boundaries)]
acceleration[ , yCell := findInterval(y, northing.boundaries)]

accelerationByCell <-
    acceleration[ , .(acc = mean(acc, na.rm = TRUE)), by = .(xCell, yCell, Species, Period2, TagCodeTrial)]

## Create list of all cells (not all fish visited all cells)
maxCell <- acceleration[ , .(xMax =  max(xCell), yMax = max(yCell))]
allCells <- data.table(expand.grid(xCell =  1:maxCell[ , xMax], yCell = 1:maxCell[ , yMax] ))
setkey(allCells, xCell, yCell)

## loop through all fish during all trial
source("MoranIfunctions.R")
registerDoParallel(cores = 4)

accOut <- foreach( i = 1:nrow(tagTrialList),
                  .combine = 'rbind') %dopar% {
                      accFunction(index = i, tagTrialList = tagTrialList,
                                  accelerationByCell = accelerationByCell,
                                  firstOrder = firstOrder, allCells = allCells)
                  }

MoranIacc <- ggplot(accOut, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIacc)
ggsave(file = "MoranI_acc.pdf", MoranIacc, width = 6, height = 4)

head(accOut)
summary(lmer(observed ~ Period + (1| TagCodeTrial), data = accOut))

## Analyze velocity by cell
velocity <- trialLocation[ , .(acc = dist/(dt)),
                               by = .(TagCodeTrial, date, Trial, Species,
                                      Period, Period2, x, y) ]

velocity[ , xCell := findInterval(x, easting.boundaries)]
velocity[ , yCell := findInterval(y, northing.boundaries)]

velocityByCell <-
    velocity[ , .(acc = mean(acc, na.rm = TRUE)), by = .(xCell, yCell, Species, Period2, TagCodeTrial)]

## Create list of all cells (not all fish visited all cells)
maxCell <- velocity[ , .(xMax =  max(xCell), yMax = max(yCell))]
allCells <- data.table(expand.grid(xCell =  1:maxCell[ , xMax], yCell = 1:maxCell[ , yMax] ))
setkey(allCells, xCell, yCell)

## loop through all fish during all trial
source("MoranIfunctions.R")
registerDoParallel(cores = 4)

velOut <- foreach( i = 1:nrow(tagTrialList),
                  .combine = 'rbind') %dopar% {
                      accFunction(index = i, tagTrialList = tagTrialList,
                                  accelerationByCell = velocityByCell,
                                  firstOrder = firstOrder, allCells = allCells)
                  }


MoranIvel <- ggplot(velOut, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIvel)
ggsave(file = "MoranI_vel.pdf", MoranIvel, width = 6, height = 4)


head(velOut)
summary(lmer(observed ~ Period + (1| TagCodeTrial), data = velOut))

## Examine first order observations
trialCellAve <- trialLocation[ , .(abs.angle = mean(abs.angle, na.rm = TRUE),
                                   rel.angle = mean(rel.angle, na.rm = TRUE),
                                   dist      = mean(dist, na.rm = TRUE)),
                              by = .(x.zone, y.zone, Species, Period2, TagCodeTrial)]


## Examine distance
distOut <- foreach( i = 1:nrow(tagTrialList),
                  .combine = 'rbind') %dopar% {
                      firstFunction(index = i, tagTrialList = tagTrialList,
                                    trialCellAve = trialCellAve,
                                    firstOrder = firstOrder, allCells = allCells,
                                    varUse = "dist")
                  }


MoranIvel <- ggplot(distOut, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIvel)
ggsave(file = "MoranI_dist.pdf", MoranIdist, width = 6, height = 4)


head(distOut)
summary(lmer(observed ~ Period + (1| TagCodeTrial), data = distOut))


## Examine relative angle
rel.angleOut <- foreach( i = 1:nrow(tagTrialList),
                  .combine = 'rbind') %dopar% {
                      firstFunction(index = i, tagTrialList = tagTrialList,
                                    trialCellAve = trialCellAve,
                                    firstOrder = firstOrder, allCells = allCells,
                                    varUse = "rel.angle")
                  }


MoranIvel <- ggplot(rel.angleOut, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIvel)
ggsave(file = "MoranI_rel.angle.pdf", MoranIrel.angle, width = 6, height = 4)


head(rel.angleOut)
summary(lmer(observed ~ Period + (1| TagCodeTrial), data = rel.angleOut))


## Examine abs angle
abs.angleOut <- foreach( i = 1:nrow(tagTrialList),
                  .combine = 'rbind') %dopar% {
                      firstFunction(index = i, tagTrialList = tagTrialList,
                                    trialCellAve = trialCellAve,
                                    firstOrder = firstOrder, allCells = allCells,
                                    varUse = "abs.angle")
                  }


MoranIvel <- ggplot(abs.angleOut, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIvel)
ggsave(file = "MoranI_abs.angle.pdf", MoranIabs.angle, width = 6, height = 4)


head(abs.angleOut)
summary(lmer(observed ~ Period + (1| TagCodeTrial), data = abs.angleOut))


sink("sessionInfo.txt")
print(sessionInfo())
sink()


