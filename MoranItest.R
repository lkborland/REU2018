## Load source file to load required packages and
## format data
source("MarkovModel.R")
library(ape)
library(dplyr)


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

MoranIacc <- ggplot(accOut, aes(x = Period, y = observed )) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    ylab("Observed Moran's I") +
    xlab("Trial Period")
print(MoranIacc)
ggsave(file = "MoranI_acc.pdf", MoranIacc, width = 6, height = 4)

head(accOut)
accLmer <- lmer(observed ~ Period + Species + (1| TagCodeTrial), data = accOut)
summary(accLmer)

accConfInt <- data.frame(cbind(fixef(accLmer),
                    confint(accLmer)[-c(1:2),]
                    ))

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


velLmer <- lmer(observed ~ Period + Species + (1| TagCodeTrial), data = velOut)
summary(velLmer)

velConfInt <- data.frame(cbind(fixef(velLmer),
                    confint(velLmer)[-c(1:2),]
                    ))

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
distLmer <- lmer(observed ~ Period + Species+ (1| TagCodeTrial), data = distOut)

summary(distLmer)

distConfInt <- data.frame(cbind(fixef(distLmer),
                     confint(distLmer)[-c(1:2),]
                     ))

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
relLmer <- lmer(observed ~ Period + Species + (1| TagCodeTrial), data = rel.angleOut)

summary(relLmer)

relAngleConfInt <- data.frame(
    cbind(fixef(relLmer),
          confint(relLmer)[-c(1:2),]
          ))


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
absLmer <- lmer(observed ~ Period + Species  +(1| TagCodeTrial), data = abs.angleOut)
summary(absLmer)

absAngleConfInt <- data.frame(cbind(fixef(absLmer),
                                    confint(absLmer)[-c(1:2),]
                                    ))


## Merge together and plot all results

head(accOut)
accOut$EndPoint <- "Acceleration"
head(velOut)
velOut$EndPoint <-  "Velocity"
head(distOut)
distOut$EndPoint <- "Distance"
head(rel.angleOut)
rel.angleOut$EndPoint <- "Relative angle"
head(abs.angleOut)
abs.angleOut$EndPoint <- "Absolute angle"

allEndPoints <- rbind(velOut,
                      accOut,
                      distOut,
                      rel.angleOut,
                      abs.angleOut)



MoranIall <-
    ggplot(allEndPoints, aes(x = Period, y = observed)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_point() +
    theme_bw() +
    geom_hline(aes(yintercept = expected), color = 'red', size = 2) +
    facet_grid( EndPoint ~ . ) + 
    ylab("Observed Moran's I") +
    xlab("Trial Period") +
    theme(
        strip.background = element_blank()
    ) 
print(MoranIall) 

ggsave(file = "MoranI_all.pdf", MoranIall, width = 6, height = 8)

accConfInt$EndPoint <- "Acceleration"
velConfInt$EndPoint <-  "Velocity"
distConfInt$EndPoint <- "Distance"
relAngleConfInt$EndPoint <- "Relative angle"
absAngleConfInt$EndPoint <- "Absolute angle"

allConfInt <- rbind(
    accConfInt,
    velConfInt,
    distConfInt,
    relAngleConfInt,
    absAngleConfInt)

head(allConfInt)
allConfInt$Coefficient <- gsub("(\\))\\d$", "\\1", rownames(allConfInt))
allConfInt$Coefficient <- gsub("(2)\\d$", "\\1",   allConfInt$Coefficient)
allConfInt$Coefficient <- gsub("\\(|\\)|Period", "", allConfInt$Coefficient)
allConfInt$Coefficient <- gsub("SVC\\d$", "SVC", allConfInt$Coefficient)

head(allConfInt)

colnames(allConfInt)[1:3] <- c("Estimate", "L95", "U95")

rownames(allConfInt) <- 1:nrow(allConfInt)



allConfInt$Coefficient <- factor( allConfInt$Coefficient,
                                 levels = rev(unique(allConfInt$Coefficient)))

MoranICoef <- ggplot(allConfInt, aes(x = Coefficient, y= Estimate, ymin = L95, ymax = U95)) +
    geom_point() +
    geom_linerange() +
    facet_grid(EndPoint~ . ) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "red") +
    theme_bw() +
        theme(
        strip.background = element_blank()
        ) +
    ylab("Regression estimate for Moran's I") 
print(MoranICoef)

ggsave("MoranICoef.pdf", MoranICoef, width = 4, height = 8)
                         
sink("sessionInfo.txt")
print(sessionInfo())
sink()


