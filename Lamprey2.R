# lamprey data
setwd("E:/REU 2018")

#setting data frame for csv file of lamprey data
lamprey <- read.csv("lampreyDataUse.csv")

#summarize data
summary(lamprey)

library(ggplot2)
library(Hmisc)
library(dplyr)
library(reshape)
library(lme4)
library(nlme)



#density plots by trial type
cont <- filter(lamprey, treatmentTrial == "control")
l.nic <- filter(lamprey, treatmentTrial == "nic" & treatSide == "Right control")
r.nic <- filter(lamprey, treatmentTrial == "nic" & treatSide == "Left control")
l.tfm <- filter(lamprey, treatmentTrial == "tfm" & treatSide == "Right control")
r.tfm <- filter(lamprey, treatmentTrial == "tfm" & treatSide == "Left control")
both <- filter(lamprey, treatmentTrial == "tfmnic")

dens.l.nic <- ggplot(l.nic, aes(x=x, y=y)) + geom_density2d() + 
  stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
  scale_fill_gradient("Relative\nDetection\nDensity",low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
  ggtitle("Left Nic")


dens.r.nic <- ggplot(r.nic, aes(x=x, y=y)) + geom_density2d() + 
  stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
  scale_fill_gradient("Relative\nDetection\nDensity",low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
  ggtitle("Right Nic")


dens.l.tfm <- ggplot(l.tfm, aes(x=x, y=y)) + geom_density2d() + 
  stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
  scale_fill_gradient("Relative\nDetection\nDensity",low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
  ggtitle("Left TFM")


dens.r.tfm <- ggplot(r.tfm, aes(x=x, y=y)) + geom_density2d() + 
  stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
  scale_fill_gradient("Relative\nDetection\nDensity",low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
  ggtitle("Right TFM")


dens.cont <- ggplot(cont, aes(x=x, y=y)) + geom_density2d() + 
  stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
  scale_fill_gradient("Relative\nDetection\nDensity",low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
  ggtitle("Control")

dens.both <- ggplot(both, aes(x=x, y=y)) + geom_density2d() + 
        stat_density2d(aes(fill = ..level..),size = 0.01, bins = 15, geom = 'polygon') + 
        scale_fill_gradient("Relative\nDetection\nDensity", low = 'lightblue', high = 'navyblue') + xlim(45,285) + ylim(45,236) + 
        ggtitle("TFM and Nic")

library(gridExtra)
grid.arrange(dens.cont, dens.l.nic, dens.r.nic, dens.l.tfm, 
             dens.r.tfm, dens.both, top = "Density of lamprey \ndetections by treatment", 
             bottom = "X-coordinate", left = "Y-coordinate")

#make animation with boxes

head(boxes)
boxes$treatmentTrial = NA
boxes$seconds = NA
pL <- ggplot(data = l.nic[1:500,]) + 
  geom_point(aes(x, y, color = treatmentTrial, frame = seconds), alpha = 0.7) + 
  ggtitle("Left Nic \nFirst 500 Observations") +
  geom_segment(data=boxes, aes(x=x, y=y, xend = x + delta_long, yend = y + delta_lat))
pL
gganimate(pL, title_frame = T, interval = 0.2, filename = "animateL.gif")

#make scatterplots of data showing 
library(RColorBrewer)
head(boxes)
boxes$Period = NA
boxes$Time = NA
pL2 <- ggplot(data = lamprey) + 
  geom_path(aes(x, y, color = treatmentTrial), alpha = 0.7) + 
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  facet_wrap(treatSide ~ treatmentTrial )
pL2

#make detection area plot
detect_lnicleft <- l.nic %>% filter(x <= 166 )