library(sp)
library(methods)
library(ade4)
library(CircStats)
library(stats)
library(plyr)
library(dplyr)
library(adehabitatHR)
library(ggplot2)


#carp data: type II
#TRIAL 1
#make variable where duplicate observations of TagCode and DateTime are removed
t1 <- trial1[!duplicated(trial1[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial1 <- t1 %>% mutate(dt1 = as.POSIXct(DateTime))


#create object of class ltraj with info about locations (Period, Period2, Species, Trial)
carp1 <- as.ltraj(xy = t1[,c("Easting","Northing")], date = dttrial1$dt1,
                  id = t1$TagCode, typeII = TRUE,
                  infolocs = t1[,7:10])
carp1

#make data regular
#date.ref1 <- as.POSIXct(c("2014-10-15 07:00:07", "2014-10-15 07:00:12", "2014-10-15 07:00:44", "2014-10-15 07:00:17", 
#                          "2014-10-15 07:00:07", "2014-10-15 07:00:13", "2014-10-15 07:00:10", "2014-10-15 07:00:20", 
#                          "2014-10-15 07:00:12", "2014-10-15 07:00:08"))
#reg_c1 <- sett0(carp1, date.ref1, 3.5)


#TRIAL 2
t2 <- trial2[!duplicated(trial2[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial2 <- t2 %>% mutate(dt2 = as.POSIXct(DateTime))


#create object of class ltraj
carp2 <- as.ltraj(xy = t2[,c("Easting","Northing")], date = dttrial2$dt2,
                  id = t2$TagCode, typeII = TRUE,
                  infolocs = t2[,7:10])
carp2

#TRIAL 3
t3 <- trial3[!duplicated(trial3[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial3 <- t3 %>% mutate(dt3 = as.POSIXct(DateTime))


#create object of class ltraj
carp3 <- as.ltraj(xy = t3[,c("Easting","Northing")], date = dttrial3$d3,
                  id = t3$TagCode, typeII = TRUE,
                  infolocs = t3[,7:10])
carp3


#summarize carp by trial with basic trajectory calculations
#make df from carp1
carp1df <- ldply(carp1)
#make order of observations in order to join with infolocs
carp1df <- mutate(carp1df, num = 1:611514)

#add infolocs to carp1df
carp1info <- ldply(infolocs(carp1))
#make order of observations in order to join with carp1df
carp1info <- mutate(carp1info, num = 1:611514)
carp1df <- left_join(carp1df, carp1info, by = "num")
#organize by Period
carp1df$Period <- factor(carp1df$Period, levels = c("PreCO2", "DuringCO2", "PostCO2"))

#separate fish 1 from carp1df as example
carp1.1 <- carp1df[1:73039,]
carp1.1$Period <- factor(carp1.1$Period, levels = c("PreCO2", "DuringCO2", "PostCO2"))

#subset carp1.1 by Period as example
carp1.1Pre <- filter(carp1.1, Period == "PreCO2")
carp1.1During <- filter(carp1.1, Period == "DuringCO2")
carp1.1Post <- filter(carp1.1, Period == "PostCO2")


#filter out bad rows
carp1dirty <- carp1df
badrows <- (which(carp1df$dt > 5))
carp1df[carp1df$badrows+1,"rel.angle"] <- NA
carp1c <- carp1df[-badrows,]



barrier <- c(2.25, 3.1)
carp1c$zone <- ifelse(carp1c$x < barrier[1] & carp1c$y > barrier[2], "High", "Low")
table(carp1c$zone[-1],carp1c$zone[-length(carp1c$zone)]) / length(carp1c$zone)

#making boxes in pool
easting.width <- 4.9
northing.width <- 10
easting.zones <- 4
northing.zones <- 8
easting.boundaries <- sort(c(0, barrier[1], easting.width/seq(1:easting.zones)))
northing.boundaries <- sort(c(0, northing.width/seq(1:northing.zones)))

#creating zones in pool, searching for locations (x values)
lenx <- length(carp1c$x)
x.zones <- NULL

for(i in 1:lenx){
  for(j in 2:length(easting.boundaries)){
    if(carp1c$x[i] > easting.boundaries[j-1] & carp1c$x[i]<= easting.boundaries[j]) {
      x.zones[i] <- paste0("e",j-1)
    }
  }
}

#creating zones in pool, searching for locations (y values)
leny <- length(carp1c$y)
y.zones <- NULL

for(i in 1:leny){
  for(j in 2:length(northing.boundaries)){
    if(carp1c$y[i] > northing.boundaries[j-1] & carp1c$y[i]<= northing.boundaries[j]) {
      y.zones[i] <- paste0("n",j-1)
    }
  }
}

#create variable for exact grid the fish is in
carp1c$x.zones <- x.zones
carp1c$y.zones <- y.zones
carp1c <- carp1c %>% mutate(grid = paste0(x.zones, y.zones))

#create boxplot displaying movement segment distance
#for fish 1 trial 1
carp1.1mvmt <- ggplot(carp1.1, aes(x = Period, y = dt, fill = Period)) + 
  geom_boxplot() +
  scale_fill_manual("Period", values = c("sienna1", "sienna3", "sienna4")) +
  theme_bw() +
  ylab("Movement Segment Distance") +
  xlab("Period of CO2 Treatment") +
  theme(
    strip.background = element_blank()
  ) 
#for all fish in trial 1
carp1mvmt <- ggplot(carp1df, aes(x = Period, y = dt, fill = Period)) + 
  geom_boxplot() +
  scale_fill_manual("Period", values = c("sienna1", "sienna3", "sienna4")) +
  theme_bw() +
  ylab(expression("Log"[10]*" Movement Segment Distance")) +
  xlab(expression("Period of CO"[2]*" Treatment")) +
  ggtitle("Distance of Segment Movement by Period\nTrial One") +
  theme(
    strip.background = element_blank()
  ) 
#scale on log10
carp1mvmt + scale_y_log10()
#testing differences
mod.mvmt <- lm(dist ~ Period, data=carp1df)
summary(mod.mvmt)
anova(mod.mvmt)


#create boxplot displaying relative angle (absolute values) trial 1
#make new variable for absolute value of relative angle
carp1ARA <- mutate(carp1df, absrelang = abs(rel.angle))
carp1ang <- ggplot(carp1ARA, aes(x = Period, y = absrelang, fill = Period)) + 
  geom_boxplot() +
  scale_fill_manual("Period", values = c("sienna1", "sienna3", "sienna4")) +
  theme_bw() +
  ylab("Absolute Value of Relative Angle") +
  xlab(expression("Period of CO"[2]*" Treatment")) +
  ggtitle("Relative Angle Movement by Period\nTrial One") +
  theme(
    strip.background = element_blank()
  ) 
#testing differences
mod.ARA <- lm(absrelang ~ Period, data=carp1ARA)
summary(mod.ARA)


#create boxplot for acceleration, trial 1
#first make variable for acceleration, reduce to accelerations
carp1acc <- mutate(carp1df, acc = dist/((dt)^2))
carp1accplot <- ggplot(carp1acc, aes(x = Period, y = acc, fill = Period)) + 
  geom_boxplot() +
  scale_fill_manual("Period", values = c("sienna1", "sienna3", "sienna4")) +
  theme_bw() +
  ylab("Relocation Acceleration") +
  xlab(expression("Period of CO"[2]*" Treatment")) +
  ggtitle("Acceleration of a Movement by Period\nTrial One") +
  theme(
    strip.background = element_blank()
  ) 
#testing differences
mod.acc <- lm(acc ~ Period, data=carp1acc)
summary(mod.acc)



#create cumulative distance plot trial 1 fish 1
carp1cdist <- ggplot(carp1.1, aes(x = date, y = cumsum(dist), color = Period)) + 
  geom_line(aes(size = Period)) +
  scale_color_manual("Period", values = c("sienna1", "sienna3", "sienna4"))+
  scale_size_manual("Period", values = c(2.5,2.5,2.5))+
  theme_bw() +
  ylab("Cumulative Distance") +
  xlab("Time") +
  ggtitle("Cumulative Distance Traveled Over Time\nTrial One, Fish One")



#null model (multiple)
#define a function to plot randomized trajectory over study area (outdoor pool)
#carp1nmm <- NMs.CRW(N=10, nlocs=50000, nrep=1)


#correlogram: needs regular data
#corrgram1 <- acfdist.ltraj(carp1, which = c("dist", "dx", "dy"), nrep = 999, lag = 1,
#              plot = TRUE, xlab = "Lag", ylab = "autocorrelation")
#corrgram1