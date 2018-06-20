library(sp)
library(methods)
library(ade4)
library(CircStats)
library(stats)
library(plyr)
library(dplyr)
library(adehabitatHR)


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
                          "2014-10-15 07:00:07", "2014-10-15 07:00:13", "2014-10-15 07:00:10", "2014-10-15 07:00:20", 
                          "2014-10-15 07:00:12", "2014-10-15 07:00:08"))
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


#summarize carp by trial ? with basic trajectory calculations
#make df from carp1
carp1df <- ldply(carp1)
#make order of observations in order to join with infolocs
carp1df <- mutate(carp1df, num = 1:611514)

#add infolocs to carp1df
carp1info <- ldply(infolocs(carp1))
#make order of observations in order to join with carp1df
carp1info <- mutate(carp1info, num = 1:611514)
carp1df <- left_join(carp1df, carp1info, by = "num")

#separate fish 1 from carp1df
carp1.1 <- carp1df[1:73039,]

#subset carp1.1 by Period
carp1.1Pre <- filter(carp1.1, Period == "PreCO2")
carp1.1During <- filter(carp1.1, Period == "DuringCO2")
carp1.1Post <- filter(carp1.1, Period == "PostCO2")


#null model (multiple)
#define a function to plot randomized trajectory over study area (outdoor pool)
#carp1nmm <- NMs.CRW(N=10, nlocs=50000, nrep=1)



#correlogram: needs regular data
#corrgram1 <- acfdist.ltraj(carp1, which = c("dist", "dx", "dy"), nrep = 999, lag = 1,
#              plot = TRUE, xlab = "Lag", ylab = "autocorrelation")
#corrgram1