library(sp)
library(methods)
library(ade4)
library(adehabitatMA)
library(CircStats)
library(stats)
library(dplyr)
library(adehabitatHR)


#carp data: type II
#TRIAL 1
#make variable where duplicate observations of TagCode and DateTime are removed
t1 <- trial1[!duplicated(trial1[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial1 <- as.POSIXct(t1$DateTime)


#create object of class ltraj with info about locations (Period, Period2, Species, Trial)
carp1 <- as.ltraj(xy = t1[,c("Easting","Northing")], date = dttrial1,
                  id = t1$TagCode, typeII = TRUE,
                  infolocs = t1[,7:10])
carp1

#TRIAL 2
t2 <- trial2[!duplicated(trial2[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial2 <- as.POSIXct(t2$DateTime)


#create object of class ltraj
carp2 <- as.ltraj(xy = t2[,c("Easting","Northing")], date = dttrial2,
                  id = t2$TagCode, typeII = TRUE,
                  infolocs = t2[,7:10])
carp2

#TRIAL 3
t3 <- trial3[!duplicated(trial3[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial3 <- as.POSIXct(t3$DateTime)


#create object of class ltraj
carp3 <- as.ltraj(xy = t3[,c("Easting","Northing")], date = dttrial3,
                  id = t3$TagCode, typeII = TRUE,
                  infolocs = t3[,7:10])
carp3

