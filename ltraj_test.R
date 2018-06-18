library(sp)
library(methods)
library(ade4)
library(adehabitatMA)
library(CircStats)
library(stats)
library(dplyr)
library(adehabitatHR)


#carp data: type II
#make variable where duplicate observations of TagCode and DateTime are removed
t1 <- trial1[!duplicated(trial1[c("TagCode","DateTime")]), ]

#change date and time to class POSIXct
dttrial1 <- as.POSIXct(t1$DateTime)


#create object of class ltraj
carp1 <- as.ltraj(xy = t1[,c("Easting","Northing")], date = dttrial1,
                  id = t1$TagCode, typeII = TRUE,
                  infolocs = t1[,7:10])
carp1
