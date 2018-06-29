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


#show where barrier is, make high and low zones of CO2 distinction
barrier <- c(2.25, 3.1)
carp1c$zone <- ifelse(carp1c$x < barrier[1] & carp1c$y > barrier[2], "High", "Low")
table(carp1c$zone[-1],carp1c$zone[-length(carp1c$zone)]) / length(carp1c$zone)

#making boxes in pool
easting.width <- 4.9
northing.width <- 10
easting.zones <- 4
northing.zones <- 8
easting.boundaries <- sort(c(0, (seq(1:easting.zones)*(easting.width/easting.zones))))
northing.boundaries <- sort(c(0, (seq(1:northing.zones)*(northing.width/northing.zones))))

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

#transition matrix for all periods all fish, trial 1
table(carp1c$grid[-1], carp1c$grid[-length(carp1c$grid)])
carp1matrix <- table(carp1c$grid[-1], carp1c$grid[-length(carp1c$grid)]) / length(carp1c$grid)

#general Markov chain for all periods
carp1mchain <- function (nn, transition.matrix, start=sample(1:nrow(transition.matrix), 1)) {
  output <- rep (NA, nn)
  output[1] <- start
  for (mvmt in 2:nn) 
    output[mvmt] <- sample(ncol(transition.matrix), 1, prob=transition.matrix[output[mvmt-1],])
  #print(table(output))
  #print(summary(output))
  output
}
#simulation with 1000 relocations for all periods
carp1mchain(1000, carp1matrix)

#divide observations by period
carp1cPre <- carp1c %>% filter(Period == "PreCO2")
carp1cDur <- carp1c %>% filter(Period == "DuringCO2")
carp1cPost <- carp1c %>% filter(Period == "PostCO2")

#transition matrix for PreCO2 Period
carp1mPre <- table(carp1cPre$grid[-1], carp1cPre$grid[-length(carp1cPre$grid)]) / length(carp1cPre$grid)
#general Markov chain for PreCO2 period simulation
prechain <- carp1mchain(200000, carp1mPre)

#transition matrix for DuringCO2 Period
carp1mDur <- table(carp1cDur$grid[-1], carp1cDur$grid[-length(carp1cDur$grid)]) / length(carp1cDur$grid)
#general Markov chain for DuringCO2 period simulation
durchain <- carp1mchain(200000, carp1mDur)

#transition matrix for PostCO2 Period
carp1mPost <- table(carp1cPost$grid[-1], carp1cPost$grid[-length(carp1cPost$grid)]) / length(carp1cPost$grid)
#general Markov chain for PostCO2 period simulation
postchain <- carp1mchain(200000, carp1mPost)


#total time spent in each grid cell
time_gridPre <- carp1cPre %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridDur <- carp1cDur %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridPost <- carp1cPost %>% group_by(grid) %>% mutate(sumtime = sum(dt))

#make box "map" to show study pool with grid cells
boxes <- data.frame(read.table(text = "0 0 10 0 
                               0 0 0 4.9 
                               10 0 0 4.9
                               10 4.9 -10 0 
                               10 2.25 -6.9 0
                               10 2.5 -6.9 0
                               3.1 2.25 0 .25
                               4 0 0 2.25", header = F))
gridlines <- data.frame(read.table(text = "0 1.225 10 0
                                   0 2.45 10 0
                                   0 3.68 10 0
                                   1.25 0 0 4.9
                                   2.5 0 0 4.9
                                   3.75 0 0 4.9
                                   5 0 0 4.9
                                   6.25 0 0 4.9
                                   7.5 0 0 4.9
                                   8.75 0 0 4.9", header = F))
gridlines <- setNames(gridlines, c("Northing", "Easting", "delta_lat", "delta_long"))
boxes <- setNames(boxes, c("Northing", "Easting", "delta_lat", "delta_long"))
ggplot() + 
  theme_bw() + 
  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat)) + 
  geom_segment(data=gridlines, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat), color = "darkblue",
               size = 1.5)

#create density 'heatmap' for counts in each grid for a simulation by Period
#PreCO2 period
grid.precount <- table(prechain)
raster.gridpre <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpre$w <- easting.width/easting.zones
raster.gridpre$z <-  factor(grid.precount)
raster.gridpre <- data.frame(raster.gridpre)
cc <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=32))
boxes$z = NA
ggplot(raster.gridpre, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Simualated Relocations\nPreCO2")
#DuringCO2 period
grid.durcount <- table(durchain)
raster.griddur <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.griddur$w <- easting.width/easting.zones
raster.griddur$z <-  factor(grid.durcount)
raster.griddur <- data.frame(raster.griddur)
ggplot(raster.griddur, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Simualated Relocations\nDuring CO2")
#Post CO2 Period
grid.postcount <- table(postchain)
raster.gridpost <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpost$w <- easting.width/easting.zones
raster.gridpost$z <-  factor(grid.postcount)
raster.gridpost <- data.frame(raster.gridpost)
ggplot(raster.gridpost, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Simualated Relocations\nPost CO2")

#create density 'heatmap' for total amount of time in each grid by Period (OG data)
#Pre
grid.pretime <- time_gridPre %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridpret <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpret$w <- easting.width/easting.zones
raster.gridpret$z <-  factor(grid.pretime$sumtime)
raster.gridpret <- data.frame(raster.gridpret)
ggplot(raster.gridpret, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nPre CO2")
#During
grid.durtime <- time_gridDur %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.griddurt <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.griddurt$w <- easting.width/easting.zones
raster.griddurt$z <-  factor(grid.durtime$sumtime)
raster.griddurt <- data.frame(raster.griddurt)
ggplot(raster.griddurt, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nDuring CO2")
#Post
grid.posttime <- time_gridPost %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridpostt <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpostt$w <- easting.width/easting.zones
raster.gridpostt$z <-  factor(grid.posttime$sumtime)
raster.gridpostt <- data.frame(raster.gridpostt)
ggplot(raster.gridpostt, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nPost CO2")

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
