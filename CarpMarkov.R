#Markov chains for asian carp data
library(plyr)
library(dplyr)
library(ggplot2)
library(hsmm)
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
#assign fish # to rows before filtering
carp1df$fish <- NA
carp1df[1:73039,"fish"] <- 1
carp1df[73040:122976,"fish"] <- 2
carp1df[122977:182586, "fish"] <- 3
carp1df[182587:248716, "fish"] <- 4
carp1df[248717:307388, "fish"] <- 5
carp1df[307389:370901, "fish"] <- 6
carp1df[370902:444279, "fish"] <- 7
carp1df[444280:502373, "fish"] <- 8
carp1df[502374:555111, "fish"] <- 9
carp1df[555112:611514, "fish"] <- 10
carp1c <- carp1df[-badrows,]


#show where barrier is, make high and low zones of CO2 distinction
barrier <- c(2.25, 3.1)
carp1c$zone <- ifelse(carp1c$x < barrier[1] & carp1c$y > barrier[2], "High", "Low")
table(carp1c$zone[-1],carp1c$zone[-length(carp1c$zone)]) / length(carp1c$zone)

#making boxes in pool
easting.width <- 4.9
northing.width <- 10
easting.zones <- 5
northing.zones <- 10
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
carp1c$x.zones <- NULL
carp1c$y.zones <- NULL

carp1c$x.zones <- x.zones
carp1c$y.zones <- y.zones
carp1c <- carp1c %>% mutate(grid = paste0(x.zones, y.zones))
carp1.1 <- NA
carp1.1 <- carp1c[carp1c$fish==1,]

#begin transition matrix work
first.position <- carp1c$grid
second.position <- lead(carp1c$grid, 1)
positions <- data.frame(first.position, second.position)

# makes a 1-st order transition matrix trial 1

carp1matrix <- table(positions$first.position, positions$second.position)/length(positions$first.position) #makes a table

#general first-order Markov chain for all periods
carp1mchain <- function(nn, transition.matrix, start=sample(1:nrow(transition.matrix), 1)) {
  output <- rep(NA, nn)
  output[1] <- start
  for (mvmt in 2:nn) 
    output[mvmt] <- sample(1:ncol(transition.matrix), 1, prob=transition.matrix[output[mvmt-1],])
  return(output)
}

# start of second order markov 
third.position <- lead(first.position,2)
first.second <- paste(first.position, second.position)
positions$first.second <- first.second
positions$third.position <- third.position

# makes a second order transition matrix for all periods

carp1matrix2 <- table(positions$first.second, positions$third.position)/length(positions$first.second) #makes a table


#second order Markov chain for all periods
##------------------------------------
carp1mchain2 <- function(nn, transition.matrix, transition.matrix2, start = sample(1:nrow(transition.matrix), 1)) {
  temp1 <- rownames(transition.matrix)
  temp2 <- colnames(transition.matrix2)
  output <- rep(NA, nn)
  output[1] <- temp1[start]
  secondmove <- sample(ncol(transition.matrix), 1, prob=transition.matrix[output[1],])
  output[2] <- temp1[secondmove]
  bothmove <- rep(NA, nn)
  bothmove[1] <- paste(output[1], output[2])
  for (mvmt2 in 3:nn){
    look <- bothmove[mvmt2-2]
    nextmove <- mvmt2-1
    output[mvmt2] <- temp2[sample(ncol(transition.matrix2), 1, prob=transition.matrix2[look,])]
    bothmove[nextmove] <- paste(output[nextmove], output[mvmt2])
  }
  return(output)
}


convert2grid <- function(chain, easting.zones, northing.zones){
  #converts string output from the markov chain to usable grid cells
  spl <- strsplit(chain, split="e|n")
  len <- length(spl)
  grid.out <- data.frame(x = rep(1:easting.zones, each=northing.zones), y = rep(1:northing.zones, easting.zones))
  grid.out$numb <- c(1:(easting.zones*northing.zones))
  x1 <- rep(NA, len)
  y1 <- rep(NA, len)
  for(num in 1:len){
    for(e in 1:easting.zones){
      if(spl[[num]][2]==e){
        x1[num] <- e
      }
    }
    for(n in 1:northing.zones){
      if(spl[[num]][3]==n){
        y1[num] <- n
      }
    }
  }
  x1 <- x1[!is.na(x1)]
  y1 <- y1[!is.na(y1)]
  gridcell <- data.frame(xzone = x1, yzone = y1)
  grid <- rep(NA, len)
  for(assign in 1:len){
    tempx <- gridcell[assign,"xzone"]
    tempy <- gridcell[assign, "yzone"]
    grid[assign] <- grid.out %>% filter(x == tempx) %>% filter(y == tempy) %>% select(numb)
  }
  grid <- unlist(grid)
  return(grid)
}

positioning2 <- function(df){
  first.position <- df$grid
  second.position <- lead(df$grid, 1)
  positions <- data.frame(first.position, second.position)
  third.position <- lead(first.position,2)
  first.second <- paste(first.position, second.position)
  positions$first.second <- first.second
  positions$third.position <- third.position
  return(positions)
}

trans.matrix2 <- function(positions){
  matrix2 <- table(positions$first.second, positions$third.position)/length(positions$first.second) #makes a table
  return(matrix2)
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


#divide observations by period2 (pre, inc, dur, dec, post)
carp1cPre2 <- carp1c %>% filter(Period2 == "PreCO2")
carp1cInc2 <- carp1c %>% filter(Period2 == "IncreasingCO2")
carp1cDur2 <- carp1c %>% filter(Period2 == "DuringCO2")
carp1cDec2 <- carp1c %>% filter(Period2 == "DecreasingCO2")
carp1cPost2 <- carp1c %>% filter(Period2 == "PostCO2")

#transition matrix for PreCO2 Period
carp1mPre2 <- table(carp1cPre2$grid[-length(carp1cPre2$grid)], carp1cPre2$grid[-1]) / length(carp1cPre2$grid)
#transition matrix for second order
carp1m2Pre2 <- trans.matrix2(positioning2(carp1cPre2))
#general Markov chain for PreCO2 period simulation
pre2chain <- carp1mchain(200000, carp1mPre2)
#second order Markov chain
pre2chain2 <- carp1mchain2(200000, carp1mPre2, carp1m2Pre2)
pre2chain2 <- convert2grid(pre2chain2, easting.zones, northing.zones)

#transition matrix for incCO2 Period
carp1mInc2 <- table(carp1cInc2$grid[-length(carp1cInc2$grid)], carp1cInc2$grid[-1]) / length(carp1cInc2$grid)
#transition matrix second order
carp1m2Inc2 <- trans.matrix2(positioning2(carp1cInc2))
#general Markov chain for PreCO2 period simulation
inc2chain <- carp1mchain(200000, carp1mInc2)
#second order Markov chain
inc2chain2 <- carp1mchain2(200000, carp1mInc2, carp1m2Inc2)
inc2chain2 <- convert2grid(inc2chain2, easting.zones, northing.zones)

#transition matrix for DuringCO2 Period
carp1mDur2 <- table(carp1cDur2$grid[-length(carp1cDur2$grid)], carp1cDur2$grid[-1]) / length(carp1cDur2$grid)
#transition matrix second order
carp1m2Dur2 <- trans.matrix2(positioning2(carp1cDur2))
#general Markov chain for DuringCO2 period simulation
dur2chain <- carp1mchain(200000, carp1mDur2)
#second order Markov chain
dur2chain2 <- carp1mchain2(200000, carp1mDur2, carp1m2Dur2)
dur2chain2 <- convert2grid(dur2chain2, easting.zones, northing.zones)

#transition matrix for PreCO2 Period
carp1mDec2 <- table(carp1cDec2$grid[-length(carp1cDec2$grid)], carp1cDec2$grid[-1]) / length(carp1cDec2$grid)
#transition matrix second order
carp1m2Dec2 <- trans.matrix2(positioning2(carp1cDec2))
#general Markov chain for PreCO2 period simulation
dec2chain <- carp1mchain(200000, carp1mDec2)
#second order Markov chain
dec2chain2 <- carp1mchain2(200000, carp1mDec2, carp1m2Dec2)
dec2chain2 <- convert2grid(dec2chain2, easting.zones, northing.zones)

#transition matrix for PostCO2 Period
carp1mPost2 <- table(carp1cPost2$grid[-length(carp1cPost2$grid)], carp1cPost2$grid[-1]) / length(carp1cPost2$grid)
#transition matrix second order
carp1m2Post2 <- trans.matrix2(positioning2(carp1cPost2))
#general Markov chain for PostCO2 period simulation
post2chain <- carp1mchain(200000, carp1mPost2)
#second order Markov chain
post2chain2 <- carp1mchain2(200000, carp1mPost2, carp1m2Post2)
post2chain2 <- convert2grid(post2chain2, easting.zones, northing.zones)

#total time spent in each grid cell by period
time_gridPre <- carp1cPre %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridDur <- carp1cDur %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridPost <- carp1cPost %>% group_by(grid) %>% mutate(sumtime = sum(dt))

#total time spent in each grid cell by Period2
time_gridPre2 <- carp1cPre2 %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridInc2 <- carp1cInc2 %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridDur2 <- carp1cDur2 %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridDec2 <- carp1cDec2 %>% group_by(grid) %>% mutate(sumtime = sum(dt))
time_gridPost2 <- carp1cPost2 %>% group_by(grid) %>% mutate(sumtime = sum(dt))

#graphics---------------------------------------------------

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
  geom_segment(data=gridlines, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat), color = "#7BCCC4",
               size = 1.5)

#create density 'heatmap' for counts in each grid for a simulation by Period
rasterdens <- function(chain, easting.boundaries, northing.boundaries){
  #creates a df that has boundaries of the grid cells and includes a factor by count density for each grid cell
  grid.count <- table(chain)
  raster.grid <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
  raster.grid$w <- easting.width/easting.zones
  raster.grid$z <- NA
  raster.grid$z <-  factor(grid.count)
  raster.grid <- data.frame(raster.grid)
  return(raster.grid)
}

rasterplot <- function(raster.grid, cc, per = "PreCO2"){
  #creates a ggplot with raster overlay based on simulated count densities
  ggplot(raster.grid, aes(x=x, y=y, fill = z)) + 
    geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
    theme_bw() + 
    xlab("Easting") + ylab("Northing") + 
    ggtitle(sprintf("Density of Simualated Relocations\n%s", per))
}

#PreCO2 period
raster.gridpre <- rasterdens(prechain, easting.boundaries, northing.boundaries)
cc <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=easting.zones*northing.zones))
boxes$z <- NA
rasterplot(raster.gridpre, cc, "PreCO2")

#DuringCO2 period
raster.griddur <- rasterdens(durchain, easting.boundaries, northing.boundaries)
rasterplot(raster.griddur, cc, "During CO2")

#Post CO2 Period
raster.gridpost <- rasterdens(postchain, easting.boundaries, northing.boundaries)
rasterplot(raster.gridpost, cc, "Post CO2")

#create density 'heatmap' for counts in each grid for a simulation by Period2
#Pre CO2 Period2
raster.gridpre2 <- rasterdens(pre2chain, easting.boundaries, northing.boundaries)
rasterplot(raster.gridpre2, cc, "PreCO2")
#Pre CO2, second order
raster.grid2pre2 <- rasterdens(pre2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2pre2, cc, "PreCO2")
#Inc CO2 Period2
raster.gridinc2 <- rasterdens(inc2chain, easting.boundaries, northing.boundaries)
rasterplot(raster.gridinc2, cc, "Increasing CO2")
#Dur CO2 Period2
raster.griddur2 <- rasterdens(dur2chain, easting.boundaries, northing.boundaries)
rasterplot(raster.griddur2, cc, "During CO2")
#Dec CO2 Period2
raster.griddec2 <- rasterdens(dec2chain, easting.boundaries, northing.boundaries)
rasterplot(raster.griddec2, cc, "Decreasing CO2")
#Post CO2 Period2
raster.gridpost2 <- rasterdens(post2chain, easting.boundaries, northing.boundaries)
rasterplot(raster.gridpost2, cc, "Post CO2")

#create density 'heatmap' for total amount of time in each grid by Period (OG data)
#Pre
grid.pretime <- time_gridPre %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridpret <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpret$w <- easting.width/easting.zones
raster.gridpret$z <-  factor(grid.pretime$sumtime)
raster.gridpret <- data.frame(raster.gridpret)
ggplot(raster.gridpret, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
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
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
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
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nPost CO2")

#create density 'heatmap' for total amount of time in each grid by Period2 (OG data)
#Pre2
grid.pre2time <- time_gridPre2 %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridpre2t <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpre2t$w <- easting.width/easting.zones
raster.gridpre2t$z <-  factor(grid.pre2time$sumtime)
raster.gridpre2t <- data.frame(raster.gridpre2t)
ggplot(raster.gridpre2t, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nPre CO2")
#Inc2
grid.inc2time <- time_gridInc2 %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridinc2t <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridinc2t$w <- easting.width/easting.zones
raster.gridinc2t$z <-  factor(grid.inc2time$sumtime)
raster.gridinc2t <- data.frame(raster.gridinc2t)
ggplot(raster.gridinc2t, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nIncreasing CO2")
#Dur2
grid.dur2time <- time_gridDur2 %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.griddur2t <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.griddur2t$w <- easting.width/easting.zones
raster.griddur2t$z <-  factor(grid.dur2time$sumtime)
raster.griddur2t <- data.frame(raster.griddur2t)
ggplot(raster.griddur2t, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nDuring CO2")
#Dec2
grid.dec2time <- time_gridDec2 %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.griddec2t <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.griddec2t$w <- easting.width/easting.zones
raster.griddec2t$z <-  factor(grid.dec2time$sumtime)
raster.griddec2t <- data.frame(raster.griddec2t)
ggplot(raster.griddec2t, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nDecreasing CO2")
#Post2
grid.post2time <- time_gridPost2 %>% select("sumtime", "grid") %>% arrange(grid) %>% distinct()
raster.gridpost2t <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
raster.gridpost2t$w <- easting.width/easting.zones
raster.gridpost2t$z <-  factor(grid.post2time$sumtime)
raster.gridpost2t <- data.frame(raster.gridpost2t)
ggplot(raster.gridpost2t, aes(x=x, y=y, fill = z)) + 
  geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc) + 
  theme_bw() + 
  xlab("Easting") + ylab("Northing") + 
  ggtitle("Density of Time Spent\nPost CO2")
