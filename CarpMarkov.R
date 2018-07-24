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

convert2coord <- function(chain){
  spl <- strsplit(chain, split="e|n")
  len <- length(spl)
  output <- data.frame(x = rep(NA, len), y = rep(NA, len))
  for(i in 1:len){
    output[i, "x"] <- spl[[i]][2]
    output[i, "y"] <- spl[[i]][3]
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
carp1mPre <- table(carp1cPre$grid[-length(carp1cPre$grid)], carp1cPre$grid[-1]) / length(carp1cPre$grid)
#general Markov chain for PreCO2 period simulation
prechain <- carp1mchain(200000, carp1mPre)

#transition matrix for DuringCO2 Period
carp1mDur <- table(carp1cDur$grid[-length(carp1cDur$grid)], carp1cDur$grid[-1]) / length(carp1cDur$grid)
#general Markov chain for DuringCO2 period simulation
durchain <- carp1mchain(200000, carp1mDur)

#transition matrix for PostCO2 Period
carp1mPost <- table(carp1cPost$grid[-length(carp1cPost$grid)], carp1cPost$grid[-1]) / length(carp1cPost$grid)
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
    ggtitle(sprintf("Density of Simualated Relocations\n%s", per)) +
    theme(plot.title = element_text(size = 30)) +
    theme(axis.title = element_text(size=25))
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
#Pre CO2, second order
raster.grid2pre2 <- rasterdens(pre2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2pre2, cc, "PreCO2")
#Inc CO2 Period2 second order
raster.grid2inc2 <- rasterdens(inc2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2inc2, cc, "Increasing CO2")
#Dur CO2 Period2
raster.grid2dur2 <- rasterdens(dur2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2dur2, cc, "During CO2")
#Dec CO2 Period2
raster.grid2dec2 <- rasterdens(dec2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2dec2, cc, "Decreasing CO2")
#Post CO2 Period2
raster.grid2post2 <- rasterdens(post2chain2, easting.boundaries, northing.boundaries)
rasterplot(raster.grid2post2, cc, "Post CO2")


#create density 'heatmap' for total amount of time in each grid by Period (OG data)
rastertime <- function(time_grid, easting.boundaries, northing.boundaries){
  grid.time <- distinct(time_grid, grid)
  grid.time <- mixedsort(grid.time$grid)
  time_grid <- time_grid %>% distinct(sumtime, grid)
  time_grid$grid <- factor(time_grid$grid, levels = grid.time)
  time_grid <- arrange(time_grid, grid)
  raster.gridt <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
  raster.gridt$w <- easting.width/easting.zones
  raster.gridt$z <-  factor(time_grid$sumtime)
  raster.gridt <- data.frame(raster.gridt)
  return(raster.gridt)
}
#function for making ggplot for density of time spent
library(gtools)
rastertimeplot <- function(raster.gridt, cc, per = "Pre CO2"){
  ggplot(raster.gridt, aes(x=x, y=y, fill = z)) + 
    geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
    theme_bw() + 
    xlab("Easting") + ylab("Northing") + 
    ggtitle(sprintf("Density of Time Spent\n%s", per)) +
    theme(plot.title = element_text(size = 30)) +
    theme(axis.title = element_text(size=25))
}

#Pre
raster.gridpret <- rastertime(time_gridPre, easting.boundaries, northing.boundaries)
rastertimeplot(raster.gridpret, cc, "Pre CO2")

#During
raster.griddurt <- rastertime(time_gridDur, easting.boundaries, northing.boundaries)
rastertimeplot(raster.griddurt, cc, "During CO2")

#Post
raster.gridpostt <- rastertime(time_gridPost, easting.boundaries, northing.boundaries)
rastertimeplot(raster.gridpostt, cc, "Post CO2")


#create density 'heatmap' for total amount of time in each grid by Period2 (OG data)
#Pre2
raster.gridpre2t <- rastertime(time_gridPre2, easting.boundaries, northing.boundaries)
rastertimeplot(raster.gridpre2t, cc, "Pre CO2")

#Inc2
raster.gridinc2t <- rastertime(time_gridInc2, easting.boundaries, northing.boundaries)
rastertimeplot(raster.gridinc2t, cc, "Increasing CO2")

#Dur2
raster.griddur2t <- rastertime(time_gridDur2, easting.boundaries, northing.boundaries)
rastertimeplot(raster.griddur2t, cc, "During CO2")

#Dec2
raster.griddec2t <- rastertime(time_gridDec2, easting.boundaries, northing.boundaries)
rastertimeplot(raster.griddec2t, cc, "Decreasing CO2")

#Post2
raster.gridpost2t <- rastertime(time_gridPost2, easting.boundaries, northing.boundaries)
rastertimeplot(raster.gridpost2t, cc, "Post CO2")


#BY FISH------------------------------------------------------------------------------------------
fishmatrix2 <- function(df, numfish){
  #creates a list of  transition matrices (all 5 periods) for a specified fish number
  df <- df %>% filter(fish == numfish)
  #divide observations by period2 (pre, inc, dur, dec, post)
  pre <- df %>% filter(Period2 == "PreCO2")
  inc <- df %>% filter(Period2 == "IncreasingCO2")
  dur <- df %>% filter(Period2 == "DuringCO2")
  dec <- df %>% filter(Period2 == "DecreasingCO2")
  post <- df %>% filter(Period2 == "PostCO2")
  
  m1 <- trans.matrix2(positioning2(pre))
  m2 <- trans.matrix2(positioning2(inc))
  m3 <- trans.matrix2(positioning2(dur))
  m4 <- trans.matrix2(positioning2(dec))
  m5 <- trans.matrix2(positioning2(post))
  list.m <- list("Pre"=m1, "Inc"=m2, "Dur"=m3, "Dec"=m4, "Post"=m5)
  return(list.m)
}

#transition matrix for each fish trial 1
fishmatrix <- function(df, numfish){
  #creates a list of transition matrices (all 5 periods) for a specified fish number
  df <- df %>% filter(fish == numfish)
  #divide observations by period2
  pre <- df %>% filter(Period2 == "PreCO2")
  inc <- df %>% filter(Period2 == "IncreasingCO2")
  dur <- df %>% filter(Period2 == "DuringCO2")
  dec <- df %>% filter(Period2 == "DecreasingCO2")
  post <- df %>% filter(Period2 == "PostCO2")
  
  m1 <- table(pre$grid[-length(pre$grid)], pre$grid[-1]) / length(pre$grid)
  m2 <- table(inc$grid[-length(inc$grid)], inc$grid[-1]) / length(inc$grid)
  m3 <- table(dur$grid[-length(dur$grid)], dur$grid[-1]) / length(dur$grid)
  m4 <- table(dec$grid[-length(dec$grid)], dec$grid[-1]) / length(dec$grid)
  m5 <- table(post$grid[-length(post$grid)], post$grid[-1]) / length(post$grid)
  list.m <- list("Pre"=m1, "Inc"=m2, "Dur"=m3, "Dec"=m4, "Post"=m5)
  return(list.m)
}

#second order transition matrices for each fish in trial 1
fish1.01 <- fishmatrix2(carp1c, 1)
fish1.02 <- fishmatrix2(carp1c, 2)
fish1.03 <- fishmatrix2(carp1c, 3)
fish1.04 <- fishmatrix2(carp1c, 4)
fish1.05 <- fishmatrix2(carp1c, 5)
fish1.06 <- fishmatrix2(carp1c, 6)
fish1.07 <- fishmatrix2(carp1c, 7)
fish1.08 <- fishmatrix2(carp1c, 8)
fish1.09 <- fishmatrix2(carp1c, 9)
fish1.10 <- fishmatrix2(carp1c, 10)
fish.01 <- fishmatrix(carp1c, 1)
fish.02 <- fishmatrix(carp1c, 2)
fish.03 <- fishmatrix(carp1c, 3)
fish.04 <- fishmatrix(carp1c, 4)
fish.05 <- fishmatrix(carp1c, 5)
fish.06 <- fishmatrix(carp1c, 6)
fish.07 <- fishmatrix(carp1c, 7)
fish.08 <- fishmatrix(carp1c, 8)
fish.09 <- fishmatrix(carp1c, 9)
fish.10 <- fishmatrix(carp1c, 10)

#second order Markov chain
chain.f <- function(matrix1, matrix2){
  #get markov chain results for each period
  prechain <- carp1mchain2(50000, matrix1[[1]], matrix2[[1]])
  incchain <- carp1mchain2(50000, matrix1[[2]], matrix2[[2]])
  durchain <- carp1mchain2(50000, matrix1[[3]], matrix2[[3]])
  decchain <- carp1mchain2(50000, matrix1[[4]], matrix2[[4]])
  postchain <- carp1mchain2(50000, matrix1[[5]], matrix2[[5]])
  list.m <- list("Pre"=prechain, "Inc"=incchain, "Dur"=durchain, "Dec"=decchain, "Post"=postchain)
  return(list.m)
}

chain1.01 <- chain.f(fish.01, fish1.01)
chain1.02 <- chain.f(fish.02, fish1.02)
chain1.03 <- chain.f(fish.03, fish1.03)
chain1.04 <- chain.f(fish.04, fish1.04)
chain1.05 <- chain.f(fish.05, fish1.05)
chain1.06 <- chain.f(fish.06, fish1.06)
chain1.07 <- chain.f(fish.07, fish1.07)
chain1.08 <- chain.f(fish.08, fish1.08)
chain1.09 <- chain.f(fish.09, fish1.09)
chain1.10 <- chain.f(fish.10, fish1.10)


convert.f <- function(chain, e.z = easting.zones, n.z = northing.zones){
  #convert markov chain results into grid cells
  prechain <- convert2grid(chain[[1]], e.z, n.z)
  incchain <- convert2grid(chain[[2]], e.z, n.z)
  durchain <- convert2grid(chain[[3]], e.z, n.z)
  decchain <- convert2grid(chain[[4]], e.z, n.z)
  postchain <- convert2grid(chain[[5]], e.z, n.z)
  list.m <- list("Pre"=prechain, "Inc"=incchain, "Dur"=durchain, "Dec"=decchain, "Post"=postchain)
  return(list.m)
}

chain1.01 <- convert.f(chain1.01)
chain1.04 <- convert.f(chain1.04)
chain1.06 <- convert.f(chain1.06)
chain1.07 <- convert.f(chain1.07)
chain1.08 <- convert.f(chain1.08)
chain1.09 <- convert.f(chain1.09)
chain1.10 <- convert.f(chain1.10)

r1 <- rasterdens.f(chain1.01, easting.boundaries, northing.boundaries, period = 3)
r4 <- rasterdens.f(chain1.04, easting.boundaries, northing.boundaries, period = 3)
r6 <- rasterdens.f(chain1.06, easting.boundaries, northing.boundaries, period = 3)
r7 <- rasterdens.f(chain1.07, easting.boundaries, northing.boundaries, period = 1)
r8 <- rasterdens.f(chain1.08, easting.boundaries, northing.boundaries, period = 1)
r9 <- rasterdens.f(chain1.09, easting.boundaries, northing.boundaries, period = 1)
r10 <- rasterdens.f(chain1.10, easting.boundaries, northing.boundaries, period = 1)


#create density 'heatmap' for counts in each grid for a simulation by specified period for one fish
rasterdens.f <- function(chain, easting.boundaries, northing.boundaries, period = 1){
  #creates a df that has boundaries of the grid cells and includes a factor by count density for each grid cell
  grid.count <- table(chain[[period]])
  raster.grid <- data.frame(x = rep(easting.boundaries[-1], each=northing.zones), y = rep(northing.boundaries[-1], easting.zones))
  raster.grid$w <- easting.width/easting.zones
  raster.grid$z <- NA
  raster.grid$z <-  factor(grid.count)
  raster.grid <- data.frame(raster.grid)
  return(raster.grid)
}

rasterplot.f <- function(raster.grid, cc, per = "PreCO2"){
  #creates a ggplot with raster overlay based on simulated count densities
  ggplot(raster.grid, aes(x=x, y=y, fill = z)) + 
    geom_raster(hjust=0, vjust=0) + scale_fill_manual(values=cc, guide=FALSE) + 
    theme_bw() + 
    xlab("Easting") + ylab("Northing") + 
    ggtitle(sprintf("Density of Simualated Relocations\n%s", per)) +
    theme(plot.title = element_text(size = 30)) +
    theme(axis.title = element_text(size=25))
}


#TESTING TRANSITION MATRICES WITH TOTAL VARIANCE----------------------------------------------
tV <- function(m1, m2) {
  #creating function for getting total variation row-wise between two matrices m1 and m2
  row.out <- 0
  col.out <- 0
  for(i in 1:nrow(m1)){
    for(j in 1:ncol(m1)){
      col.out[j] <- abs(m1[i, j] - m2[i, j])
    }
    row.out[i] <- max(col.out)
  }
  return(row.out)
}

#get proportion of time spent in each cell by period and fish
fish.gridtime <- function(df, f = 1, per = "Pre"){
  #creates a list of transition matrices (all 5 periods) for a specified fish number
  df <- df %>% filter(fish == f)
  #divide observations by period2
  pre <- df %>% filter(Period2 == "PreCO2") %>% group_by(grid) %>% mutate(sumtime = sum(dt)) %>% arrange(grid)
  inc <- df %>% filter(Period2 == "IncreasingCO2") %>% group_by(grid) %>% mutate(sumtime = sum(dt)) %>% arrange(grid)
  dur <- df %>% filter(Period2 == "DuringCO2") %>% group_by(grid) %>% mutate(sumtime = sum(dt)) %>% arrange(grid)
  dec <- df %>% filter(Period2 == "DecreasingCO2") %>% group_by(grid) %>% mutate(sumtime = sum(dt)) %>% arrange(grid)
  post <- df %>% filter(Period2 == "PostCO2") %>% group_by(grid) %>% mutate(sumtime = sum(dt)) %>% arrange(grid)
  #get just the sum time spent in each cell as a vector arranged by grid cell name
  v1 <- pre %>% dplyr::select(sumtime) %>% distinct(sumtime) %>% pull(sumtime)
  v2 <- inc %>% dplyr::select(sumtime) %>% distinct(sumtime) %>% pull(sumtime)
  v3 <- dur %>% dplyr::select(sumtime) %>% distinct(sumtime) %>% pull(sumtime)
  v4 <- dec %>% dplyr::select(sumtime) %>% distinct(sumtime) %>% pull(sumtime)
  v5 <- post %>% dplyr::select(sumtime) %>% distinct(sumtime) %>% pull(sumtime) %>% na.omit()
  print(v5)
  #calculate proportion of time 
  v1.sum <- sum(v1)
  v2.sum <- sum(v2)
  v3.sum <- sum(v3)
  v4.sum <- sum(v4)
  v5.sum <- sum(v5)
  print(v5.sum)
  v1 <- v1 / v1.sum
  v2 <- v2 / v2.sum
  v3 <- v3 / v3.sum
  v4 <- v4 / v4.sum
  v5 <- v5 / v5.sum
  #return list of vectors of proportion of times for each grid cell
  list.m <- list("Pre"=v1, "Inc"=v2, "Dur"=v3, "Dec"=v4, "Post"=v5)
  return(list.m[[per]])
}



#MORANS I TESTING----------------------------------------------------------------------

mvmtavg <- function(df){
  #creates new variables for averages of movement descriptors by grid cell
  output <- df
  output <- output %>% group_by(grid) %>% mutate(avgdist = mean(dist))
  output <- output %>% group_by(grid) %>% mutate(avgt = mean(dt))
  output <- output %>% group_by(grid) %>% mutate(avgdispl = mean(R2n))
  output <- output %>% group_by(grid) %>% mutate(avgrang = mean(rel.angle))
  output <- output %>% group_by(grid) %>% mutate(avgaang = mean(abs.angle))
  output <- output %>% group_by(grid) %>% mutate(vel = dist/(dt))
  output <- output %>% group_by(grid) %>% mutate(avgv = mean(vel))
  output <- output %>% dplyr::select("grid", "x.zones", "y.zones", "avgdist", "avgt", "avgdispl", "avgrang", "avgaang", "avgv")
  output <- output %>% distinct(grid, .keep_all = TRUE)
  output <- output %>% arrange(grid)
  output$n <- seq.int(nrow(output))
  emp_coord <- lapply(output$grid, as.character)
  emp_coord <- unlist(emp_coord)
  emp_coord <- convert2coord(emp_coord)
  emp_coord$n <- seq.int(nrow(emp_coord))
  output <- left_join(output, emp_coord, by = "n")
  return(output)
}

library(ape)
pre <- carp1mchain2(200000, carp1mPre2, carp1m2Pre2)
pre <- convert2coord(pre)
pre <- pre %>% add_count(x,y)
pre <- pre[!duplicated(pre),]
pre <- pre %>% arrange(y) %>% arrange(x)
pre.dist <- as.matrix(dist(cbind(pre$x, pre$y)))
pre.dist <- 1/pre.dist
diag(pre.dist) <- 0

#get averages of movement metrics for each grid cell
#pre co2 period
pre.emp <- mvmtavg(carp1cPre2)
#Moran's I a different way, bootstrapping (monte carlo)
library(spdep)
w <- cell2nb(nrow = northing.zones, ncol = easting.zones, type="queen", torus=FALSE)
ww <- nb2listw(w, style= "U")
moran.mc(pre.emp$avgdist, ww, nsim=99) #pvalue=0.1
moran.mc(pre.emp$avgdispl, ww, nsim=99) #pvalue=0.01
moran.mc(pre.emp$avgrang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.28
moran.mc(pre.emp$avgaang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.94
moran.mc(pre.emp$avgt, ww, nsim=99) #pvalue=0.04
moran.mc(pre.emp$avgv, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.06

#During CO2 period
#get averages of movement metrics for each grid cell
dur.emp <- mvmtavg(carp1cDur2)
#monte carlo morans i test
moran.mc(dur.emp$avgdist, ww, nsim=99) #pvalue = 0.01, stat = 0.36568
moran.mc(dur.emp$avgdispl, ww, nsim=99) #pvalue=0.01, stat = 0.24934
moran.mc(dur.emp$avgrang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.46
moran.mc(dur.emp$avgaang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.96
moran.mc(dur.emp$avgt, ww, nsim=99) #pvalue=0.01, stat = 0.14679
moran.mc(dur.emp$avgv, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.01, stat=0.36796

#Post CO2 period
#get averages of movement metrics for each grid cell
post.emp <- mvmtavg(carp1cPost2)
#monte carlo morans i test
moran.mc(post.emp$avgdist, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue = 0.01, stat = 0.2485
moran.mc(post.emp$avgdispl, ww, nsim=99) #pvalue=0.01, stat = 0.25339
moran.mc(post.emp$avgrang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.67
moran.mc(post.emp$avgaang, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.01, stat=0.47991 (A LOT REMOVED THOUGH)
moran.mc(post.emp$avgt, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.08
moran.mc(post.emp$avgv, ww, nsim=99, na.action = na.omit, zero.policy = TRUE) #pvalue=0.01, stat=0.2601