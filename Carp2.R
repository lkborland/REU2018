# sets working directory in my E drive for asian carp trial data
setwd("E:/REU 2018/Archive")

# setting data frames for each of three trials, reading from csv files
trial1 <- read.csv("trial1data.csv")
trial2 <- read.csv("trial2data.csv")
trial3 <- read.csv("trial3data.csv")


# make animation showing trial 1 fish 1, first 10 points
library(plyr)
library(dplyr)
library(ggmap)
library(ggplot2)
library(gganimate)
library(magick)
library(animation)
library(devtools)
library(chron)

#make box "map" to show study pool
boxes <- data.frame(read.table(text = "0 0 10 0 
                               0 0 0 5 
                               10 0 0 5
                               10 5 -10 0 
                               10 2.25 -6.9 0
                               10 2.5 -6.9 0
                               3.1 2.25 0 .25
                               4 0 0 2.25", header = F))
boxes <- setNames(boxes, c("Northing", "Easting", "delta_lat", "delta_long"))
ggplot() + 
  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat))


#make animation for trial1, day1, fish1
# first try (DOES NOT WORK THIS WAY)
#head(boxes)
#boxes$Period = NA
#boxes$Time = NA
#p <- ggplot(data = fish1.1.1pre[1:250,], aes(Easting, Northing, color = Period, frame = Time)) + 
#  geom_point(alpha = 0.7) + 
#  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat))
#p
#gganimate(p, title_frame = T, interval = 0.2, filename = "animate1.gif")

#attempt 2 USE THIS ONE
#set theme for animation to plain black and white
theme_set(theme_bw())

#create variable for fish1 trial1
fish1.1 <- trial1 %>% 
  filter(TagCode == 2059.535)

#set boxes Period and Time = NA so they work in the ggplot
boxes$Period = NA
boxes$Time = NA

#create plot with ggplot, geom_point for telemetry data and geom_segment for pool location
p <- ggplot(data = fish1.1[1:500,]) + 
  geom_point(aes(x=Easting, y=Northing, frame = Time, color=Period), alpha = 0.7) + 
  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat))
p
#animate the plot p by time interval 0.2, save as a gif called animate2: 
#500 obs takes about 2 min on my laptop
gganimate(p, title_frame = T, interval = 0.2, filename = "animate2.gif")

#make correct order of periods
order1.1 <- c("PreCO2", "IncreasingCO2", "DuringCO2", "DecreasingCO2", "PostCO2")
fish1.1 <- fish1.1 %>%
            mutate(Period2 = factor(Period2, levels=order1.1)) %>%
            arrange()

#separate points by period for fish 1 trial 1
library(RColorBrewer)
head(boxes)
boxes$Period = NA
boxes$Time = NA
p1.1 <- ggplot(data = fish1.1) + 
  geom_point(aes(Easting, Northing, color = Period2), alpha = 0.7) + 
  scale_color_manual(values=c("sienna1", "sienna2", "sienna3",  "sienna4", "saddlebrown"), guide=FALSE) +
  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat)) +
  facet_wrap(~ Period2, ncol=5)
p1.1


#create variable for fish2 trial1
fish1.2 <- trial1 %>% 
  filter(TagCode == 2325.423)

#make correct order of periods
fish1.2 <- fish1.2 %>%
  mutate(Period2 = factor(Period2, levels=order1.1)) %>%
  arrange()

#separate points by period for fish 2 trial 1
p1.2 <- ggplot(data = fish1.2) + 
  geom_point(aes(Easting, Northing, color = Period2), alpha = 0.7) + 
  scale_color_manual(values=c("palegreen1", "palegreen2", "palegreen3",  "seagreen4", "palegreen4"), guide=FALSE) +
  geom_segment(data=boxes, aes(x=Easting, y=Northing, xend = Easting + delta_long, yend = Northing + delta_lat)) +
  facet_wrap(~ Period2, ncol=5)
p1.2


