#Violin Boxplot for absolute angles by period in zone's low and high
Absang_Viol <- ggplot(fish1.1, aes(x=Period, y=abs.angle, fill = Period)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("green", "green3", "darkgreen")) +
  scale_shape("Period") + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Absolute Angles by Period in High/low Co2 Zones")

#Violin Boxplot for distance traveled by period in zone's low and high
Dist_Viol <- ggplot(fish1.1, aes(x=Period, y=dist, fill=Period)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("green", "green3", "darkgreen")) +
  scale_shape("Period") + 
  scale_y_log10() +
  ylab(expression("Log"[10]*"Distance"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distance by Period in High/Low Zones") 

#Violin Boxplot for R2n by period in zone's low and high
R2n <- ggplot(fish1.1, aes(x=Period, y=R2n, fill= Period)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("green", "green3", "darkgreen")) +
  scale_shape("Period") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("R2n by Period in High/Low Co2 Zones")

#density of high and low zones of CO2
R2n_density <- ggplot(fish1.1, aes(x=x, y=y, fill= Period)) + geom_density() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("green", "green3", "darkgreen"))  +
  scale_shape("Period") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("R2n by Period in High/Low Co2 Zones")

R2n_density2 <- ggplot(carp1c, aes(x=x, y=y)) + geom_density_2d() +
  stat_density2d(aes(x = x, y = y)) + 
  facet_grid( ~ zone, scales = "free_x") + 
  scale_shape("Period") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("R2n by Period in High/Low Co2 Zones")

ggplot(carp1c, aes(x=x, y=R2n)) + geom_density_2d(aes(color = zone))
ggplot(carp1c, aes(x = x, y = R2n)) +
  geom_density2d() + 
  stat_density2d(aes(fill = ..level..), geom = "polygon")

order2.2 <- c("PreCO2", "IncreasingCO2", "DuringCO2", "DecreasingCO2", "PostCO2")
fish2.2 <- carp1c %>%
  mutate(Period2 = factor(Period2, levels=order2.2)) %>%
  arrange()

#Violin Boxplot for absolute angles by period in zone's low and high
Absang_V5 <- ggplot(fish2.2, aes(x=Period2, y=abs.angle, fill = Period2)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("palegreen", "palegreen3", "green", "green3", "darkgreen")) +
  scale_shape("Period") + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Absolute Angles by Period in High/low Co2 Zones")

#Violin Boxplot for distance traveled by period in zone's low and high
Dist_V5 <- ggplot(fish2.2, aes(x=Period2, y=dist, fill=Period2)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period2", values = c("palegreen", "palegreen3", "green", "green3", "darkgreen")) +
  scale_shape("Period2") + 
  scale_y_log10() +
  ylab(expression("Log"[10]*"Distance/Time"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distance by Period in High/Low Zones") 

#Violin Boxplot for R2n by period in zone's low and high
R2nV5<- ggplot(fish2.2, aes(x=Period2, y=R2n, fill= Period2)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period2", values = c("palegreen","palegreen3", "green", "green4", "darkgreen")) +
  scale_shape("Period2") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Net Squared Displacement of Trajectory") +
  ggtitle("Net Squared Displacement by Period in High/Low Co2 Zones")

#Violin Boxplot for Reltive Angles in zones low and high
RelangV5 <- ggplot(fish2.2, aes(x=Period2, y = rel.angle, fill = Period2)) + geom_violin() +
  facet_grid( ~ zone, scales = "free_x") +
  scale_fill_manual("Period2", values = c("palegreen", "palegreen3", "green", "green4", "darkgreen")) +
  scale_shape("Period2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative Angle") + 
  ggtitle("Relative Angle by Period in High/Low Co2 Zones")
