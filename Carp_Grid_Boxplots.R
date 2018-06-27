#Boxplot for absolute angles by period in zone's low and high
ggplot(carp1c, aes(x=Period, y=abs.angle, fill = Period)) + geom_boxplot() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("orange", "yellow", "red")) +
  scale_shape("Period") + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Absolute Angles by Period in High/low Co2 Zones")

#Boxplot for distance traveled by period in zone's low and high
ggplot(carp1c, aes(x=Period, y=dist, fill=Period)) + geom_boxplot() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("orange", "yellow", "red")) +
  scale_shape("Period") + 
  scale_y_log10() +
  ylab(expression("Log"[10]*"Distance"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distance by Period in High/Low Zones") 

#Boxplot for R2n by period in zone's low and high
ggplot(carp1c, aes(x=Period, y=R2n, fill= Period)) + geom_boxplot() +
  facet_grid( ~ zone, scales = "free_x") + 
  scale_fill_manual("Period", values = c("orange", "yellow", "red")) +
  scale_shape("Period") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("R2n by Period in High/Low Co2 Zones")
