#summarize lamprey by treatmentTrial ? with basic trajectory calculations
#make df from carp1
lampreydf <- ldply(lamprey_ltraj)
#make order of observations in order to join with infolocs
lampreydf <- mutate(lampreydf, num = 1:62545)

#add infolocs to lampreydf
lampreyinfo <- ldply(infolocs(lamprey_ltraj))
#make order of observations in order to join with lampreydf
lampreyinfo <- mutate(lampreyinfo, num = 1:62545)
lampreydf <- left_join(lampreydf, lampreyinfo, by = "num")

#separate fish 1 from lampreydf
lamprey.1 <- lampreydf[1:3137,]
lamprey.2 <- lampreydf[1:3603,]

#subset lamprey.1 by treatmentTrial
lamp_cont <- filter(lamprey.1, treatmentTrial == "control")
lamp_l.nic <- filter(lamprey.1, treatmentTrial == "nic" & treatSide == "Right control")
lamp_r.nic <- filter(lamprey.1, treatmentTrial == "nic" & treatSide == "Left control")
lamp_l.tfm <- filter(lamprey.1, treatmentTrial == "tfm" & treatSide == "Right control")
lamp_r.tfm <- filter(lamprey.1, treatmentTrial == "tfm" & treatSide == "Left control")
lamp_both <- filter(lamprey.1, treatmentTrial == "tfmnic")

#create boxplot displaying movement segment distance Lamprey
lamprey1.1mvmt <- ggplot(lampreydf, aes(x = treatmentTrial, y = dist, fill = treatmentTrial)) + 
  geom_boxplot() +
  scale_fill_manual("treatmentTrial", values = c("sienna1", "sienna3", "sienna4","orange")) +
  theme_bw() +
  ylab(expression("Log"[10]*" Movement Segment Distance"))  +
  xlab("Period of Lampricide Treatment") +
  scale_y_log10() +
  ggtitle("Distance of Segment Movement by treatmentTrial") +
  theme(
    strip.background = element_blank()
  ) 

#create boxplot displaying relative angle (absolute values)
lampreyARA <- mutate(lampreydf, absrelang = abs(rel.angle))

lamprey1ang <- ggplot(lampreyARA, aes(x = treatmentTrial, y = absrelang, fill = treatmentTrial)) + 
  geom_boxplot() + 
  scale_fill_manual("treatmentTrial", values = c("sienna1", "sienna3", "sienna4","orange")) +
  theme_bw() +
  ylab("Absolute Value of Relative Angle") +
  xlab("Period of Lampricide Treatment") +
  ggtitle("Relative Angle Movement by treatmentTrials") +
  theme(
    strip.background = element_blank()
  ) 

lamprey1.1mvmt <- ggplot(lamprey.2, aes(x = treatSide, y = dist, fill = treatSide)) + 
  geom_boxplot() +
  scale_fill_manual("treatSide", values = c("sienna1", "sienna3", "sienna4","orange","red","yellow")) +
  theme_bw() +
  ylab(expression("Log"[10]*" Movement Segment Distance"))  +
  xlab("Period of Lampricide Treatment") +
  scale_y_log10() +
  ggtitle("Distance of Segment Movement by treatmentTrial") +
  theme( strip.background = element_blank()
  )

Lamp_dist <- ggplot(lamprey.2, aes(x = treatSide, y = dist, color = side, 
                                   shape = treatment,
                                   fill = treatment)) + 
  geom_boxplot() +
  geom_point() + 
  facet_grid( ~ treatmentTrial, scales = "free_x") + 
  scale_color_manual("Side", values = c("black", "blue")) +
  scale_fill_manual("Treatment", values = c("red", "orange", "yellow")) +
  scale_shape("Treatment") + 
  theme_bw() +
  scale_y_log10() +
  ylab(expression("Log"[10]*" Movement Segment Distance")) +
  xlab("Treatment/side comination") +
  ggtitle("Distance of Segment Movement by Treatment and Sides") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    strip.background = element_blank()
  )

#summarize lamprey by treatmentTrial ? with basic trajectory calculations
#make df from carp1
lampreydf <- ldply(lamprey_ltraj)
#make order of observations in order to join with infolocs
lampreydf <- mutate(lampreydf, num = 1:62545)

#add infolocs to lampreydf
lampreyinfo <- ldply(infolocs(lamprey_ltraj))
#make order of observations in order to join with lampreydf
lampreyinfo <- mutate(lampreyinfo, num = 1:62545)
lampreydf <- left_join(lampreydf, lampreyinfo, by = "num")
lamprey.2 <- lampreydf[1:3603,]
lamprey.2ARA <- mutate(lamprey.2, absrelang = abs(rel.angle))

Lamp_ARA <- ggplot(lamprey.2ARA, aes(x = treatSide, y = absrelang, color = side, 
                                     shape = treatment,
                                     fill = treatment)) + 
  geom_boxplot() +
  geom_point() + 
  facet_grid( ~ treatmentTrial, scales = "free_x") + 
  scale_color_manual("Side", values = c("black", "blue")) +
  scale_fill_manual("Treatment", values = c("red", "orange", "yellow")) +
  scale_shape("Treatment") + 
  theme_bw() +
  ylab("Absolute Value of Relative Angle") +
  xlab("Treatment/side combination") +
  ggtitle("Relative Angle Movement by Treatment Sides") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    strip.background = element_blank()
  )
