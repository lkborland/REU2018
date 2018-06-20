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

#subset lamprey.1 by treatmentTrial
lamp_cont <- filter(lamprey.1, treatmentTrial == "control")
lamp_l.nic <- filter(lamprey.1, treatmentTrial == "nic" & treatSide == "Right control")
lamp_r.nic <- filter(lamprey.1, treatmentTrial == "nic" & treatSide == "Left control")
lamp_l.tfm <- filter(lamprey.1, treatmentTrial == "tfm" & treatSide == "Right control")
lamp_r.tfm <- filter(lamprey.1, treatmentTrial == "tfm" & treatSide == "Left control")
lamp_both <- filter(lamprey.1, treatmentTrial == "tfmnic")