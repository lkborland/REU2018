## Load packages
library(ggplot2)
library(data.table)
library(lmerTest)
library(e1071)

## define colorblind friendly colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Load data
trialLocation <- fread("trialLocation.csv")
acceleration <- fread("acceleration.csv")


## Set factor orders
trialLocation[ , Trial := factor(Trial)]
trialLocation[ , Period := factor(Period, levels = unique(Period))]
trialLocation[ , unique(Period)]

trialLocation[ , Period2 := factor(Period2, levels = unique(Period2))]
trialLocation[ , unique(Period2)]

## Relative angle
rel.angle <- trialLocation[ ,
                           .(median = quantile(rel.angle,
                                               probs = c(0.5), na.rm = TRUE),
                             u75 = quantile(rel.angle,
                                            probs = c(0.75), na.rm = TRUE)),
                           by = .(TagCodeTrial, Species, Period,
                                  Period2, Trial)]

raAllData <- ggplot(data = trialLocation, aes(x = Period2, y = rel.angle, color = Species)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Absolute relative angle") +
    scale_color_manual(values = cbPalette) +
    xlab(expression("Period of CO"[2]*" treatment"))
print(raAllData)

ggsave('raAllData.pdf', raAllData, width = 6, height = 4)

## Create plot
u75lmerRel <- lmer(u75 ~ Period2 + Species + (1|Trial), data = rel.angle)
summary(u75lmerRel)
u75lmerRelCI <- data.frame(cbind(fixef(u75lmerRel),
                                 confint(u75lmerRel)[ - c(1:2), ]))


u75lmerRelCI$Parameter <-
    gsub("\\(|\\)|Period2", "", rownames(u75lmerRelCI))
colnames(u75lmerRelCI)[1:3] <- c("Coefficient", "l95", "u95")
u75lmerRelCI$Endpoint <- "Relative angle"



u75plotRel <- ggplot(data = rel.angle, aes(x = Period2, y = u75, color = Species)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of absolute relative angle") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw()
u75plotRel
ggsave("u75plotRel.pdf", u75plotRel, width = 6, height = 4)


## absolute relative angle
abs.rel.angle <- trialLocation[ ,
                               .(median = quantile(abs.rel.angle,
                                                   probs = c(0.5), na.rm = TRUE),
                                 u75 = quantile(abs.rel.angle,
                                                probs = c(0.75), na.rm = TRUE)),
                               by = .(TagCodeTrial, Species, Period,
                                      Period2, Trial)]
araAllData <- ggplot(data = trialLocation, aes(x = Period2, y = abs.rel.angle, color = Species)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Absolute relative angle") +
    xlab(expression("Period of CO"[2]*" treatment"))
print(araAllData)
ggsave('araAllData.pdf', araAllData, width = 6, height = 4)

## Create plot
u75lmer <- lmer(u75 ~ Period2 + (1|Trial) + Species, data = abs.rel.angle)
summary(u75lmer)
confint(u75lmer)



u75lmerAbsCI <- data.frame(cbind(fixef(u75lmer),
                                 confint(u75lmer)[ - c(1:2), ]))


u75lmerAbsCI$Parameter <-
    gsub("\\(|\\)|Period2", "", rownames(u75lmerAbsCI))
colnames(u75lmerAbsCI)[1:3] <- c("Coefficient", "l95", "u95")
u75lmerAbsCI$Endpoint <- "Absolute relative angle"




u75plot <- ggplot(data = abs.rel.angle, aes(x = Period2, y = u75, color = Species)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of absolute relative angle") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw()
u75plot
ggsave("u75plot.pdf", u75plot, width = 6, height = 4)


## Acceleration
acceleration
acc <- acceleration[ ,
                    .(median = quantile(acc,
                                        probs = c(0.5), na.rm = TRUE),
                      u75 = quantile(acc,
                                     probs = c(0.75), na.rm = TRUE)),
                    by = .(TagCodeTrial, Species, Period,
                           Period2, Trial)]

accAll <- ggplot(acceleration, aes(x = Period2, y = acc, color = Species)) +
    geom_boxplot() +
    ylab(expression("Fish acceleration m"*s^-2)) +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme(
        strip.background = element_blank()
    ) +
    theme_minimal() +
    scale_y_sqrt()
print(accAll)
ggsave("accAll.pdf", accAll, width = 6, height = 4)


acc$Period2 <- factor(acc$Period2)
acc$Period2 <- factor(acc$Period2, levels = levels(acc$Period2)[c(5,1:4)])

acc_u75lmer <- lmer(u75 ~ Period2 + (1|Trial) + Species, data = acc)
summary(acc_u75lmer)

cbind(fixef(acc_u75lmer), confint(acc_u75lmer)[ -c(1:2),])


u75lmerACCCI <- data.frame(cbind(fixef(acc_u75lmer),
                                 confint(acc_u75lmer)[ - c(1:2), ]))


u75lmerACCCI$Parameter <-
    gsub("\\(|\\)|Period2", "", rownames(u75lmerACCCI))
colnames(u75lmerACCCI)[1:3] <- c("Coefficient", "l95", "u95")
u75lmerACCCI$Endpoint <- "Acceleration"


u75lmerAbsCI
u75lmerRelCI
u75lmerACCCI

ACC_u75plot <-
    ggplot(data = acc, aes(x = Period2, y = u75, color = Species)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of acceleration") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw()
ACC_u75plot
ggsave("acc_u75plot.pdf", ACC_u75plot, width = 6, height = 4)


## look at distance traveled
distT <- trialLocation[ ,
                      .(median = quantile(dist,
                                          probs = c(0.5), na.rm = TRUE),
                        u75 = quantile(dist,
                                       probs = c(0.75), na.rm = TRUE)),
                      by = .(TagCodeTrial, Species, Period,
                             Period2, Trial)]

distAllData <- ggplot(data = trialLocation, aes(x = Period2, y = dist, color = Species)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Distance traveled (m)") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    scale_color_manual(values = cbPalette)
distAllData
ggsave('distAllData.jpg', distAllData, width = 6, height = 4)
ggsave('Figure_1.jpg', distAllData, width = 6, height = 4)


## Create plot
distT$Period2
dist_u75lmer <- lmer(u75 ~ Period2 + Species +  (1|Trial), data = distT)
summary(dist_u75lmer)

u75lmerDistCI <- data.frame(
    cbind(
        fixef(dist_u75lmer),
        confint(dist_u75lmer)[ -c(1:2),]
    )
)


u75lmerDistCI$Parameter <-
    gsub("\\(|\\)|Period2", "", rownames(u75lmerDistCI))
colnames(u75lmerDistCI)[1:3] <- c("Coefficient", "l95", "u95")
u75lmerDistCI$Endpoint <- "Distance"



dist_u75plot <- ggplot(data = distT, aes(x = Period2, y = u75, color = Species)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of distance traveled") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw()
dist_u75plot
ggsave("dist_u75plot.pdf", dist_u75plot, width = 6, height = 4)


u75lmer <- rbind(
    u75lmerAbsCI,
    u75lmerRelCI,
    u75lmerACCCI,
    u75lmerDistCI
    )


u75lmer$Parameter <- factor(u75lmer$Parameter)


u75lmer$Parameter <- factor(u75lmer$Parameter,
                            levels =
                                rev(c("Intercept",
                                      "SpeciesSVC",
                                      "IncreasingCO2",
                                      "DuringCO2",
                                      "DecreasingCO2",
                                      "PostCO2")),
                            labels =
                                rev(c("Intercept",
                                      "Silver carp",
                                      "Increasing CO2",
                                      "During CO2",
                                      "Decreasing CO2",
                                      "Post CO2")),
                            )

print(u75lmer[ grepl("Abs", u75lmer$Endpoint), ], digits = 2)
print(u75lmer[ grepl("Rel", u75lmer$Endpoint), ], digits = 2)
print(u75lmer[ grepl("Dist", u75lmer$Endpoint), ], digits = 2)

ggAllLmer <- ggplot(u75lmer, aes(x = Parameter, y = Coefficient, ymin = l95, ymax = u95)) +
    geom_point()+
    geom_hline(yintercept = 0, color = 'red') +
    coord_flip() +
    facet_grid( Endpoint ~ . ) +
    geom_linerange() +
    theme_bw() +
    theme(strip.background = element_blank())


print(ggAllLmer)
ggsave("ggAllLmer.pdf", width = 6, height = 6)
ggsave("ggAllLmer.jpg", width = 6, height = 6)
ggsave("Figure_2.jpg", width = 6, height = 6)

## Examine skewness
distT_skew <- trialLocation[ ,
                      .(median = quantile(dist,
                                          probs = c(0.5), na.rm = TRUE),
                        skew = skewness(dist, na.rm = TRUE),
                        u75 = quantile(dist,
                                       probs = c(0.75), na.rm = TRUE)),
                      by = .(TagCodeTrial, Species, Period,
                             Period2, Trial)]

distT_skew[ , mean(skew), by = .(Species, Period2)]


