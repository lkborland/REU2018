## Load packages
library(ggplot2)
library(data.table)
library(lmerTest)

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

raAllData <- ggplot(data = trialLocation, aes(x = Period2, y = rel.angle)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Absolute relative angle") +
    xlab(expression("Period of CO"[2]*" treatment")) 
print(raAllData)

ggsave('raAllData.pdf', raAllData, width = 6, height = 4)

## Create plot 
u75lmerRel <- lmer(u75 ~ Period2 + (1|Trial), data = rel.angle)
summary(u75lmerRel)
cbind(fixef(u75lmerRel), 
      confint(u75lmerRel)[ - c(1:2), ])

u75plotRel <- ggplot(data = rel.angle, aes(x = Period2, y = u75)) +
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
araAllData <- ggplot(data = trialLocation, aes(x = Period2, y = abs.rel.angle)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Absolute relative angle") +
    xlab(expression("Period of CO"[2]*" treatment")) 
print(araAllData)
ggsave('araAllData.pdf', araAllData, width = 6, height = 4)

## Create plot 
u75lmer <- lmer(u75 ~ Period2 + (1|Trial), data = abs.rel.angle)
summary(u75lmer)
confint(u75lmer)

u75plot <- ggplot(data = abs.rel.angle, aes(x = Period2, y = u75)) +
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

accAll <- ggplot(acceleration, aes(x = Period2, y = acc)) +
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


acc_u75lmer <- lmer(u75 ~ Period2 + (1|Trial), data = acc)
summary(acc_u75lmer)

cbind(fixef(acc_u75lmer), confint(acc_u75lmer)[ -c(1:2),])

acc_u75plot <- ggplot(data = acc, aes(x = Period2, y = u75)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of acceleration") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw() 
acc_u75plot
ggsave("acc_u75plot.pdf", acc_u75plot, width = 6, height = 4)    


## look at distance traveled 

distT <- trialLocation[ ,
                      .(median = quantile(dist,
                                          probs = c(0.5), na.rm = TRUE),
                        u75 = quantile(dist,
                                       probs = c(0.75), na.rm = TRUE)),
                      by = .(TagCodeTrial, Species, Period,
                             Period2, Trial)]

distAllData <- ggplot(data = trialLocation, aes(x = Period2, y = dist)) +
    geom_boxplot()  + theme_minimal() +
    ylab("Distance traveled") +
    xlab(expression("Period of CO"[2]*" treatment"))
distAllData
ggsave('distAllData.pdf', distAllData, width = 6, height = 4)

## Create plot 
dist_u75lmer <- lmer(u75 ~ Period2 + (1|Trial), data = distT)
summary(dist_u75lmer)

cbind(
    fixef(dist_u75lmer),
    confint(dist_u75lmer)[ -c(1:2),]
 )

dist_u75plot <- ggplot(data = distT, aes(x = Period2, y = u75)) +
    geom_violin(draw_quantiles = c(0.5), fill = NA) +
    scale_color_manual(values = cbPalette) +
    ylab("Upper 75th quantile\n of distance traveled") +
    xlab(expression("Period of CO"[2]*" treatment")) +
    theme_bw() 
dist_u75plot
ggsave("dist_u75plot.pdf", dist_u75plot, width = 6, height = 4)    

