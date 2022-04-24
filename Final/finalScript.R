library(survival)
library(ggplot2)
library(survminer)

br_surv <- read.table("br_surv.txt", header = T, na.strings = ".")
ncol(br_surv)

# NA represent missing data, check missing data
# need to check which variable is significant
sum(is.na(br_surv)) / nrow(br_surv)
colSums(is.na(br_surv))

time <- br_surv$survival
status <- br_surv$dead

# KM estimate
fit<-survfit(Surv(time, status)~ age, conf.type="log-log", data = br_surv)
plot(fit, xlab = "Time", ylab = "Estimated Survival", 
     mark.time = T, conf.int = T, main="kM plot")

fit.weib <- survreg(Surv(time, status)~.-ID -disfree+ log(disfree), data = br_surv)
sig1<-fit.weib$scale
mu1<-fit.weib$coef
summary(fit.weib)

fit_arm<-survfit(Surv(time, status)~perform, conf.type="log-log", data = br_surv)
temp <- length(unique(br_surv$perform))
col <- c("#E7B800", "#2E9FDF", "#04e700", "#e700d0")
name = "perform"
ggsurvplot(fit_arm,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = col[1:temp])
