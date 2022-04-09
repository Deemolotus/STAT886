#2
#a
dat<-read.table('p42.txt',header = TRUE)
library(survival)
delta<-dat$status
x<-dat$time
AG<-dat$AG
wbc<-dat$wbc
#scatter plot for AG
plot(AG[delta==1],log(x[delta==1]),pch=4,xlab="cell type (positive=1,negative=0)",
     ylab="log survival time",main="scatter plot:leukemia patients")
points(AG[delta==0],log(x[delta==0]),pch=1)
legend(0.4,2,c("failure","censoring"),pch=c(4,1))
#scatter plot for wbc
plot(wbc[delta==1],log(x[delta==1]),pch=4,xlab="white blood cell count",
     ylab="log survival time",main="scatter plot:leukemia patients")
points(wbc[delta==0],log(x[delta==0]),pch=1)
legend(70,2,c("failure","censoring"),pch=c(4,1))
#km weib and lognormal
fit.km=survfit(Surv(time,status)~1,data = dat)
fit.weib=survreg(Surv(time,status)~1,data = dat)
fit.log=survreg(Surv(time,status)~1,data=dat,dist='lognormal')
mu=fit.weib$coefficients
mu.log=fit.log$coefficients
sig=fit.weib$scale
#b
fit1<-survreg(Surv(time,status)~wbc+as.factor(AG),data=dat)
summary(fit1)
#since p-value of wbc>0.05, we can't reject null hypothesis
#wbc has no significant effect on model
fit2<-survreg(Surv(time,status)~as.factor(AG),data=dat)
summary(fit2)
mu<-fit2$coefficients[1]+fit2$coefficients[2]*1
sig<-fit2$scale

#c
#model2: C-S residual(modified for censoring)
r<-exp(log(x)-mu/sig)+1-delta
plot(AG,r,pch=16,xlab="white blood cell characteristics positive=1 negative=0"
     ,ylab="Residuals",main="Cox-Snell Residuals: leukemia patients")
lines(seq(0,1),rep(1,2),lty=2)
legend(0.4,4,"mean of exp(1)",lty = 2)

#Model2: original C-S
r1<-exp((log(x)-mu)/sig)
fit.res<-survfit(Surv(r1,delta)~1)
plot(fit.res,fun="cumhaz",conf.int=F,xlim=c(0,2),ylim=c(0,2.5),
     xlab="residuals(R)",ylab="K-M Estimates of H(r)",
     main="Estimated H(r) for Original Cox-Snell Residuals")
lines(0:2,0:2,lty=2)

#d
fit2$coefficients

#3
dat1<-read.table("bladder.txt",header = TRUE)
attach(dat1)
time=rep(0,nrow(dat1))
status=rep(0,nrow(dat1))
for(i in 1:nrow(dat1)){
  if(dat1$recurrence_time[i]==0){
    time[i]=dat1$futime[i]
  }
  else{
    time[i]=dat1$recurrence_time[i]
    status[i]=1
  }
}
dat2<-data.frame(time,status,Group,number,size)
#model1
fit.cox1=coxph(Surv(time,status)~number+size+as.factor(Group),data = dat2)
summary(fit.cox1)
#LR statistics for H0
2*(fit.cox1$loglik[2]-fit.cox1$loglik[1])
#model2
fit.cox2=coxph(Surv(time,status)~number+as.factor(Group),data = dat2)
summary(fit.cox2)
#LR test compare model 1 and model 2
2*(fit.cox1$loglik[2]-fit.cox2$loglik[2])#observed LR statistic
1-pchisq(2*(fit.cox1$loglik[2]-fit.cox2$loglik[2]),3)#p-value
#model3
fit.cox3=coxph(Surv(time,status)~number,data = dat2)
summary(fit.cox3)
#LR test compare model 2 and model 3
2*(fit.cox2$loglik[2]-fit.cox3$loglik[2])#observed LR statistic
1-pchisq(2*(fit.cox2$loglik[2]-fit.cox3$loglik[2]),3)#p-value

#km plot for comparing the treatment groups
fit.km=survfit(Surv(time,status)~as.factor(Group),data=dat2)
plot(fit.km,xlab="Remission Time",ylab="Estimated S(t)",lty = 1:2)
legend(40,0.85,c("group1","group2"),lty=1:2)

#deviance residuals
fit3.resids<-resid(fit.cox3,type="deviance")
PS<-dat2$status
grp<-dat2$Group
par(mfrow=c(2,2))
#Residuals vs covariate
plot(PS,fit3.resids,ylab="Deviance Residuals")
plot(grp,fit3.resids,ylab="Deviance Residuals")
#Residuals vs risk score
plot(predict(fit.cox3),fit3.resids,ylab="Deviance Residuals",xlab="Risk Score")
#Normal Q-Q plot
qqnorm(fit3.resids,ylab="Deviance Residuals",xlab = "N(0,1) Quantiles")
abline(0,1,lty=2)

#4



