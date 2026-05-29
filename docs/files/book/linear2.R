#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 8
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Dental data motivating plots and summaries
#
library(nlme)
data(Orthodont)
dist <- Orthodont$distance[Orthodont$Sex=="Female"]
gdistance <- matrix(Orthodont$distance[Orthodont$Sex=="Female"],
       nrow=11,ncol=4,byrow=T)
age <- Orthodont$age[Orthodont$Sex=="Female"]
gage <- matrix(Orthodont$age[Orthodont$Sex=="Female"],nrow=11,ncol=4,byrow=T)
gSubject <- Orthodont$Subject[Orthodont$Sex=="Female"]
#
# Fig 8.1(a)
#
pdf("linear2dent1a.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(gage[1,],gdistance[1,],ylim=c(min(gdistance),max(gdistance)),
     type="n",ylab="Length (mm)",xlab="Age (years)")
modmarg <- lm(dist~age)
for (i in 1:11){
        text(gage[i,],gdistance[i,],labels=i)
        lines(gage[i,],gdistance[i,])
        modcond <- lm(gdistance[i,]~gage[i,])
}
dev.off()
#
# Fig 8.1(b)
#
pdf("linear2dent1b.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(gage[1,],gdistance[1,],ylim=c(min(gdistance),max(gdistance)),
     type="n",ylab="Length (mm)",xlab="Age (years)")
modmarg <- lm(dist~age)
abline(modmarg$coef[1],modmarg$coef[2],lty=1)
for (i in 1:11){
        modcond <- lm(gdistance[i,]~gage[i,])
        abline(modcond$coef[1],modcond$coef[2],lty=2)
}
dev.off()
#
# Residual variances and correlations in (8.3)
#
residmarg <- matrix(0,nrow=11,ncol=4)
for (i in 1:11){
	residmarg[i,] <- gdistance[i,]-modmarg$fit[1:4]
}
for (i in 1:4){
    cat("SD ",i,sqrt(var(residmarg[,i])),"\n")
    for (j in 1:4){
    	if (j>i) cat("Correlation ",i,j,cov(residmarg[,i],residmarg[,j])/(sqrt(var(residmarg[,i]))*sqrt(var(residmarg[,j]))),"\n") 
    }
}
#
# Section 8.5 dental data analyses
#
# Combined fit -seperate intercepts and slopes, correlated random 
# intercepts and slopes
#
data(Orthodont)
remlMF <- lme( distance ~ I(age-11)+Sex+I(age-11):Sex, data = Orthodont,
                 random=~I(age-11) )
summary(remlMF) # p-value to test whether gender slopes needed from here
#
# Fit full and reduced models under MLE to obtain LR test statistic as
# 2*(216.4176-213.903)
#
mlMF <- lme( distance ~ I(age-11)+Sex+I(age-11):Sex, data = Orthodont,
                 method="ML",random=~I(age-11) )
#
mlMFnull <- lme( distance ~ I(age-11)+Sex, data = Orthodont,
                 method="ML",random=~I(age-11) )
2*(mlMF$logLik-mlMFnull$logLik)
# Easy way to do it:
anova(mlMF,mlMFnull)
#
# Implied marginal standard deviations and correlations 
#
Dreml <- matrix(c(1.84^2,0.21*1.84*0.18,0.21*1.84*0.18,0.18^2),nrow=2)
vart <- c(Dreml[1,1]+4*Dreml[2,2]-4*Dreml[1,2] + remlMF$sigma^2,
          Dreml[1,1]+Dreml[2,2]-2*Dreml[1,2] + remlMF$sigma^2,
          Dreml[1,1]+Dreml[2,2]+2*Dreml[1,2] + remlMF$sigma^2,
          Dreml[1,1]+4*Dreml[2,2]+4*Dreml[1,2] + remlMF$sigma^2)
vart1 <- rep(1.82^2+1.39^2,4)
time <- c(-2,-1,1,2)
corr <- c(
(Dreml[1,1]+(time[1]+time[2])*Dreml[1,2]+time[1]*time[2]*Dreml[2,2])/
  sqrt(vart[1]*vart[2]), # (1,2)
(Dreml[1,1]+(time[1]+time[3])*Dreml[1,2]+time[1]*time[3]*Dreml[2,2])/
  sqrt(vart[1]*vart[3]), # (1,3)
(Dreml[1,1]+(time[1]+time[4])*Dreml[1,2]+time[1]*time[4]*Dreml[2,2])/
  sqrt(vart[1]*vart[4]), # (1,4)
(Dreml[1,1]+(time[2]+time[3])*Dreml[1,2]+time[2]*time[3]*Dreml[2,2])/
  sqrt(vart[2]*vart[3]), # (2,3)
(Dreml[1,1]+(time[2]+time[4])*Dreml[1,2]+time[2]*time[4]*Dreml[2,2])/
  sqrt(vart[2]*vart[4]), # (2,4)
(Dreml[1,1]+(time[3]+time[4])*Dreml[1,2]+time[3]*time[4]*Dreml[2,2])/
  sqrt(vart[3]*vart[4])) # (3,4)
#
# Empirical SDs and correlations for girls are above, now for boys
#
distM <- Orthodont$distance[Orthodont$Sex=="Male"]
bdistance <- matrix(Orthodont$distance[Orthodont$Sex=="Male"],
       nrow=16,ncol=4,byrow=T)
ageM <- Orthodont$age[Orthodont$Sex=="Male"]
bage <- matrix(Orthodont$age[Orthodont$Sex=="Male"],nrow=16,ncol=4,byrow=T)
bSubject <- Orthodont$Subject[Orthodont$Sex=="Male"]
bmodmarg <- lm(distM~ageM)
#
bresidmarg <- matrix(0,nrow=16,ncol=4)
for (i in 1:16){
	bresidmarg[i,] <- bdistance[i,]-bmodmarg$fit[1:4]
}
for (i in 1:4){
    cat("SD ",i,sqrt(var(bresidmarg[,i])),"\n")
    for (j in 1:4){
    	if (j>i) cat("Correlation ",i,j,cov(bresidmarg[,i],
   bresidmarg[,j])/(sqrt(var(bresidmarg[,i]))*sqrt(var(bresidmarg[,j]))),"\n") 
    }
}
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# Wishart simulation in Section 8.6
#
# Function to simulate from a Wishart density
#
rwishart <- function(df, p = nrow(SqrtSigma), SqrtSigma = diag(p)) {
  if((Ident <- missing(SqrtSigma)) && missing(p))
    stop("either p or SqrtSigma must be specified")
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
  if(p > 1) {
    pseq <- 1:(p-1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  if(Ident)
    crossprod(Z)
  else
    crossprod(Z %*% SqrtSigma)
}
d <- 4
R <- matrix(c(1.0/d,0,0,0.1/d),nrow=2,ncol=2,byrow=T)
nsim <- 50000
SqrtSigma <- chol(solve(d*R))
rv11 <- rv12 <- rv22 <- irv11 <- irv12 <- irv22 <- rep(0,nsim)
for (i in 1:nsim){
        X <- rwishart(d,2,SqrtSigma)
        rv11[i] <- X[1,1]
        rv12[i] <- X[1,2]
        rv22[i] <- X[2,2]
        temp <- solve(X)
        irv11[i] <- temp[1,1]
        irv12[i] <- temp[1,2]
        irv22[i] <- temp[2,2]
}
library(MASS)
sd1 <- sqrt(irv11)
sd2 <- sqrt(irv22)
rho <- irv12/(sd1*sd2)
sd0min <- 0; sd1min <- 0
sd0max <- 1.5; sd1max <- .5; rhomin <- -1; rhomax <- 1
#
# Fig 8.2(a)
#
pdf("wishrand1.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(density(sd1,adjust=2,n=2048),xlim=c(sd0min,sd0max),
     xlab=expression(sigma[0]),ylab="",main="")
dev.off()
#
# Fig 8.2(b)
#
pdf("wishrand2.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(density(rho,adjust=2.5,from=-1,to=1,kernel="gaussian"),
     xlim=c(rhomin,rhomax),xlab=expression(rho),ylab="",main="")
dev.off()
#
# Fig 8.2(c)
#
pdf("wishrand3.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(density(sd2,adjust=2,n=2048),main="",xlim=c(sd1min,sd1max),
     xlab=expression(sigma[1]),ylab="")
dev.off()
#
# Fig 8.2(d)
#
pdf("wishrand4.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
nval <- 60
f2 <- kde2d(sd1,sd2,h=c(.6,.2),n=nval,lims=c(sd0min,sd0max,sd1min,sd1max))
contour(f2$x,f2$y,f2$z,drawlabels=F,xlab=expression(sigma[0]),
        ylab=expression(sigma[1]))
dev.off()
#
# For dental example with girls data using Bayes in Section 8.6, see
# WinBUGS file DentalWinBUGSGirls.txt
#
# GEE for dental data in Section 8.7, in particular the entries in Table 8.2
#
lmind <- lm(distance~I(age-11)*Sex,data=Orthodont) # Gives model-based ses
                                                   # for independence
modgeeind <- geese(distance~I(age-11)*Sex,id=Subject,data=Orthodont,
          corstr="independence")
summary(modgeeind)
modgee <- geese(distance~I(age-11)*Sex,id=Subject,data=Orthodont,
          corstr="exchangeable")
summary(modgee)
#
# Model-based standard errors:
#
sqrt(modgee$vbeta.naiv[1,1])
sqrt(modgee$vbeta.naiv[2,2])
sqrt(modgee$vbeta.naiv[3,3])
sqrt(modgee$vbeta.naiv[4,4])
#
# Reduced dental dataset example in Section 8.7
#
library(geepack)
zred <- read.table("dental-reduced.txt")
trelldat <- groupedData(Length ~ Age | subred, data=zred)
#
# Fig 8.3
#
pdf("dental_red1.pdf",height=6.5,width=6.5)
plot(trelldat,xlab="Age (years)",ylab="Length (mm)")
dev.off()
#
# GEE
#
# First exchangeable:
#
modgeered <- geese(Length~I(Age-11)*sexred,id=subred,data=zred,
             corstr="exchangeable")
summary(modgeered)
#
# Model is parameterized differently to Table 8.3 so we need to convert
# the estimates and se's
#
cat("Exchangeable Model Results for Reduced Data\n")
cat("beta0 est and se:",modgeered$beta[1]+modgeered$beta[3],
  sqrt(modgeered$vbeta[1,1]+modgeered$vbeta[3,3]+2*modgeered$vbeta[1,3]),"\n")
cat("beta1 est and se:",modgeered$beta[2]+modgeered$beta[4],
  sqrt(modgeered$vbeta[2,2]+modgeered$vbeta[4,4]+2*modgeered$vbeta[2,4]),"\n")
cat("beta2 est and se:",-modgeered$beta[3],sqrt(modgeered$vbeta[3,3]),"\n")
cat("beta3 est and se:",-modgeered$beta[4],sqrt(modgeered$vbeta[4,4]),"\n")
#
# Now independence:
#
modgeered2 <- geese(Length~I(Age-11)*sexred,id=subred,data=zred,
             corstr="independence")
summary(modgeered2)
cat("Independence Model Results for Reduced Data\n")
cat("beta0 est and se:",modgeered2$beta[1]+modgeered2$beta[3],
 sqrt(modgeered2$vbeta[1,1]+modgeered2$vbeta[3,3]+2*modgeered2$vbeta[1,3]),"\n")
cat("beta1 est and se:",modgeered2$beta[2]+modgeered2$beta[4],
 sqrt(modgeered2$vbeta[2,2]+modgeered2$vbeta[4,4]+2*modgeered2$vbeta[2,4]),"\n")
cat("beta2 est and se:",-modgeered2$beta[3],sqrt(modgeered2$vbeta[3,3]),"\n")
cat("beta3 est and se:",-modgeered2$beta[4],sqrt(modgeered2$vbeta[4,4]),"\n")
#
# LMEM: some problems here with convergence, but answers look about right.
#
lmeControl(maxIter=500,msMaxIter=500)
remlMFred1 <- lme(Length ~ I(Age-11)+sexred+I(Age-11):sexred,data = trelldat,
                 random=~1 )
remlMFred <- lme(Length ~ I(Age-11)+sexred+I(Age-11):sexred,data = trelldat,
        random=~I(Age-11),control=list(returnObject=T,tolerance=1e-6,
        niterEM=100,tmaxIter=500,msMaxIter=500,msVerbose=T) )
update(remlMFred1,random=~I(Age-11))
cat("LMEM  Results for Reduced Data\n")
cat("beta0 est and se:",remlMFred$coef$fixed[1]+remlMFred$coef$fixed[3],
 sqrt(remlMFred$varFix[1,1]+remlMFred$varFix[3,3]+2*remlMFred$varFix[1,3]),"\n")
cat("beta1 est and se:",remlMFred$coef$fixed[2]+remlMFred$coef$fixed[4],
 sqrt(remlMFred$varFix[2,2]+remlMFred$varFix[4,4]+2*remlMFred$varFix[2,4]),"\n")
cat("beta2 est and se:",-remlMFred$coef$fixed[3],
     sqrt(remlMFred$varFix[3,3]),"\n")
cat("beta3 est and se:",-remlMFred$coef$fixed[4],
     sqrt(remlMFred$varFix[4,4]),"\n")
#
# Bayes analysis -- see WinBUGS file
#
# Section 8.8 variogram figure
#
times <- seq(0,10,l=10)
Variogram <- seq(0,5,l=10)
#
# Fig 8.4
#
pdf("variogram-example.pdf")
plot(Variogram~times,type="n",axes=F,xlab="Time")
axis(1)
axis(2)
sigmaeps2 <- 1
sigmadel2 <- 4
text(x=9.5,y=sigmaeps2-.2,expression(sigma[epsilon]^2),cex=1.3)
text(x=9.5,y=sigmaeps2+sigmadel2-.2,
     expression(sigma[epsilon]^2+sigma[delta]^2),cex=1.3)
corr <- .3
tseq2 <- seq(0,10,.1)
lines(tseq2,rep(5,length(tseq2)),lty=2)
lines(tseq2,rep(1,length(tseq2)),lty=2)
ymod <- sigmaeps2 + sigmadel2*(1-corr**tseq2)
lines(ymod~tseq2)
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FEV example at the end of Section 8.8
#
#-------------------------------------------------------------------
# Function to evaluate empirical variogram, written by Patrick Heagerty
#
lda.variogram <- function( id, y, x ){
#
# INPUT:  id = (nobs x 1) id vector
#          y = (nobs x 1) response (residual) vector
#          x = (nobs x 1) covariate (time) vector
#
# RETURN:  delta.y = vec( 0.5*(y_ij - y_ik)^2 )
#          delta.x = vec( abs( x_ij - x_ik ) )
#
uid <- unique( id )
m <- length( uid ) 
delta.y <- NULL
delta.x <- NULL
did <- NULL 
for( i in 1:m ){
  yi <- y[ id==uid[i] ]
  xi <- x[ id==uid[i] ]
  n <- length(yi)
  expand.j <- rep( c(1:n), n )
  expand.k <- rep( c(1:n), rep(n,n) )
  keep <- expand.j > expand.k
  if( sum(keep)>0 ){
        expand.j <- expand.j[keep]
        expand.k <- expand.k[keep]
        delta.yi <- 0.5*( yi[expand.j] - yi[expand.k] )^2
        delta.xi <- abs( xi[expand.j] - xi[expand.k] )
        didi <- rep( uid[i], length(delta.yi) )
        delta.y <- c( delta.y, delta.yi )
        delta.x <- c( delta.x, delta.xi )
        did <- c( did, didi )
  }
}
out <- list( id = did, delta.y = delta.y, delta.x = delta.x )
out
}
#
#------------------------------------------------------------------------
# Plotting function
#
myfollowup <- function(id, time, outcome, by = NULL){
    plot(time, outcome, xlab = "Time (years)", ylab = "FEV1",type="n",
         ylim=c(1,5.5))
    id1 <- id
    time1 <- time
    by1 <- by
    outcome1 <- outcome
          id <- id[order(id1, time1, by1)]
            id.factor <- factor(id)
            time <- time[order(id1, time1, by1)]
            outcome <- outcome[order(id1, time1, by1)]
            by <- by[order(id1, time1, by1)]
            by.factor <- factor(by)
            for (i in 1:length(levels(by.factor))) {
                for (j in 1:length(levels(id.factor))) {
                  lines(time[id.factor == levels(id.factor)[j] & 
                    by.factor == levels(by.factor)[i]], outcome[id.factor == 
                    levels(id.factor)[j] & by.factor == levels(by.factor)[i]], 
                    lty = i)
                }
            }
}
fev <- read.table("fev.txt",col.names=c("id","smoker","time","fev1"))
library(xtable)
#
# Tabs 8.4 and 8.5 were constructed by combining the following two:
#
xtable(tapply(fev$fev1,list(fev$time,fev$smoker),FUN=mean))
xtable(tapply(fev$fev1,list(fev$time,fev$smoker),FUN=length),
       caption="Sample sizes by time and smoking status")
#
attach(fev)
tabsum <- matrix(tapply(fev$fev1,list(fev$time,fev$smoker),FUN=mean),
          nrow=7,ncol=2)
#
# Fig 8.5
#
pdf("fev1_fig1.pdf",height=6.5,width=6.5)
plot(c(0,3,6,9,12,15,19),tabsum[,1],ylim=range(tabsum),type="n",
     xlab="Time (years)",ylab="Average FEV1")
lines(c(0,3,6,9,12,15,19),tabsum[,1],lty=1)
lines(c(0,3,6,9,12,15,19),tabsum[,2],lty=2)
legend("topright",legend=c("Former smokers","Current smokers"),
       lty=c(1,2),bty="n",cex=1.5)
dev.off()
#
# Fig 8.6
#
pdf("fev1_fig3.pdf",height=6.5,width=6.5)
myfollowup(id=fev$id,time=fev$time,outcome=fev$fev1,by=fev$smoker)
legend("topright",legend=c("Former smokers","Current smokers"),bty="n",lty=1:2)
dev.off()
#
# Model fits appearing in Tab 8.5
#
# Next two models do not allow for dependence of responses on the same
# individual -- included for comparison only
#
lmmod <- lm(fev1~time+smoker,data=fev)
lmmod2 <- lm(fev1~time+smoker+smoker:time,data=fev) # Last line in Tab 8.5
#
# Now mixed models
#
library(nlme)
trelldat <- groupedData(fev1~time|id,data=fev)
plot(trelldat) # Not very informative as lots of subjects
#
lmemod <- lme(fev1~time,data=fev,random=~1|id) # 1st line in Tab 8.5
lmemod1 <- lme(fev1~time+smoker,data=fev,random=~1|id) # 2nd line in Tab 8.5
lmemod2 <- lme(fev1~time+smoker+smoker:time,data=fev,random=~1|id) # 3rd line
#
# Fit with MLE to allow LR tests
#
lmemodmle <- lme(fev1~time,data=fev,random=~1|id,method="ML")
lmemod1mle <- lme(fev1~time+smoker,data=fev,random=~1|id,method="ML")
lmemod2mle <- lme(fev1~time+smoker+smoker:time,data=fev,random=~1|id,
              method="ML")
anova(lmemodmle,lmemod1mle)
anova(lmemod1mle,lmemod2mle)
#
# Now GEE model: results in Tab 8.6
#
library(geepack)
geemod <- geese(fev1 ~ time + smoker, id=id, data=fev,
                  corstr="exch")
#
# To get Bayes results use WinBUGS file
#
fevsmoke0 <- fev[fev$smoker==0,]
fevsmoke1 <- fev[fev$smoker==1,]
#
# Fit models separately to current and former smokers to examine random
# effects assumption of common variance
#
lmemodsmoker0 <- lme(fev1~time,data=fevsmoke0,random=~1|id)
lmemodsmoker1 <- lme(fev1~time,data=fevsmoke1,random=~1|id)
#
# Allow for random slopes as well as random intercepts and then compare
# using a LR test
lmemodtime <- lme(fev1~time+smoker,data=fev,random=~1+time|id)
anova(lmemod1,lmemodtime)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Assessment of assumptions
#
resid1 <- resid(lmemod1,level=1)
#
# Fig 8.7(a)
#
pdf("assess1a.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
qqnorm(resid1,main="")
dev.off()
#
# Fig 8.7(b)
#
pdf("assess1b.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(fev$time,resid1,xlab="Time (years)",ylab="Residuals")
lines(lowess(fev$time,resid1))
dev.off()
#
# Fig 8.7(c)
#
pdf("assess1c.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
boxplot(resid1~fev$smoker,xlab="Smoking status",ylab="Residuals")
dev.off()
#
# Fig 8.7(d)
#
pdf("assess1d.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(fitted(lmemod1,level=1),abs(resid1),xlab="Fitted values",
     ylab="| residuals |")
lines(lowess(fitted(lmemod1,level=1),abs(resid1)))
dev.off()
uid <- unique( id )
m <- length( uid ); smokerID <- NULL
for( i in 1:m ){
  smokeri <- smoker[ id==uid[i] ]
  smokerID[i] <- smokeri[1]
}
#
# Asses second stage assumptions.
#
lmfits <- lmList(fev1~time|id,fev)
#
# Fig 8.8(a)
#
pdf("assess2a.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
qqnorm(coef(lmfits)[,1],main="")
dev.off()
#
# Fig 8.8(b)
#
pdf("assess2b.pdf",height=4.5,width=4.5)
qqnorm(coef(lmfits)[,2],main="")
dev.off()
#
# Fig 8.8(c)
#
pdf("assess2c.pdf",height=4.5,width=4.5)
plot(coef(lmfits)[,1],coef(lmfits)[,2],xlab="LS intercepts",ylab="LS slopes")
dev.off()
#
# Fig 8.9(a)
#
pdf("assess4a.pdf",height=4.5,width=4.5)
variiout <- lda.variogram( id, resid1, fev$time )
plot(variiout$delta.x,variiout$delta.y,xlab="Time (years)",ylab="Variogram")
lines(lowess(variiout$delta.x,variiout$delta.y))
dev.off()
#
# Fig 8.9(a)
#
pdf("assess4b.pdf",height=4.5,width=4.5)
plot(variiout$delta.x,variiout$delta.y,xlab="Time (years)",ylab="Variogram",
     ylim=c(0,0.08))
lines(lowess(variiout$delta.x,variiout$delta.y))
dev.off()
#
# Assess the need for AR error terms
#
lmemodAR <- lme(fev1~time+smoker,data=fev,random=~1|id,correlation=corAR1())
anova(lmemod1,lmemodAR)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Demonstrate cohort and longitundinal effects in Section 8.9
#
age1 <- c(40,45,50,55,60,65,70); birth1 <- 1930
age2 <- c(40,45,50,55,60); birth2 <- 1940
age3 <- c(40,45,50); birth3 <- 1950
betaC <- .15; betaL <- -0.5; beta0 <- 40
beta1 <- birth1*betaC
beta2 <- birth2*betaC
beta3 <- birth3*betaC
#
# Fig 8.10
#
pdf("crosslong1.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
# Now plot versus calendar time
calendar1 <- age1  + (2000-70)
calendar2 <- age2  + (2000-60)
calendar3 <- age3  + (2000-50)
y1 <- beta0 + birth1*betaC + age1*betaL + rnorm(length(age1),0,0.06)
y2 <- beta0 + birth2*betaC + age2*betaL + rnorm(length(age2),0,0.06)
y3 <- beta0 + birth3*betaC + age3*betaL + rnorm(length(age3),0,0.06)
plot(calendar1,y1,ylim=c(min(beta0 + birth1*betaC + age1*betaL)-.2,
max(beta0 + birth3*betaC + age3*betaL)+.2),xlab="Calendar Year",
    ylab="FEV1",type="n")
points(calendar1,y1,pch=4)
modcal1 <- lm(y1~calendar1)
abline(modcal1$coeff,lty=1)
points(calendar2,y2,pch=4)
modcal2 <- lm(y2~calendar2)
abline(modcal2$coeff,lty=2)
points(calendar3,y3,pch=4)
modcal3 <- lm(y3~calendar3)
abline(modcal3$coeff,lty=3)
legend("bottomleft",legend=c("Cohort 1","Cohort 2","Cohort 3"),
       bty="n",lty=1:3)
dev.off()
#
# Fig 8.11
#
pdf("crosslongbook2.pdf",height=4.5,width=4.5)
plot(age1,beta0 + birth1*betaC + age1*betaL,type="n",xlab="Age (years)",
     ylab="FEV1",ylim=c(min(beta0 + birth1*betaC + age1*betaL)-.2,
     max(beta0 + birth3*betaC + age3*betaL)+.2))
lines(age1,beta0 + birth1*betaC + age1*betaL,lty=1)
points(age1,y1,pch=4)
lines(age2,beta0 + birth2*betaC + age2*betaL,lty=2)
points(age2,y2,pch=4)
lines(age3,beta0 + birth3*betaC + age3*betaL,lty=3)
points(age3,y3,pch=4)
ycross <- c(y1[length(age1)],y2[length(age2)],y3[length(age3)])
xcross <- c(age1[length(age1)],age2[length(age2)],age3[length(age3)])
points(xcross,ycross,pch=4)
modcross <- lm(ycross~xcross)
abline(modcross$coeff,lty=5,lwd=1.3)
legend("bottomleft",legend=c("Cohort 1","Cohort 2","Cohort 3"),bty="n",lty=1:3)
dev.off()






