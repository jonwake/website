#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 6.
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# See the file PKall.R for the PK data analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Scottish lip cancer data in Section 6.5
#
source("scotlip.txt")
#
Y <- scotdat$Y
X <- scotdat$X
Z <- X-mean(X)
E <- scotdat$E
x1 <- E*(1-X)
x2 <- X*E
#
# Poisson linear link model
#
pind <- glm( Y ~ x1+x2-1, family=poisson(link=identity))
# The summary command gives the estimate and se for beta0
summary(pind)
# Now use the delta method for beta1/beta0, ie the relative risk
pRRx <- pind$coeff[2]/pind$coeff[1]
pbaserisk <- pind$coeff[1]
pVmat <- vcov(pind)
pgderiv <- matrix(c(-pind$coeff[2]/pind$coeff[1]^2,
                   1/pind$coeff[1]),nrow=1,ncol=2)
pvarRRx <- pgderiv %*% pVmat %*% t(pgderiv)
cat("Non-spatial model\n baseline risk = ",pbaserisk,
    ": relative risk (se) = ",pRRx,sqrt(pvarRRx),"\n")
#
# Fig 6.2
#
pdf("scotxplot.pdf",height=4.5,width=4.5)
plot(X,Y/E,xlab="Proportion exposed",ylab="SMR")
abline(pind$coef[1],pind$coef[2])
dev.off()
#
# Minnesota Lung cancer data example at the end of Section 6.5
#
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i] <- mean(radon[radon$county==i,2])
}
rad.avg[26] <- 0
rad.avg[63] <- 0
x <- rad.avg
x[26] <- NA
x[63] <- NA
ynew <- Obs[is.na(x)==F]
Expnew <- Exp[is.na(x)==F]
xnew <- x[is.na(x)==F]
#
# Fit the null model and the model with radon 
#
poismod0 <- glm(ynew~offset(log(Expnew)),family="poisson")
poismod <- glm(ynew~offset(log(Expnew))+xnew,family="poisson")
#
# Empirically calculate the LR statistic -- could do this by looking at the difference
# in deviances which is available from summary(poismod). The MLE and se can be read off
# from this output also
#
empLLS <- 2*sum( fitted(poismod0)-fitted(poismod)+ynew*log(fitted(poismod)/fitted(poismod0)))
#
# Example of quasi-likelihood at the end of Section 6.6
#
quasimod <- glm(ynew~offset(log(Expnew))+xnew,family=quasipoisson(link="log"))
#
# overdispersion parameter alpha estimate
#
alphahat <- sum((ynew-fitted(quasimod))^2/fitted(quasimod))/quasimod$df.res
quasiLLS <- empLLS/alphahat
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sandwich estimation for PK GLM model -- in PKall.R
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Bayesian impropriety plot in Section 6.8
#
gamma <- seq(-4,10,.01)
logpost <- term <- NULL
for (j in 1:length(gamma)){
  for (i in 1:length(E)){
    term[i] <- Y[i]*log( E[i]*( (1-X[i])+X[i]*exp(gamma[j]) ) )
  }
  sumterm <- -sum(Y)*log(sum( E*((1-X)+X*exp(gamma[j]) ) ))
  logpost[j] <- sum(term) + sumterm
}
loglimit <- -sum(Y)*log(sum(E*(1-X)))
for (i in 1:length(E)){
    loglimit <- loglimit + Y[i]*log(E[i]*(1-X[i]))
}
#
# Fig 6.3
#
pdf("improper.pdf",height=6,width=6)
plot(gamma,logpost-max(logpost),type="l",xlab=expression(gamma),ylab="Log Likelihood")
abline(loglimit-max(logpost),0,lty=2)
pind <- glm( Y ~ x1+x2-1, family=poisson(link=identity))
#
# gamma hat is log(pRRx)
#
pRRx <- pind$coeff[2]/pind$coeff[1]
abline(v=log(pRRx),lty=3)
dev.off()
#
# Bayesian analysis of radon and lung cancer data at the end of Section 6.8.4
#
source("http://www.math.ntnu.no/inla/givemeINLA.R")
library(INLA)
inla.upgrade()
#
# INLA analysis of the radon data
#
raddat <- data.frame(ynew,Expnew,xnew)
inla.fit = inla(ynew~offset(log(Expnew))+xnew, data=raddat,family="poisson",
      control.fixed=list(prec=0,prec.intercept=0))
summary(inla.fit)
inla.expectation(function(x) exp(x), inla.fit$marginals.fixed[[2]])
exp(inla.qmarginal(c(.025,.975), inla.fit$marginals.fixed[[2]]))
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bayesian analysis of PK GLM model  -- in PKall.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Lung cancer and radon example
#
library(MASS)
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i] <- mean(radon[radon$county==i,2])
}
rad.avg[26] <- 0
rad.avg[63] <- 0
x <- rad.avg
x[26] <- NA
x[63] <- NA
poismod <- glm(Obs~offset(log(Exp))+rad.avg,family="poisson")
quasimod <- glm(Obs~offset(log(Exp))+rad.avg,family="quasipoisson")
negbin <- glm.nb(Obs~offset(log(Exp))+x)
fp <- Exp*exp(poismod$coef[1]+poismod$coef[2]*rad.avg)
museq <- seq(min(fp),max(fp),.01)
varseq <- museq*(1+museq/negbin$theta)
#
# Fig 6.7
#
pdf("LungCancResFig1.pdf",height=5,width=5)
par(mfrow=c(1,1))
plot(museq,varseq,type="l",xlab=expression(hat(mu)),
     ylab=expression(hat(mu)(1+hat(mu)/hat(b))))
abline(0,2.76,lty=3)
points(fp,jitter(rep(122000,length(fp)),amount=5000),pch=4)
dev.off()
#
# Simulation to demonstrate use of Poisson residuals
#
poissim <- function(alpha,beta0,beta1,expected,x,plots=seq(1,4)){
   n <- length(expected)
   mu <- expected*exp(beta0+beta1*x)
   theta <- rgamma(n,alpha,alpha)
   y <- rpois(n,mu*theta)
   mod <- glm(y~x+offset(log(expected)),family=poisson)
   resid <- (y-fitted(mod))/sqrt(fitted(mod))
   museq <- seq(min(mu),max(mu),.01)
   varseq <- museq*(1+museq/alpha)
   if (plots[1]==1) {
      plot(museq,varseq,type="l",xlab=expression(mu),
      ylab=expression(mu+mu/theta))
      abline(0,1,lty=2)}
   if (plots[2]==1) plot(x,log(y+.1))
   if (plots[3]==1) { 
      plot(fitted(mod),abs(resid),xlab=expression(mu),
           ylab="| Pearson residuals |")
      lines(lowess(fitted(mod),abs(resid),f=.95))}
   if (plots[4]==1) {
      plot(sqrt(fitted(mod)),abs(resid),xlab=expression(sqrt(mu)),
           ylab="| Pearson residuals |")
      lines(lowess(sqrt(fitted(mod)),abs(resid),f=.95))}
}
#
# Fig 6.8
#
par(mfrow=c(2,2))
pdf("LungCancSimFig1.pdf",height=6,width=6)
poissim(alpha=negbin$theta,beta0=negbin$coeff[1],beta1=negbin$coeff[2],expected=Exp,x=rad.avg,plots=c(0,0,0,1))
dev.off()
pdf("LungCancSimFig2.pdf",height=6,width=6)
poissim(alpha=negbin$theta,beta0=negbin$coeff[1],beta1=negbin$coeff[2],expected=Exp,x=rad.avg,plots=c(0,0,0,1))
dev.off()
pdf("LungCancSimFig3.pdf",height=6,width=6)
poissim(alpha=negbin$theta,beta0=negbin$coeff[1],beta1=negbin$coeff[2],expected=Exp,x=rad.avg,plots=c(0,0,0,1))
dev.off()
pdf("LungCancSimFig4.pdf",height=6,width=6)
poissim(alpha=negbin$theta,beta0=negbin$coeff[1],beta1=negbin$coeff[2],expected=Exp,x=rad.avg,plots=c(0,0,0,1))
dev.off()
#
# Now produce the plot for the original data
#
poisresid <- (Obs-fitted(poismod))/sqrt(fitted(poismod))
par(mfrow=c(1,1))
#
# Fig 6.9
#
pdf("poismodresids.pdf",height=5,width=5)
plot(sqrt(fitted(poismod)),abs(poisresid),xlab=expression(sqrt(mu)),
     ylab="| Pearson residuals |")
lines(lowess(sqrt(fitted(poismod)),abs(poisresid),f=.95))
dev.off()

#
# Test assumptions for lung cancer and radon example in Section 6.9
#
respears <- (ynew-fitted(quasimod))/sqrt(fitted(quasimod)*kappahat)
plot(xnew,respears,ylim=c(-max(abs(respears)),max(abs(respears))),
     ylab="Pearson residuals",xlab="Average Radon")
lines(lowess(xnew,respears,f=0.9),lty=3)
abline(0,0,lty=2)
plot(2*sqrt(fitted(quasimod)),respears,ylab="| Pearson residuals |",xlab="Fitted values")
lines(lowess(2*sqrt(fitted(quasimod)),abs(respears),f=0.95),lty=3)
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code to produce geometry figure in Section 6.15
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
hw <- 4
x <- c(1,2)
y <- c(.2,.7)
xup <- max(x)+.1
yup <- 1
#
# Fig 6.11(a)
#
pdf("nonlingeom1.pdf",height=hw,width=hw)
par(mfrow=c(1,1),pty="s")
plot(x,y,xlim=c(0,xup),ylim=c(0,yup),pch=4)
linmod <- lm(y~x+0)
pretend <- data.frame(x,y)
nonlinmod <- nls(y~exp(theta*x),data=pretend,start=list(theta=1))
nonlinmod2 <- nls(y~theta2^x,data=pretend,start=list(theta2=1))
nonlinmod3 <- nls(y~x/(x+phi),data=pretend,start=list(phi=1)) # MM
nonlinmod4 <- nls(y~x/(x+exp(phi)),data=pretend,start=list(phi=1)) # MM 2
abline(0,linmod$coeff,lty=1)
xseq <- seq(0,xup,.01); fitnls <- xseq/(coef(nonlinmod3)+xseq)
lines(xseq,fitnls,lty=3)
legend("topleft",legend=c("Zero intercept linear","Michaelis-Menten"),
       lty=c(1,2),bty="n")
dev.off()
#
# Fig 6.11(b)
#
pdf("nonlingeom2.pdf",height=hw,width=hw);par(mfrow=c(1,1),pty="s")
betalin <- seq(0,3,.1)
plot(y[1],y[2],xlim=c(0,yup),ylim=c(0,yup),pch=4,xlab=expression(mu[1]),
     ylab=expression(mu[2]))
lines(betalin*x[1],betalin*x[2])
betapts <- seq(linmod$coef-3,linmod$coef+3,.1)
for (i in 1:length(betapts)){
    points(betapts[i]*x[1],betapts[i]*x[2])
}
points(linmod$coeff*x[1],linmod$coeff*x[2],pch=19)
lines(x=c(y[1],linmod$coeff*x[1]),y=c(y[2],linmod$coeff*x[2]),lty=3)
phival <- seq(0,100,.1)
dev.off()
#
# Fig 6.11(c)
#
pdf("nonlingeom3.pdf",height=hw,width=hw);par(mfrow=c(1,1),pty="s")
plot(y[1],y[2],xlim=c(0,yup),ylim=c(0,yup),pch=4,xlab=expression(mu[1]),
    ylab=expression(mu[2]))
lines(x[1]/(phival+x[1]),x[2]/(phival+x[2]))
phipts <- seq(0,1000,1)
for (i in 1:length(betapts)){
    points(x[1]/(phipts[i]+x[1]),x[2]/(phipts[i]+x[2]))
}
lines(x=c(y[1],x[1]/(coef(nonlinmod3)+x[1])),
      y=c(y[2],x[2]/(coef(nonlinmod3)+x[2])),lty=3)
points(x[1]/(coef(nonlinmod3)+x[1]),x[2]/(coef(nonlinmod3)+x[2]),pch=19) 
legend(x=.59,y=1.08,legend=expression(paste(theta,"=0")),bty="n")
legend(x=-.08,y=.15,legend=expression(paste(theta,"=",infinity)),bty="n")
dev.off()
#
# Fig 6.11(d)
#
pdf("nonlingeom4.pdf",height=hw,width=hw);par(mfrow=c(1,1),pty="s")
plot(y[1],y[2],xlim=c(0,yup),ylim=c(0,yup),pch=4,xlab=expression(mu[1]),
   ylab=expression(mu[2]))
lines(x[1]/(phival+x[1]),x[2]/(phival+x[2]))
phi2pts <- seq(-10,10,1)
for (i in 1:length(betapts)){
    points(x[1]/(exp(phi2pts[i])+x[1]),x[2]/(exp(phi2pts[i])+x[2]))
}
lines(x=c(y[1],x[1]/(exp(coef(nonlinmod4))+x[1])),
      y=c(y[2],x[2]/(exp(coef(nonlinmod4))+x[2])),lty=3)
points(x[1]/(exp(coef(nonlinmod4))+x[1]),x[2]/(exp(coef(nonlinmod4))+x[2]),
       pch=19) 
legend(x=.5,y=1.08,legend=expression(paste(phi,"= -",infinity)),bty="n")
legend(x=-.08,y=.15,legend=expression(paste(phi,"=",infinity)),bty="n")
dev.off()



