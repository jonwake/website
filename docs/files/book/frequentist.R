#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 2
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Lung cancer example
#
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i]<-mean(radon[radon$county==i,2])
}
x <- rad.avg
rad.avg[26] <- 0
rad.avg[63] <- 0
x[26] <- NA
x[63] <- NA
#
# Fig 2.1
#
pdf("lungcancerlogass.pdf",width=4,height=3.5)
plot(log(Obs/Exp)~x,xlab="Average Radon (pCi/liter)",ylab="Log SMR")
lines(lowess(log(Obs/Exp)[is.na(x)==F]~x[is.na(x)==F]))
dev.off()
#
# Lung cancer MLE example in Section 2.4
#
poismod <- glm(Obs~offset(log(Exp))+x,family="poisson",na.action=na.omit)
quasimod <- glm(Obs~offset(log(Exp))+x,family=quasipoisson(link=log),
            na.action=na.omit)
summary(poismod)
sqrt(vcov(poismod)[1,1]) # se beta0
sqrt(vcov(poismod)[2,2]) # se beta1
vcov(poismod)[1,2]/(sqrt(vcov(poismod)[1,1])*sqrt(vcov(poismod)[2,2])) # correl
poismod$coeff[2]-1.96*sqrt(vcov(poismod)[2,2]) # 95%
poismod$coeff[2]+1.96*sqrt(vcov(poismod)[2,2]) # CI for beta1
exp(quasimod$coeff[2]-1.96*sqrt(vcov(poismod)[2,2]))  # 95%
exp(quasimod$coeff[2]+1.96*sqrt(vcov(poismod)[2,2]))  # CI for exp(beta1)
#
# Lung cancer QMLEs in Section 2.5
#
quasimod <- glm(Obs~offset(log(Exp))+x,family=quasipoisson(link="log"))
exp(quasimod$coeff[2]-1.96*sqrt(vcov(quasimod)[2,2]))  # 95%
exp(quasimod$coeff[2]+1.96*sqrt(vcov(quasimod)[2,2]))  # CI for exp(beta1)
#
# Negative binomial model in Section 2.5
#
library(MASS)
newy <- Obs[is.na(x)==F]
newx <- x[is.na(x)==F]
newE <- Exp[is.na(x)==F]
modnegbinmle <- glm.nb(newy~offset(log(newE))+newx)
exp(modnegbinmle$coeff[2]-1.96*sqrt(vcov(modnegbinmle)[2,2]))  # 95%
exp(modnegbinmle$coeff[2]+1.96*sqrt(vcov(modnegbinmle)[2,2]))  # CI for exp(beta1)
#
bhat <- 61.3 # Obtained from summary(modnegbinmle)
#
# Now find the method of moments estimator of b
#
bnew <- 0
for (i in 1:5){	
	fit <- glm(newy~offset(log(newE))+newx,family=negative.binomial(bhat))
	mu <- fit$fitted
	bnew <- 1/(sum(((newy-mu)^2-mu)/mu^2)/(length(newy)-2)) 
	bhat <- bnew
	cat("Iteration ",i,bhat,"\n")
}
#
# Fig 2.2
#
pdf("poisquad.pdf",width=4,height=4)
par(mfrow=c(1,1))
z <- seq(15,3000,.5)
plot(z,z+z*z/bhat,type="n",log="xy",xlab="Log Mean",ylab="Log Variance",axes=F)
box()
axis(1)
axis(2,at=c(5,50,500,5000,50000),labels=c(5,50,500,5000,50000))
lines(z,2.81*z)
lines(z,z+z*z/bhat,lty=2)
points(newy,jitter(rep(21,length(newy))))
legend("topleft", c("Linear", "Quadratic"), lty=1:2,bty="n")
dev.off()
#
# Poisson mean example in Section 2.6
#
# Specify mean mu and overdispersion alpha of data
#
n <- c(5,10,15,20,25,50,100)
mu <- 10
alpha <- c(1,2,3)
S <- 5000 # No of simulations
countm <- countq <- counts <- matrix(0,nrow=length(n),ncol=length(alpha))
for (k in 1:length(n)){
  for (l in 1:length(alpha)){
    for (j in 1:S) {
	if (alpha[l]>1) 
             {b <- 1/(alpha[l]-1);a <- mu*b;theta <- rgamma(n[k],a,b)}
	else 
             {theta <- mu}
	y <- rpois(n[k],theta)
	thetahat <- mean(y)
#	Model-based estimate of variance
	modelv <- thetahat/n[k]
	if ((mu > thetahat - qnorm(0.975)*sqrt(modelv) )&&
           (mu < thetahat + qnorm(0.975)*sqrt(modelv) ) ) 
           countm[k,l] <- countm[k,l]+1
#	Quasi-likelihood
	phihat = sum((y-thetahat)**2/thetahat)/(n[k]-1)
	quasiv <- phihat*thetahat/n[k]
	if ((mu > thetahat - qnorm(0.975)*sqrt(quasiv) )&&
           (mu < thetahat + qnorm(0.975)*sqrt(quasiv) ) ) 
           countq[k,l] <- countq[k,l]+1
#	Sandwich
	an <- n[k]/thetahat
	bn <- sum((y-thetahat)**2)/thetahat**2
	sandv <- (1/an)*bn*(1/an)
        if ((mu > thetahat - qnorm(0.975)*sqrt(sandv) )&&
           (mu < thetahat + qnorm(0.975)*sqrt(sandv) ) ) 
           counts[k,l] <- counts[k,l]+1
     }
   }
}
#
# Table 2.3 results
#
for (k in 1:length(n)){
    cat("n=",n[k],round(100*countm[k,1]/S,digits=0),round(100*countq[k,1]/S,digits=0),round(100*counts[k,1]/S,digits=0),round(100*countm[k,2]/S,digits=0),round(100*countq[k,2]/S,digits=0),round(100*counts[k,2]/S,digits=0),round(100*countm[k,3]/S,digits=0),round(100*countq[k,3]/S,digits=0),round(100*counts[k,3]/S,digits=0),"\n")
}
#
# Sandwich estimation for lung cancer in Section 2.6
#
mod1 <- glm(newy~offset(log(newE))+newx,family=poisson)
resid1 <- (newy-mod1$fit)/sqrt(mod1$fit)
# Reproduce regular standard errors.
Dmat <- matrix(0,ncol=2,nrow=length(newy))
Dmat <- cbind( mod1$fit, mod1$fit*newx)
V <- diag(mod1$fit)
invinf <- solve( t(Dmat) %*% solve(V) %*% Dmat )
sqrt(invinf[1,1])
sqrt(invinf[2,2])
invinf[1,2]/sqrt(invinf[1,1]*invinf[2,2])
#  Sandwich standard errors.
ingred <- diag(resid1^2/mod1$fit)
lambhat <- matrix(0,ncol=2,nrow=2)
lambhat <-  t(Dmat) %*% ingred %*% Dmat 
sandcov <- invinf  %*% lambhat %*% invinf
sqrt(sandcov[1,1])
sqrt(sandcov[2,2])
sandcov[1,2]/sqrt(sandcov[1,1]*sandcov[2,2])
sqrt(sandcov[1,1]/vcov(poismod)[1,1]) # 
sqrt(sandcov[2,2]/vcov(poismod)[2,2]) #
sqrt(sandcov[1,1]/vcov(quasimod)[1,1]) # 
sqrt(sandcov[2,2]/vcov(quasimod)[2,2]) #
exp(mod1$coeff[2]-1.96*sqrt(sandcov[2,2]))  # 95%
exp(mod1$coeff[2]+1.96*sqrt(sandcov[2,2]))  # CI for exp(beta1)
#
# Bootstrap for lung cancer in Section 2.7, to give Table 2.4
#
#	Non-parametric bootstrap
#
B <- 1000
n <- length(newy)
yb <- rep(0,n)
betab <- rep(0,B)
yb <- rep(0,n)
for (i in 1:B){	
	indices <- sample(1:n,size=n,replace=TRUE)
	yb <- newy[indices]
        xb <- newx[indices]
	Eb <- newE[indices]
	mod <- glm(yb~offset(log(Eb))+xb,family=poisson)
	betab[i] <- mod$coef[2]
}
quasimod <- glm(newy~offset(log(newE))+newx,family=quasipoisson)
#
# Fig 2.3
#
pdf("lungboot.pdf",width=4,height=4)
hist(betab,prob=TRUE,xlab=expression(hat(beta)[1]),ylim=c(0,75),main="")
curve(dnorm(x,quasimod$coef[2],sqrt(vcov(quasimod)[2,2])),add=T,lty=2)
curve(dnorm(x,poismod$coef[2],sqrt(vcov(poismod)[2,2])),add=T,lty=1)
dev.off()
#
cat("NON-PARAMETRIC\n")
cat("Mean and se:",mean(betab),sqrt(var(betab)),"\n")
cat("95% Normal CI:",exp(mean(betab)-1.96*sqrt(var(betab))),exp(mean(betab)+1.96*sqrt(var(betab))),"\n")
cat("95% Percentile CI:",exp(quantile(betab,probs=c(.025,.975))),"\n")
pivlo <- exp(2*mean(betab)-quantile(betab,probs=c(.975)))
pivhi <- exp(2*mean(betab)-quantile(betab,probs=c(.025)))
cat("95% Pivotal CI:",pivlo,pivhi,"\n")
#
# Fig 2.4
#
#  12 counts in total, from 20 Poisson experiments.
#
n <- 20
ybar <- 0.6
lambda0 <- 1
lambda <- seq(0.3,1.3,.001)
const <- -n*ybar + n*ybar*log(ybar)
loglik <- -n*lambda + n*ybar*log(lambda) - const
llybar <- -n*ybar + n*ybar*log(ybar) - const
lllambda0 <- -n*lambda0 + n*ybar*log(lambda0) - const
pdf("poissontest.pdf",width=5.5,height=6)
plot( lambda, loglik, type="n", xlab=expression(lambda), ylab="Log Likelihood")
lines( lambda, loglik )
arrows(ybar,-n*1.0 + n*ybar*log(1.0) - const,lambda0,-n*1.0 + n*ybar*log(1.0) - const,length=.1,code=3)
#lines(c(ybar,lambda0),c(-n*1.0 + n*ybar*log(lambda0) - const,-n*1.0 + n*ybar*log(lambda0) - const))
text(ybar+(lambda0-ybar)/2,lllambda0+ .2,"Wald Test",adj=0.5)
arrows(1.0+.05,llybar,1.0+.05,lllambda0,length=.1,code=3)
text(1.08,-.9," LR Test",adj=0)
lines(c(1,1),c(min(loglik),-n*1.0 + n*ybar*log(1.0) - const ))
b <- n*(ybar-lambda0)/lambda0
a <- lllambda0 - b*lambda0
arrows(1.0,lllambda0,0.77,a+b*0.77,length=0.1)
text(.91,-.1,"Score Test",adj=0.5)
lowval <- -4
lines(c(ybar,ybar),c(min(loglik),-n*ybar + n*ybar*log(ybar) - const ))
#par(cex=1.5)
text(ybar+.03,lowval,expression(hat(lambda)),adj = 0)
text(1.0+.03, lowval-.07,expression(lambda[0]),adj = 0)
dev.off()
#
# Test statistics for example in Section 2.9
#
cat("Score and Wald Tests: ",n*(lambda0-ybar)^2/lambda0,"\n")
cat("LR Test: ",2*n*((lambda0-ybar) + ybar*(log(ybar)-log(lambda0))),"\n")
# Reparameterization
cat("Reparametrized Wald Test: ",n*log(ybar)^2/1,"\n")
cat("Reparametrized Score Test: ",n*(ybar-1)^2/1,"\n")