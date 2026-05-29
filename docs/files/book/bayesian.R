#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 3
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Specify two quantiles and two associated probs and the function
# returns the mu and sigma of the lognormal
#
LogNormalPriorCh <- function(theta1,theta2,prob1,prob2){
  zq1 <- qnorm(prob1)
  zq2 <- qnorm(prob2)
  mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
  sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
  list(mu=mu,sigma=sigma)
}
#
# Lung cancer prior specification in Section 3.4
#
lungradprior0 <- LogNormalPriorCh(1.5,2/3,.975,.025)
cat("2/3,1.5 wp 0.95: ",lungradprior0$mu,lungradprior0$sigma,"\n")
lungradprior1 <- LogNormalPriorCh(1.2,0.8,.975,.025)
cat("0.8,1.2 wp 0.95: ",lungradprior1$mu,lungradprior1$sigma,"\n")
#
# Function to be used with optim: evaluates the squared difference between
# two probs and the theoretical tail areas corresponding to two candidate
# values of a and b (the beta parameters)
#
priorch <- function(x,q1,q2,p1,p2){
        (p1-pbeta(q1,x[1],x[2]))^2 + (p2-pbeta(q2,x[1],x[2]))^2 }
#
# Fig 3.1
#
p1 <- 0.05
p2 <- 0.95
q1 <- 0.1
q2 <- 0.6
opt <- optim(par=c(1,1),fn=priorch,q1=q1,q2=q2,p1=p1,p2=p2)
cat("a1 and a2 are ",opt$par,"\n")
probvals <- seq(0,1,.001)
pdf("priorchbeta.pdf",width=4,height=4)
plot(probvals,dbeta(probvals,shape1=opt$par[1],shape2=opt$par[2]),type="l",
     xlab="Probability, p",ylab="Beta Density")
abline(v=q1,lty=2)
abline(v=q2,lty=2)
dev.off()
#
# Fig 3.2
#
n <- 9
sigma2 <- n
limits <- 3
mu <- seq(-limits,limits,.01)
m1 <- 0
v1 <- 1
w1 <- n*v1/(n*v1+sigma2)
y1 <- w1*sigma2/n + (1-w1)^2*(m1-mu)^2
m2 <- 0
v2 <- 3
w2 <- n*v2/(n*v2+sigma2)
y2 <- w2*sigma2/n + (1-w2)^2*(m2-mu)^2
pdf("msebayes.pdf",height=4,width=4)
plot(mu,y1,type="n",xlab=expression(mu),ylab="Mean Squared Error",
        ylim=c(0,max(y1,y2)))
lines(mu,y1,lty=2)
lines(mu,y2,lty=3)
abline(sigma2/n,0)
dev.off()
#
# Poisson Lognormal integration example in Section 3.7.
# Summaries of these calculations
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This function evaluates log[ theta^r p*(theta|y) ] where p*(.|.) is the 
# unnormalized posterior.
# This function is used by most of the implementation methods described
# below.
#
gfun <- function(theta,y,Expect,a,b,r) {
        const <- Expect^y*exp(-lgamma(y+1))/sqrt(2*pi*b*b)
	if (r==0) gfun <- log(const) -Expect*exp(theta)+y*theta-
                  (theta-a)^2/(2*b*b)
	if (r!=0) {
# Integrals are defined wrt exp( r*log(theta) + log p*(theta|y) ) so
# we must exclude theta<0 pts if we write in this form -- not a problem 
# for these data (only used by Laplace and GH approaches)
	if (theta < 0) gfun <- -10000 
        else gfun <- log(const) + r*log(theta)-Expect*exp(theta)+y*theta-
                        (theta-a)^2/(2*b*b)
	}
	gfun
}
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to evaluate the function that is to be minimized when the 
# Laplace approximation is carried out.
#
laproot <- function(theta,y,Expect,a,b,r) {
	laproot <- r/theta-Expect*exp(theta)+y-(theta-a)/(b*b)
	laproot
}
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Log ratio for normal prior, used by importance sampling routine.
#
lognfun <- function(theta,y,Expect,a,b,r,m,sd) {
	lognfun <- gfun(theta,y=y,Expect=Expect,a=a,b=b,r=r) - 
		   log(dnorm(theta,m,sd))
	lognfun
}
#
# Implementation comparison example, answers appear in Table 3.2
#
# Read in data
y <- 4; Expect <- 0.25
#
# INLA results (not in book, accuracy is about the same as Laplace, perhaps
# due to the small expected number).
# The following code was supplied by Havard Rue.
#
source("http://www.math.ntnu.no/inla/givemeINLA.R")
inla.upgrade()
library(INLA)
a <- 0; stopifnot(a==0)
b <- 1.38
pdat <- list(y = y, E = Expect)
mod.pois <- inla(y ~ 1, 
                 family='poisson',
                 E=E,
                 data=pdat,
                 control.compute=list(mlik=T),
                 control.fixed=list(mean.intercept=a, prec.intercept=1/b^2),
                 control.inla = list(strategy = "laplace"))
summary(mod.pois)
#
# Now approximate the posterior variance
m <- mod.pois$marginals.fixed[[1]]
t1 <- inla.expectation(function(x) x, m)
t2 <- inla.expectation(function(x) x^2, m)
v1 <- t2-t1^2
#
fac = 10
## compute the true posterior for the intercept
xx = seq(t1 - fac*sqrt(v1), t1 + fac*sqrt(v1), length = 10000)
dx = diff(xx)[1]
ff = dpois(y, lambda = Expect * exp(xx)) * dnorm(xx, mean = a,  sd = b)
ff = ff/max(ff)
ff = ff / sum(ff) / dx
## and its properties
mm = inla.smarginal(list(x = xx, y = ff))
tt1 <- inla.expectation(function(x) x, mm)
tt2 <- inla.expectation(function(x) x^2, mm)
vv1 <- tt2-tt1^2
#
cat("*** INLA posterior mean and var: ", t1, v1,"\n")
cat("*** True posterior mean and var: ", tt1, vv1,"\n")
#
# Obtain the result exactly
#
integral = integrate(function(x, y, Expect, a, b) dnorm(x, mean=a, sd = b) * dpois(y, lambda = Expect*exp(x)),
        lower = t1 - fac*sqrt(v1), upper = t1 + fac*sqrt(v1),
        y=y, Expect=Expect, a=a, b=b)
true.mlik = log(integral$value)
#
print(paste("*** true mlik", true.mlik, "mod.pois$mlik", mod.pois$mlik[[1]]))
#
plot(m, type="l")
lines(mm, col="blue")
title("INLA (black) and the true posterior (blue) for the intercept")
# End of Havard code
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LAPLACE
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Laplace approximation for normalizing constant (r=0).
#
thetatilde <- uniroot(laproot,interval=c(log(4),log(16)),y=y,Expect=Expect,
	      a=a,b=b,r=0)$root
g20 <- -Expect*exp(thetatilde)-1/(b*b)
vtilde <- -1/g20
const <- Expect^y*exp(-lgamma(y+1))/sqrt(2*pi*b*b)
lapapprox <- exp(gfun(thetatilde,y=y,Expect=Expect,a=a,b=b,r=0))*
	     sqrt(2*pi*vtilde)
#
# 	Laplace approximation for posterior mean (r=1).
#
theta1tilde <- uniroot(laproot,interval=c(log(4),log(16)),y=y,Expect=Expect,
	    a=a,b=b,r=1)$root
g21 <- -1/theta1tilde^2-Expect*exp(theta1tilde)-1/(b*b)
v1tilde <- -1/g21
lap1approx <- exp(gfun(theta1tilde,y=y,Expect=Expect,a=a,b=b,r=1))*
	   sqrt(2*pi*v1tilde)
#
# 	Laplace approximation for r=2 which can be used to find the
#       posterior variance.
#
theta2tilde <- uniroot(laproot,interval=c(log(4),log(16)),y=y,Expect=Expect,
	    a=a,b=b,r=2)$root
g22 <- -2/theta2tilde^2-Expect*exp(theta2tilde)-1/(b*b)
v2tilde <- -1/g22
lap2approx <- exp(gfun(theta2tilde,y=y,Expect=Expect,a=a,b=b,r=2))*
	     sqrt(2*pi*v2tilde)
cat("Laplace approx for p(y) = ",lapapprox,"\n")
pmean <- lap1approx/lapapprox
cat("Laplace approx for post mean = ",pmean,"\n")
cat("Laplace approx for post var = ",lap2approx/lapapprox-pmean^2,"\n")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++
# Gauss-Hermite rules
#+++++++++++++++++++++++++++++++++++++++++++++++++
#
library(statmod)
m <- 5 # no of numerical intergration points to use, change these for 
       # different entries in table (m=5,10,15,20)
quad <- gauss.quad(m,kind="hermite")
ghest0 <- 0
ghest1 <- 0
ghest2 <- 0
for (i in 1:m){
	theta <- thetatilde + sqrt(2*vtilde)*quad$nodes[i]
	ghest0 <- ghest0 + 
                  quad$weights[i]*(1/dnorm(theta,thetatilde,sqrt(vtilde)))*
                  exp(gfun(theta,y,Expect,a,b,r=0))/sqrt(pi)
	ghest1 <- ghest1 + 
                  quad$weights[i]*(1/dnorm(theta,thetatilde,sqrt(vtilde)))*
                  exp(gfun(theta,y,Expect,a,b,r=1))/sqrt(pi)
	ghest2 <- ghest2 + 
                  quad$weights[i]*(1/dnorm(theta,thetatilde,sqrt(vtilde)))*
                  exp(gfun(theta,y,Expect,a,b,r=2))/sqrt(pi)
}
cat("G-H no of int pts p(y): ",m,ghest0,"\n")
ghmean <- ghest1/ghest0
cat("G-H no of int pts mean: ",m,ghmean,"\n")
cat("G-H no of int pts: ",m,ghest2/ghest0-ghmean^2,"\n")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#	Importance sampling from normal proposal.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
n <- 5000 # pts to use to evaluate integrals
thsampn <- rnorm(n,thetatilde,sqrt(vtilde))
n0vals <- exp( lognfun(thsampn,y=y,Expect=Expect,a=a,b=b,r=0,
          m=thetatilde,sd=sqrt(vtilde)) )
n1vals <- thsampn*exp( lognfun(thsampn,y=y,Expect=Expect,a=a,b=b,r=0,
          m=thetatilde,sd=sqrt(vtilde)) )
n2vals <- thsampn^2*exp( lognfun(thsampn,y=y,Expect=Expect,a=a,b=b,r=0,
          m=thetatilde,sd=sqrt(vtilde)) )
impn0est <- mean(n0vals)
varn0est <- var(n0vals)/n
cat("Normal Importance Sampling\n")
cat("p(y) Imp samp = ", impn0est, ":  95% Interval: ",impn0est-1.96*sqrt(varn0est),impn0est+1.96*sqrt(varn0est),"\n")
impn1est <- mean(n1vals)
varn1est <- var(n1vals)/n
covn0n1est <- cov(n0vals,n1vals)/n
covn0n2est <- cov(n0vals,n2vals)/n
covn1n2est <- cov(n1vals,n2vals)/n
pmeannest <- impn1est/impn0est
pmeannse <- sqrt( varn0est*impn1est^2/impn0est^4 + varn1est/impn0est^2 -   2*covn0n1est*impn1est/impn0est^3 )
cat("Post mean Imp samp = ", pmeannest, ":  95% Interval for r=1: ",pmeannest-1.96*pmeannse,pmeannest+1.96*pmeannse,"\n")
impn2est <- mean(n2vals)
varn2est <- var(n2vals)/n
pvarnest <- impn2est/impn0est - impn1est^2/impn0est^2
# Derviatives for delta method on var of var
deriv0 <- -impn2est/impn0est^2  + 2*impn1est^2/impn0est^3
deriv1 <- -2*impn1est/impn0est^2
deriv2 <- 1/impn0est
#
pvarnse <- sqrt( deriv0^2*varn0est + deriv1^2*varn1est + deriv2^2*varn2est +2*deriv0*deriv1*covn0n1est + 2*deriv0*deriv2*covn0n2est + 2*deriv1*deriv2*covn1n2est )
cat("Post var Imp samp = ", pvarnest, ":  95% Interval for r=2: ",pvarnest-1.96*pvarnse,pvarnest+1.96*pvarnse,"\n")
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Rejection algorithm, sampling from prior.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
thetahat <- log(y/Expect); a <- 0; b <- 1.38
logM <- -Expect*exp(thetahat) + y*thetahat 
accept <- 0
counter <- 0
n <- 5000 # No of samples to generate from the posterior
thetasamp <- rep(0,n)
#
sum1 <- sum2 <- 0
while (accept < n) {
	counter <- counter + 1
	theta <- rnorm(1,a,b)
	u <- runif(1)
	test <- -Expect*exp(theta) + y*theta - logM
        sum1 <- sum1 + exp( -Expect*exp(theta) + y*theta )
        sum2 <- sum2 + exp( -Expect*exp(theta) + y*theta )^2
	if (log(u) < test ) {
		accept <- accept + 1
		thetasamp[accept] <- theta
	}
}
#
# Fig 3.3(a)
#
pdf("chap3poisnormfig1.pdf",height=3,width=4.0)
par(mfrow=c(1,1),mar=c(4,4,2,1)+.1)
hist(thetasamp,main="",xlab=expression(theta),prob=T,nclass=20)
curve(dnorm(x,a,b),lty=1,add=T)
dev.off()
#
# Fig 3.3(b) (Add ,ylim=c(0,.5) to the hist command if you want to see
# the prior not truncated)
#
pdf("chap3poisnormfig2.pdf",height=3,width=4.0)
par(mfrow=c(1,1),mar=c(4,4,2,1)+.1)
hist(exp(thetasamp),main="",xlim=c(0,42),xlab=expression(exp *~( theta )),prob=T,nclass=20)
curve(dlnorm(x,a,b),lty=1,add=T)
dev.off()
const2 <- Expect^y*exp(-lgamma(y+1))
paccept <- n/counter
cat("Acceptance prob = ",paccept,"\n")
varpa <- const2^2*exp(logM)^2*paccept^2*(1-paccept)/n
pyest <- const2*exp(logM)*paccept # negative binomial estimate
rejpy <- const2*sum1/counter # Importance sampling estimate
rejpyvar <- (const2^2*sum2/counter - rejpy^2)/counter # IS var est
cat("Rejection IS p(y) = ",rejpy,"95% Interval: ",rejpy - 1.96*sqrt(rejpyvar),
     rejpy + 1.96*sqrt(rejpyvar),"\n")
cat("Rejection, neg bin p(y) = ",pyest,"95% Interval: ",
     pyest - 1.96*sqrt(varpa),pyest + 1.96*sqrt(varpa),"\n")
rejmeanest <- mean(thetasamp)
rejvarest <- var(thetasamp)
rejmeanse <- sqrt( rejvarest/n )
rejvarse <- sqrt( (mean((thetasamp-rejmeanest)^4) - rejvarest^2)/n )
geyervarse <- sqrt( mean( (thetasamp^2-mean(thetasamp^2)-2*mean(thetasamp)*
              (thetasamp-mean(thetasamp)))^2)/n )
cat("Rejection, post mean = ",rejmeanest,":  95% Interval: ",
     rejmeanest-1.96*rejmeanse,rejmeanest+1.96*rejmeanse,"\n")
cat("Rejection, post var = ",rejvarest,":  95% Interval: ",
     rejvarest-1.96*rejvarse,rejvarest+1.96*rejvarse,"\n")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Metropolis-Hastings algorithm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
n <- 51000 # no if iterations
thetaMH <- NULL
thetahat <- log(y/Expect)
sethetahat <- sqrt( (thetahat/Expect)*(1/thetahat)^2 )
thetaMH[1] <- thetahat
accept <- 0
for (i in 2:n){
       thetaMH[i] <- thetaMH[i-1]
       thetaprop <- rnorm(1,m=thetaMH[i-1],s=sqrt(3)*sethetahat)
       if (log(runif(1)) < gfun(thetaprop,y,Expect,a,b,r=0) -  
             gfun(thetaMH[i-1],y,Expect,a,b,r=0)) {
    		   thetaMH[i] <- thetaprop; accept <- accept+1
        }
}
burnin <- 1000; indcalc <- seq(burnin+1,n)
# Estimate se's using batching
cat("MH:\n")
cat("Acceptance prob = ",accept/n,"\n")
B <- 1000 # length of each batch 
pmeanMH <- mean(thetaMH[indcalc])
pvarMH <- var(thetaMH[indcalc])
K <- length(indcalc)/B # K is the number of batches
batchmean <- batchvar <- NULL
for (i in 1:K){
       batchmean[i] <- mean(thetaMH[seq(burnin+1+(i-1)*K,burnin+i*K)])
       batchvar[i] <- var(thetaMH[seq(burnin+1+(i-1)*K,burnin+i*K)])
}
batchmeanvar <- var(batchmean)
batchvarvar <- var(batchvar)
batchmeanse <- sqrt(batchmeanvar/K)
batchvarse <- sqrt(batchvarvar/K)
par(mfrow=c(1,2))
#plot(thetaMH); abline(v=burnin)
plot(batchmean)
plot(batchvar)
cat("Posterior mean estimate: ",pmeanMH,pmeanMH-1.96*batchmeanse,pmeanMH+1.96*batchmeanse,"\n")
cat("Posterior var estimate: ",pvarMH,pvarMH-1.96*batchvarse,pvarMH+1.96*batchvarse,"\n")
#cat("Posterior variance estimate from batching: ",sum(batchest)/K,"\n")
#
# Lung cancer and radon example in Section 3.8
#
poisloglik <- function(x,y,Expected,beta0,beta1){
	mu <- rep(0,length(y))
	mu <- Expected*exp(beta0+beta1*x)
	poisloglik <- -sum(mu)+sum(y*log(mu))
	return(poisloglik)
}
# Read in data
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
        rad.avg[i] <- mean(radon[radon$county==i,2])
}
x <- rad.avg
rad.avg[26]<-0
rad.avg[63]<-0
x[26] <- NA
x[63] <- NA
newy <- Obs[is.na(x)==F]
newx <- x[is.na(x)==F]
newE <- Exp[is.na(x)==F]
library(MASS)
#	Univariate random walk M-H 
nit <- 10000
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
mod <- glm(newy~offset(log(newE))+newx,family="poisson")
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
sd0 <- sqrt(vcov(mod)[1,1])
sd1 <- sqrt(vcov(mod)[2,2])
c <- .1
prevloglik <- poisloglik(x=newx,y=newy,Expected=newE,beta0[1],beta1[1])
count1 <- count2 <- count3 <- 0
for( i in 2:nit){
	beta0new <- rnorm(1,beta0[i-1],sd0*c)
	beta1new <- rnorm(1,beta1[i-1],sd1*c)
	newloglik <- poisloglik(newx,newy,newE,beta0new,beta1new)
	u <- runif(1)
	if (log(u) < newloglik-prevloglik){
		beta0[i] <- beta0new
		beta1[i] <- beta1new
		prevloglik <-  newloglik
		count1 <- count1 + 1
	}
	else {
		beta0[i] <- beta0[i-1]
		beta1[i] <- beta1[i-1]
	}
}
#
# Now run again with a different proposal (larger c)
#
iteration <- seq(1,1000,1)
gamma0 <- beta0
gamma1 <- beta1
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
mod <- glm(newy~offset(log(newE))+newx,family="poisson")
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
sd0 <- sqrt(vcov(mod)[1,1])
sd1 <- sqrt(vcov(mod)[2,2])
c <- 2
prevloglik <- poisloglik(newx,newy,newE,beta0[1],beta1[1])
for( i in 2:nit){
	beta0new <- rnorm(1,beta0[i-1],sd0*c)
	beta1new <- rnorm(1,beta1[i-1],sd1*c)
	newloglik <- poisloglik(newx,newy,newE,beta0new,beta1new)
	u <- runif(1)
	if (log(u) < newloglik-prevloglik){
		beta0[i] <- beta0new
		beta1[i] <- beta1new
		prevloglik <-  newloglik
		count2 <- count2 + 1
	}
	else {
		beta0[i] <- beta0[i-1]
		beta1[i] <- beta1[i-1]
	}
}
delta0 <- beta0
delta1 <- beta1
#
#	Bivariate random walk M-H
#
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
mod <- glm(newy~offset(log(newE))+newx,family="poisson")
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
vmat <- vcov(mod)
c <- 2
prevloglik <- poisloglik(newx,newy,newE,beta0[1],beta1[1])
for( i in 2:nit){
	temp <- mvrnorm(1,mu=c(beta0[i-1],beta1[i-1]),Sigma=c*c*vmat)
	beta0new <- temp[1]
	beta1new <- temp[2]
	newloglik <- poisloglik(newx,newy,newE,beta0new,beta1new)
	u <- runif(1)
	if (log(u) < newloglik-prevloglik){
		beta0[i] <- beta0new
		beta1[i] <- beta1new
		prevloglik <-  newloglik
		count3 <- count3 + 1
	}
	else {
		beta0[i] <- beta0[i-1]
		beta1[i] <- beta1[i-1]
	}
}
#
# Fig 3.4(a) - different samples so won't look the same
#
pdf("lungradonMHrandom1.pdf",width=4.5,height=4)
par(mfrow=c(1,1))
plot(iteration,gamma0[iteration],ylab=(expression(beta[0])),xlab="Iteration",
  type="n",main="",ylim=c(min(beta0,gamma0,delta0),max(beta0,gamma0,delta0)))
lines(iteration,gamma0[iteration])
cat("gamma0 summary:", quantile(gamma0,c(.025,.5,.975)),"\n")
cat("gamma1 summary:", quantile(gamma1,c(.025,.5,.975)),"\n")
cat("expgamma0 summary:", quantile(exp(gamma0),c(.025,.5,.975)),"\n")
cat("expgamma1 summary:", quantile(exp(gamma1),c(.025,.5,.975)),"\n")
dev.off()
#
# Fig 3.4(b) - different samples so won't look the same
#
pdf("lungradonMHrandom2.pdf",width=4.5,height=4)
plot(iteration,delta0[iteration],ylab=(expression(beta[0])),xlab="Iteration",
  type="n",ylim=c(min(beta0,gamma0,delta0),max(beta0,gamma0,delta0)),main="")
lines(iteration,delta0[iteration])
dev.off()
#
# Fig 3.4(c) - different samples so won't look the same
#
pdf("lungradonMHrandom3.pdf",width=4.5,height=4)
plot(iteration,beta0[iteration],ylab=(expression(beta[0])),xlab="Iteration",
  type="n",ylim=c(min(beta0,gamma0,delta0),max(beta0,gamma0,delta0)),main="")
lines(iteration,beta0[iteration])
dev.off()
#
# Fig 3.4(d) - different samples so won't look the same
#
pdf("lungradonMHrandom4.pdf",width=4.5,height=4)
plot(iteration,gamma1[iteration],ylab=(expression(beta[1])),xlab="Iteration",
  type="n",ylim=c(min(beta1,gamma1,delta1),max(beta1,gamma1,delta1)),main="")
lines(iteration,gamma1[iteration])
dev.off()
#
# Fig 3.4(e) - different samples so won't look the same
#
pdf("lungradonMHrandom5.pdf",width=4.5,height=4)
plot(iteration,delta1[iteration],ylab=(expression(beta[1])),xlab="Iteration",
  type="n",ylim=c(min(beta1,gamma1,delta1),max(beta1,gamma1,delta1)),main="")
lines(iteration,delta1[iteration])
dev.off()
#
# Fig 3.4(f) - different samples so won't look the same
#
pdf("lungradonMHrandom6.pdf",width=4.5,height=4)
plot(iteration,beta1[iteration],ylab=(expression(beta[1])),xlab="Iteration",
  type="n",ylim=c(min(beta1,gamma1,delta1),max(beta1,gamma1,delta1)),main="")
lines(iteration,beta1[iteration])
dev.off()
cat("POISSON MODEL\n")
cat("beta0 summary:", quantile(beta0,c(.025,.5,.975)),"\n")
cat("beta1 summary:", quantile(beta1,c(.025,.5,.975)),"\n")
cat("expbeta0 summary:", quantile(exp(beta0),c(.025,.5,.975)),"\n")
cat("expbeta1 summary:", quantile(exp(beta1),c(.025,.5,.975)),"\n")
cat("Acceptance rates: ", count1/(nit-1),count2/(nit-1),count3/(nit-1),"\n")
#
# Fig 3.5(a) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(gamma0,main="",ylab=expression(Autocorrelation*~beta[0]))
dev.off()
#
# Fig 3.5(b) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(delta0,main="",ylab=expression(Autocorrelation*~beta[0]))
dev.off()
#
# Fig 3.5(c) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(beta0,main="",ylab=expression(Autocorrelation*~beta[0]))
dev.off()
#
# Fig 3.5(d) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(gamma1,main="",ylab=expression(Autocorrelation*~beta[1]))
dev.off()
#
# Fig 3.5(e) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(delta1,main="",ylab=expression(Autocorrelation*~beta[1]))
dev.off()
#
# Fig 3.5(f) - different samples so won't look identical
#
pdf("lungradonMHrandom2fig1.pdf",width=4.5,height=4)
acf(beta1,main="",ylab=expression(Autocorrelation*~beta[1]))
dev.off()
#
# Fig 3.6 - different samples so won't the same
#
pdf("lungradonMHrandom3.pdf",width=4.5,height=4)
par(mfrow=c(3,2),mar=c(4,4,2,1)+.1)
image(kde2d(beta0,beta1,n=50),col=grey(seq(1,0,length=50)),
 xlab=expression(beta[0]),ylab=expression(beta[1]),main="(a)")
box()
ltheta0 <- beta0+beta1*mean(newx)
ltheta1 <- beta1
image(kde2d(ltheta0,ltheta1,n=50),col=grey(seq(1,0,length=50)),main="(d)",
  ylab=(expression(log * ~theta[1])),xlab=(expression(log * ~theta[0])),);
box()
hist(beta0,xlab=expression(beta[0]),main="(b)")
hist(beta1,xlab=expression(beta[1]),main="(c)")
hist(exp(ltheta0),xlab=expression(theta[0]),main="(e)")
hist(exp(ltheta1),xlab=expression(theta[1]),main="(f)")
dev.off()
#
# Negative binomial example in Section 3.8
#
#++++++++++++++++++++++++++++++++
# Variance function of a negative binomial
#
varfun <- function(alpha,mu){
	varfun <- mu+mu^2/alpha
	return(varfun)
}
#
# Fig 3.7(a)
#
pdf("negbin1.pdf",width=4.5,height=3.5)
par(mfrow=c(1,1),mar=c(4,4,2,1)+.1)
mu <- seq(0,200,.1)
alpha <- seq(50,200,10) # alpha is b!!!!
plot(mu,mu+mu*mu/min(alpha),type="n",xlab="Mean",ylab="Variance",main="(a)")
for (i in 1:length(alpha)){
	lines(mu,varfun(alpha[i],mu))
}
text(.87*max(mu),.95*max(mu+mu*mu/min(alpha)),"b=50")
text(.93*max(mu),.68*max(mu+mu*mu/max(alpha)),"b=200")
lines(mu,mu,lty=2)
dev.off()
muval <- 158 # mean of Y
c1 <- 1.5; firstpt <- muval/(c1-1)
c2 <- 5; secondpt <- muval/(c2-1)
#
# The function LogNormalPriorCh() is at the top of the file
#
bprior <- LogNormalPriorCh(secondpt,firstpt,0.5,.95)
cat("5 50 95 points: ",qlnorm(c(0.05,.5,.95),bprior$mu,bprior$sigma),"\n")
x <- seq(.01,400,.01)
#
# Fig 3.7(b)
#
pdf("negbin2.pdf",width=4.5,height=3.5)
par(mfrow=c(1,1),mar=c(4,4,2,1)+.1)
plot(x,dlnorm(x,bprior$mu,bprior$sigma),type="n",xlab="b",ylab="Prior Density",
  main="(b)")
lines(x,dlnorm(x,bprior$mu,bprior$sigma))
dev.off()
# 
# Function to evaluate the unnormalized log posterior when the likelhood
# is negative binomial and the prior on alpha (the scale parameter) is
# lognormal.
#
negbinloglik <- function(x,y,Expected,beta0,beta1,alpha,prmean,prsd){
	lnmean <- prmean
	lnsigma <- prsd
	mu <- rep(0,length(y))
	mu <- Expected*exp(beta0+beta1*x)
	negbinloglik <- sum(lgamma(y+alpha))-length(y)*lgamma(alpha)+
	sum(y*log(mu))+length(y)*alpha*log(alpha)-sum((y+alpha)*log(mu+alpha))
	negbinloglik <- negbinloglik + dlnorm(alpha,meanlog=lnmean,
            sdlog=prsd,log=T) 
	return(negbinloglik)
}
library(MASS)
mod <- glm.nb(newy~offset(log(newE))+newx)
nit <- 10000
beta0 <- rep(0,nit)
beta1 <- rep(0,nit)
alpha <- rep(0,nit)
beta0[1] <- mod$coeff[1]
beta1[1] <- mod$coeff[2]
alpha[1] <- 61.3
sealpha <- 17.2
vmat <- vcov(mod)
count <- 0
c <- sqrt(3)
prevloglik <- negbinloglik(newx,newy,newE,beta0[1],beta1[1],alpha[1],
    prmean=bprior$mu,prsd=bprior$sigma)
for( i in 2:nit){
	temp <- mvrnorm(1,mu=c(beta0[i-1],beta1[i-1]),Sigma=c*c*vmat)
	alphanew <- rnorm(1,alpha[i-1],c*sealpha)
	beta0new <- temp[1]
	beta1new <- temp[2]
	if (alphanew < 0){
		newloglik <- -100000000
	}
	else {
		newloglik <- negbinloglik(newx,newy,newE,beta0new,beta1new,
                alphanew,prmean=bprior$mu,prsd=bprior$sigma)
	}
# cat("It log post: ",newloglik,beta0new,beta1new,alphanew,"\n")
	u <- runif(1)
	if (log(u) < newloglik-prevloglik){
		beta0[i] <- beta0new
		beta1[i] <- beta1new
		alpha[i] <- alphanew
		prevloglik <-  newloglik
		count <- count + 1
	}
	else {
		beta0[i] <- beta0[i-1]
		beta1[i] <- beta1[i-1]
		alpha[i] <- alpha[i-1]
	}
}
iteration <- seq(1,nit,1)
#
# The following traceplots are not in the book but are used as a convergence 
# check
#
par(mfrow=c(1,3))
plot(iteration,beta0[iteration],type="n")
lines(iteration,beta0[iteration])
plot(iteration,beta1[iteration],type="n")
lines(iteration,beta1[iteration])
plot(iteration,alpha[iteration],type="n")
lines(iteration,alpha[iteration])
#
# Numerical summaries
#
cat("NEGATIVE BINOMIAL MODEL\n")
cat("expbeta0 summary:", quantile(exp(beta0),c(.025,.5,.975)),"\n")
cat("expbeta1 summary:", quantile(exp(beta1),c(.025,.5,.975)),"\n")
cat("alpha summary:", quantile(alpha,c(.025,.5,.975)),"\n")
cat("Acceptance rate: ", count/(nit-1),"\n")
cor(beta0,alpha)
cor(beta1,alpha)
#
# Fig 3.8
#
pdf("negbin3.pdf",width=4.5,height=4.5)
par(mfrow=c(2,3),mar=c(4,4,2,1)+.1)
hist(beta0,xlab=expression(beta[0]),main="(a)")
hist(beta1,xlab=expression(beta[0]),main="(b)")
hist(alpha,xlab="b",main="(c)")
image(kde2d(beta0,beta1,n=50),col=grey(seq(1,0,length=50)),
   xlab=expression(beta[0]),ylab=expression(beta[1]),main="(c)")
box()
image(kde2d(beta0,alpha,n=50),col=grey(seq(1,0,length=50)),ylim=c(30,110),
   xlab=expression(beta[0]),ylab="b",main="(d)")
box()
image(kde2d(beta1,alpha,n=50),col=grey(seq(1,0,length=50)),ylim=c(30,110),
   xlab=expression(beta[1]),ylab="b",main="(f)")
box()
dev.off()
# CI for b
exp(log(61.3)+1.96*17.2/61.3)
exp(log(61.3)-1.96*17.2/61.3)
