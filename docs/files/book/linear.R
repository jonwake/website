#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 5
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Prostate cancer example
#
library(lasso2)
data(Prostate)
attach(Prostate)
y <- Prostate$lpsa
lcavol <- Prostate$lcavol
lweight <- Prostate$lweight
age <- Prostate$age
lbph <- Prostate$lbph
svi <- Prostate$svi
lcp <- Prostate$lcp
gleason <- Prostate$gleason
pgg45 <- Prostate$pgg45
#
# Fig 5.1
#
pdf("prostatefig1.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(lcavol,y,xlab="log(cancer volume)",ylab="log(PSA)")
lines(lowess(lcavol,y))
dev.off()
#
# Fig 5.2
#
x <- runif(50,.1,6)
beta0 <- 78
beta1 <- 6.4
beta2 <- 3.25
beta3 <- -.89
y <- beta0 + beta1*x - exp(beta2+beta3*x) + rnorm(length(x),0,1.2)
pdf("jenss.pdf",height=4.5,width=4.5)
plot(x,y,xlab="Age (years)",ylab="Length (cms)")
xvals <- seq(min(x),max(x),.1)
lines(xvals,beta0 + beta1*xvals - exp(beta2+beta3*xvals),lty=1)
abline(lm(y[x>4]~x[x>4])$coeff,lty=3)
abline(lm(y[x<1]~x[x<1])$coeff,lty=2)
legend("bottomright",bty="n",legend=c("Full curve","Young Children",
       "Older Children"),lty=c(1,2,3))
dev.off()
#
# Fig 5.3(a)
#
pdf("interaction1.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,1))
points(c(0,1),c(0.5,0.5),pch=4)
points(c(0,1),c(0.5,0.5),pch=21)
lines(c(0,1),c(0.5,0.5),lwd=2)
dev.off()
#
# Fig 5.3(b)
#
pdf("interaction2.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,01))
points(c(0,1),c(0.2,0.8),pch=4)
points(c(0,1),c(0.2,0.8),pch=21)
lines(c(0,1),c(0.2,0.8),lwd=2)
dev.off()
#
# Fig 5.3(c)
#
pdf("interaction3.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,1))
points(c(0,1),c(0.1,0.1),pch=4)
points(c(0,1),c(0.9,0.9),pch=21)
lines(c(0,1),c(0.1,0.1),lwd=2)
lines(c(0,1),c(0.9,0.9))
dev.off()
#
# Fig 5.3(d)
#
pdf("interaction4.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,1))
points(c(0,1),c(0.1,0.6),pch=4)
points(c(0,1),c(0.5,1),pch=21)
lines(c(0,1),c(0.1,0.6),lwd=2)
lines(c(0,1),c(0.5,1),lwd=2)
dev.off()
#
# Fig 5.3(e)
#
pdf("interaction5.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,1))
points(c(0,1),c(0.1,0.2),pch=4)
points(c(0,1),c(0.3,1),pch=21)
lines(c(0,1),c(0.1,0.2),lwd=2)
lines(c(0,1),c(0.3,1),lwd=2)
lines(c(0,1),c(0.3,0.4),lty=2,lwd=2)
points(c(0,1),c(0.3,.4),pch=21)
dev.off()
#
# Fig 5.3(f)
#
pdf("interaction6.pdf",height=3,width=3)
x <- seq(0,1,.01)
y <- seq(0,1,.01)
plot(x,y,type="n",xlab=expression(x[1]),
   ylab=quote("E[Y|" *x[1]*","*x[2]*"]" ),axes=F)
axis(1,at=c(0,1))
points(c(0,1),c(0.1,0.8),pch=4)
points(c(0,1),c(0.3,0.2),pch=21)
lines(c(0,1),c(0.1,0.8),lwd=2)
lines(c(0,1),c(0.3,0.2),lwd=2)
lines(c(0,1),c(0.3,1.0),lty=2,lwd=2)
points(c(0,1),c(0.3,1.0),pch=21)
dev.off()
#
# Sandwich estimation simulation which produces Table 5.2
#
library(sandwich)
nsims <- 1000
beta0 <- 1
beta1 <- 1
sigma2 <- .001
counts0 <- array(0,dim=c(3,2,7))
counts1 <- array(0,dim=c(3,2,7))
counts2 <- array(0,dim=c(3,2,7))
counts3 <- array(0,dim=c(3,2,7))
nobs <- c(5,10,25,50,100,250,500)
#
#       edist = 1/2/3 for const, mu, mu^2 variance.
#       xdist =1/2 for normal, gamma errors.
#
for (edist in 1:3){
  for (xdist in 1:2){
    for (j in 1:length(nobs)){
       for (i in 1:nsims){
        if (xdist==1) x <- rnorm( nobs[j], 4, 1 )
        if (xdist==2) x <- rgamma( nobs[j], 1, 1 )
        mu <- beta0 + beta1*x
        if (edist == 1 ) error <- rnorm( nobs[j], 0, sqrt(sigma2) )
        else if (edist==2) error <- rnorm( nobs[j], 0, sqrt(mu*sigma2))
        else if (edist==3) error <- rnorm( nobs[j], 0, mu*sqrt(sigma2))
        y <- mu + error
        mod <- lm( y ~ x )
#
#       Form variance of hat{beta}_1 and confidence interval coverage.
#
        b1hat <- mod$coeff[2]
        seb1 <- sqrt(vcov(mod)[2,2])
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1 ) && 
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1) ) 
               counts0[edist,xdist,j] <- counts0[edist,xdist,j]+1 
#
#       Now for the pastrami (the filling of the sandwich...)
#
	xmat <- matrix(cbind(rep(1,nobs[j]),x),ncol=2,nrow=nobs[j])
        xtx <- t(xmat) %*% xmat
        xtxinv <- solve(xtx)
        Ainv <- xtxinv
        Bmat <- matrix(0,nrow=2,ncol=2)
        eps2 <- (y-mod$fit)^2
        Bmat[1,1] <- sum(eps2)
        Bmat[1,2] <- sum(eps2*x)
        Bmat[2,1] <- Bmat[1,2]
        Bmat[2,2] <- sum(eps2*x*x) 
        sandmat1 <- Ainv %*% Bmat %*% Ainv
        seb1sand1 <- sqrt( sandmat1[2,2] )
        sandmat2 <- sandmat1*nobs[j]/(nobs[j]-2)
        seb1sand2 <- sqrt( sandmat2[2,2] )
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1sand1 ) && 
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1sand1) )
               counts1[edist,xdist,j] <- counts1[edist,xdist,j]+1 
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1sand2 ) && 
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1sand2) ) 
               counts2[edist,xdist,j] <- counts2[edist,xdist,j]+1 
        Hmat <- xmat %*% xtxinv %*% t(xmat)
        diagh <- diag(Hmat)
        eps2h <- eps2/(1-diagh)
        Bmat2 <- matrix(0,nrow=2,ncol=2)
        Bmat2[1,1] <- sum(eps2h)
        Bmat2[1,2] <- sum(eps2h*x)
        Bmat2[2,1] <- Bmat2[1,2]
        Bmat2[2,2] <- sum(eps2h*x*x) 
        sandmat3 <- Ainv %*% Bmat2 %*% Ainv
        seb1sand3 <- sqrt( sandmat3[2,2] )
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1sand3 ) && 
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1sand3) ) 
               counts3[edist,xdist,j] <- counts3[edist,xdist,j]+1
	}
     }
   } 
}
#
# Results in the order of Table 5.2
for (edist in 1:3){
  for (xdist in 1:2){
    for (j in 1:length(nobs)){
        cat("OLS SAND1 SAND2 SAND3: ",100*counts0[edist,xdist,j]/nsims,
           100*counts1[edist,xdist,j]/nsims,100*counts2[edist,xdist,j]/nsims,
           100*counts3[edist,xdist,j]/nsims,"\n")
    }
  }
}
#
# Prostate cancer at the end of Section 5.6, Table 5.3 results
#
firstfit <- glm(y~lcavol,family=gaussian)
sqrt(sandwich(firstfit)[1,1])
sqrt(sandwich(firstfit)[2,2])
#
# Fig 5.4
#
pdf("prostatefig2.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1))
plot(lcavol,y,xlab="log(cancer volume)",ylab="log(PSA)")
abline(coef=firstfit$coef)
x0 <- seq(min(lcavol),max(lcavol),.01)
x0vec <- cbind(1,x0)
pred0 <- x0vec %*% vcov(firstfit) %*% t(x0vec)
sigmahat <- sum(firstfit$res^2)/firstfit$df.res
lines(x0,firstfit$coef[1] + firstfit$coef[2]*x0 + 1.96*sqrt(diag(pred0)),lty=2)
lines(x0,firstfit$coef[1] + firstfit$coef[2]*x0 - 1.96*sqrt(diag(pred0)),lty=2)
lines(x0,firstfit$coef[1] + firstfit$coef[2]*x0 + 
      1.96*(sigmahat+sqrt(diag(pred0))),lty=2)
lines(x0,firstfit$coef[1] + firstfit$coef[2]*x0 - 
      1.96*(sigmahat+sqrt(diag(pred0))),lty=2)
dev.off()
#
# Fig 5.5(a)
#
pdf("prostatefig3a.pdf",width=3.5,height=3.5)
par(mar=c(4,4,2,1)+.1)
plot(lcavol,exp(y),,xlab="log(cancer volume)",ylab="PSA")
xval <- seq(min((lcavol)),max((lcavol)),.01)
lines(xval,exp(firstfit$coef[1]+firstfit$coef[2]*(xval)))
dev.off()
#
# Fig 5.5(b)
#
pdf("prostatefig3b.pdf",width=3.5,height=3.5)
par(mar=c(4,4,2,1)+.1)
plot(exp(lcavol),exp(y),xlab="cancer volume",ylab="PSA")
exval <- seq(min(exp(lcavol)),max(exp(lcavol)),.01)
lines(exval,exp(firstfit$coef[1]+firstfit$coef[2]*log(exval)))
dev.off()
#
# Dyestuff data example in Section 5.8
#
dyedata <- read.table("dye.txt",header=T)
dyemod <- lm(dyedata$yield~as.factor(dyedata$batch))
anova(dyemod)
#
# Clotting data example in Section 5.8
#
y <- c(8.4,9.4,9.8,12.2,12.8,15.2,12.9,14.4,9.6,9.1,11.2,9.8,9.8,8.8,9.9,
  12.0,8.4,8.2,8.5,8.5,8.6,9.9,9.8,10.9,8.9,9.0,9.2,10.4,7.9,8.1,8.2,10.0)
tmt <- as.factor(rep(seq(1:4),8))
subject <- as.factor(rep(1:8,each=4))
twoway <- lm(y~tmt+subject)
anova(twoway)
#
# Fig 5.6
#
pdf("twowaytreat.pdf",width=4.5,height=3.5)
xvals <- seq(-1.3,3.6,.01)
mean2 <- twoway$coeff[2]
mean3 <- twoway$coeff[3]
mean4 <- twoway$coeff[4]
se2 <-sqrt(vcov(twoway)[2,2])
se3 <-sqrt(vcov(twoway)[3,3])
se4 <-sqrt(vcov(twoway)[4,4])
dens2 <- dt((xvals-mean2)/se2,df=21)/se2
dens3 <- dt((xvals-mean3)/se3,df=21)/se3
dens4 <- dt((xvals-mean4)/se4,df=21)/se4
plot(xvals,dens2,type="l",xlab="Treatment Difference",ylab="Posterior Density")
abline(v=0,lty=4)
lines(xvals,dens3,lty=2)
lines(xvals,dens4,lty=3)
legend("topright",legend=c("Tmt 2", "Tmt 3","Tmt 4"),bty="n",lty=1:3)
dev.off()
#
# Post probs that average responses are greater than tmt 1 for tmts 2,3,4
#
pt((mean2-0)/se2,df=21,lower=F)
pt((mean3-0)/se3,df=21,lower=F)
pt((mean4-0)/se4,df=21,lower=F)
#
# Bias-variance tradeoff example in Section 5.9
#
secondfit <- glm(y~lcp,family=gaussian)
thirdfit <- glm(y~lcavol+lcp,family=gaussian)
fourthfit <- glm(lcp~lcavol,family=gaussian)
checkcoef <- thirdfit$coeff[2] + fourthfit$coef[2]*thirdfit$coef[3]
cat("Predicted and observed = ",checkcoef,firstfit$coef[2],"\n")
#
# Fig 5.7(a)
#
pdf("prostatefig4a.pdf",height=3.5,width=3.5)
par(mfrow=c(1,1))
plot(lcavol,lcp,xlab="log(cancer volume)",ylab="log(cap pen)")
abline(fourthfit$coef)
dev.off()
#
# Fig 5.7(b)
#
pdf("prostatefig4b.pdf",height=3.5,width=3.5)
par(mfrow=c(1,1))
plot(lcp,y,xlab="log(cap pen)",ylab="log(PSA)")
abline(secondfit$coef)
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#	Simulation study that leads to Table 5.12 with errors generated from 
#       differnet errors distributions and then fitted with OLS and examine 
#       the coverage probs.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
dtgen <- function( x, df, mean, scale ){
        y <- (x-mean)/scale
        dtgen <- dt(y,df)
}
getdata <- function( nobs, sigma2, beta0, beta1, edist, xdist ){
        if (xdist == 1 ) x <- runif( nobs, 0, 1)
        else if (xdist == 2 ) x <- rgamma(nobs,1,1)
        if (edist == 1 ) error <- rnorm( nobs, 0, sqrt(sigma2))
        else if (edist == 2 ) {
                error <- rexp(nobs,1)
                z <- 2*rbinom(nobs,1,0.5)-1
                err <- z*error
        }       
        else if (edist == 3 ) error <- rt(nobs,3)
        else if (edist == 4 ) error <- rlnorm(nobs)
        y <- beta0 + beta1*x + error
        list( y=y, x=x )
}
nsims <- 10000
beta0 <- 1
beta1 <- 0
sigma2 <- 1
count <- array(0,dim=c(4,2,3))
seb1vec <- rep(0,nsims)
nobs <- c(5,20)
for (edist in 1:4){
  for (xdist in 1:2){
    for (j in 1:length(nobs)){
      for (i in 1:nsims){
        simdat <- getdata( nobs[j], sigma2, beta0, beta1, edist, xdist )
        x <- simdat$x
        y <- simdat$y
        mod <- lm( y ~ x )
#       Form se of hat{beta}_1 and confidence interval coverage.
        b1hat <- mod$coeff[2]
	seb1 <- sqrt(vcov(mod)[2,2])
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1 )&&
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1) ) 
             count[edist,xdist,j] <- count[edist,xdist,j]+1 
      }
# Table 5.12 results prooduced from this
      cat(edist,xdist,nobs[j],round(100*count[edist,xdist,j]/nsims,
         digits=0),"\n")
    }
  }
}
#
# Now examine the effect of correlated observations in Section 5.10
#
getdata2 <- function( nobs, sigma2, beta0, beta1, rho, xdist ){
        if (xdist == 1 ) x <- runif( nobs, 0, 1)
        else if (xdist == 2 ) x <- rgamma(nobs,1,1)
	x <- sort(x)
        cove <- matrix( 1, nrow=nobs, ncol=nobs )
        for (i in 1:nobs){
                for (j in 1:nobs){ cove[i,j]<-sigma2*rho^(abs(i-j))}
        }
        error <- as.vector(mvrnorm( 1, mu=rep(0,nobs), Sigma=cove) )
        y <- beta0 + beta1*x + error
        list( y=y, x=x )
}
library(MASS)
nsims <- 10000
nobs <- c(5,20,50)
beta0 <- 1
beta1 <- 1
sigma2 <- 1
rho <- c(0.1,0.5,0.95)
count <- array(0,dim=c(2,3,3))
for (xdist in 1:2){
  for (corr in 1:3){
    for (j in 1:length(nobs)){
       for (i in 1:nsims){
        simdat <- getdata2( nobs[j], sigma2, beta0, beta1, rho[corr], xdist )
        x <- simdat$x
        y <- simdat$y
        mod <- lm( y ~ x )
#       Find confidence interval coverage.
        b1hat <- mod$coeff[2]
	seb1 <- sqrt(vcov(mod)[2,2])
        if ( (beta1 > b1hat - qt(0.975,nobs[j]-2)*seb1 ) && 
             (beta1 < b1hat + qt(0.975,nobs[j]-2)*seb1) ) 
             count[xdist,corr,j] <- count[xdist,corr,j]+1 
      }
# Table 5.13 results prooduced from this command
     cat(xdist,rho[corr],nobs[j],round(100*count[xdist,corr,j]/nsims,
     digits=0),"\n")
   }
  }
}
#
# Fig 5.8 (QQ plots)
#
pdf("qq.pdf",height=4.5,width=4.5)
par(mfrow=c(4,4))
par(mar=c(1,1,1,1)+.1)
b0 <- 0
b1 <- 1
#	Normal
n <- 10; x<-rnorm(n,0,1); y <- b0 + b1*x + rnorm(n,0,1)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="n=10"); box()
n <- 25
x<-rnorm(n,0,1)
y <- b0 + b1*x + rnorm(n,0,1)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="n=25");box()
n <- 50; x<-rnorm(n,0,1); y <- b0 + b1*x + rnorm(n,0,1)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="n=50");box()
n <- 200
x<-rnorm(n,0,1)
y <- b0 + b1*x + rnorm(n,0,1)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="n=200");box()
#	Laplacian
n <- 10
x<-rnorm(n,0,1); err <- rexp(n,1); z<-2*rbinom(n,1,0.5)-1; err <- z*err
y <- b0 + b1*x + err
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 25
x<-rnorm(n,0,1); err <- rexp(n,1); z<-2*rbinom(n,1,0.5)-1; err <- z*err
y <- b0 + b1*x + err
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 50
x<-rnorm(n,0,1)
err <- rexp(n,1)
z<-2*rbinom(n,1,0.5)-1
err <- z*err
y <- b0 + b1*x + err
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 200
x<-rnorm(n,0,1)
err <- rexp(n,1)
z<-2*rbinom(n,1,0.5)-1
err <- z*err
y <- b0 + b1*x + err
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
#	t_3
n <- 10; x<-rnorm(n,0,1); y <- b0 + b1*x + rt(n,3)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 25; x<-rnorm(n,0,1); y <- b0 + b1*x + rt(n,3)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 50
x<-rnorm(n,0,1)
y <- b0 + b1*x + rt(n,3)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 200
x<-rnorm(n,0,1)
y <- b0 + b1*x + rt(n,3)
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
#	Lognormal
n <- 10; x<-rnorm(n,0,1); y <- b0 + b1*x + exp(rnorm(n,0,1))
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 25; x<-rnorm(n,0,1); y <- b0 + b1*x + exp(rnorm(n,0,1))
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 50; x<-rnorm(n,0,1); y <- b0 + b1*x + exp(rnorm(n,0,1))
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
n <- 200; x<-rnorm(n,0,1); y <- b0 + b1*x + exp(rnorm(n,0,1))
mod <- lm(y~x); eps <- (y - mod$fit)/sqrt(sum(mod$residuals^2)/(n-2))
qqnorm(eps,xlab="",ylab="",axes=F,main="");box()
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Prostate cancer example in Section 5.12
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Table 5.14 results 
#
y <- Prostate$lpsa
lcavol <- Prostate$lcavol
lweight <- Prostate$lweight
age <- Prostate$age
lbph <- Prostate$lbph
svi <- Prostate$svi
lcp <- Prostate$lcp
gleason <- Prostate$gleason
pgg45 <- Prostate$pgg45
fullmod <- lm(y~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45)
summary(fullmod)
sqrt(diag(sandwich(fullmod)))
#
resid <- (y-fullmod$fit)/(sum((y-fullmod$fit)^2)/fullmod$df.resid)
#
# Fig 5.9(a)
#
pdf("prostatediaga.pdf",height=2.5,width=2.5)
par(mfrow=c(1,1),mar=c(5,4,1,1)+.1)
plot(fullmod$fit,abs(resid),xlab="Fitted Values",ylab="|Residuals|")
lines(lowess(fullmod$fit,abs(resid)))
dev.off()
#
# Fig 5.9(b)
#
pdf("prostatediagb.pdf",height=2.5,width=2.5)
par(mfrow=c(1,1),mar=c(5,4,1,1)+.1)
qqnorm(resid,ylab="Observed",xlab="Expected",main="")
dev.off()
#
# Fig 5.9(c)
#
pdf("prostatediagc.pdf",height=2.5,width=2.5)
par(mfrow=c(1,1),mar=c(5,4,1,1)+.1)
plot(lcavol,resid,xlab="log(can vol)",ylab="Residuals")
lines(lowess(lcavol,resid)); abline(0,0,lty=2)
dev.off()
#
# Fig 5.9(d)
#
pdf("prostatediagd.pdf",height=2.5,width=2.5)
par(mfrow=c(1,1),mar=c(5,4,1,1)+.1)
plot(lweight,resid,xlab="log(weight)",ylab="Residuals")
lines(lowess(lweight,resid,f=.95)); abline(0,0,lty=2)
dev.off()
#
# Now let's standardize the x variables for the Bayesian analysis
#
standx <- function(x) z<-(x-range(x)[1])/(range(x)[2]-range(x)[1])
par(mfrow=c(4,2))
x1 <- standx(lcavol); hist(x1)
x2 <- standx(lweight); hist(x2)
x3 <- standx(age); hist(x3)
x4 <- standx(lbph); hist(x4)
x5 <- standx(svi); hist(x5)
x6 <- standx(lcp); hist(x6)
x7 <- standx(gleason); hist(x7)
x8 <- standx(pgg45); hist(x8)
fullmod2 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)
xvals <- seq(-2.3,4.5,.01)
mean1 <- fullmod2$coeff[2]
mean2 <- fullmod2$coeff[3]
mean3 <- fullmod2$coeff[4]
mean4 <- fullmod2$coeff[5]
mean5 <- fullmod2$coeff[6]
mean6 <- fullmod2$coeff[7]
mean7 <- fullmod2$coeff[8]
mean8 <- fullmod2$coeff[9]
se1 <-sqrt(vcov(fullmod2)[2,2])
se2 <-sqrt(vcov(fullmod2)[3,3])
se3 <-sqrt(vcov(fullmod2)[4,4])
se4 <-sqrt(vcov(fullmod2)[5,5])
se5 <-sqrt(vcov(fullmod2)[6,6])
se6 <-sqrt(vcov(fullmod2)[7,7])
se7 <-sqrt(vcov(fullmod2)[8,8])
se8 <-sqrt(vcov(fullmod2)[9,9])
dens1 <- dt((xvals-mean1)/se1,df=fullmod2$df.residual)/se1
dens2 <- dt((xvals-mean2)/se1,df=fullmod2$df.residual)/se2
dens3 <- dt((xvals-mean3)/se1,df=fullmod2$df.residual)/se3
dens4 <- dt((xvals-mean4)/se1,df=fullmod2$df.residual)/se4
dens5 <- dt((xvals-mean5)/se1,df=fullmod2$df.residual)/se5
dens6 <- dt((xvals-mean6)/se1,df=fullmod2$df.residual)/se6
dens7 <- dt((xvals-mean7)/se1,df=fullmod2$df.residual)/se7
dens8 <- dt((xvals-mean8)/se1,df=fullmod2$df.residual)/se8
#
# Fig 5.10 -- posterior distributions under a flat prior
#
pdf("prostbayes.pdf",width=7.5,height=7.5)
par(mfrow=c(1,1),mar=c(5, 4, 4, 2) +0.1)
plot(xvals,dens1,type="l",xlab="Coefficient",ylab="Posterior Density",
     ylim=c(0,1.9))
abline(v=0,lty=4)
lines(xvals,dens2,lty=2)
lines(xvals,dens3,lty=3)
lines(xvals,dens4,lty=4)
lines(xvals,dens5,lty=5)
lines(xvals,dens6,lty=6)
lines(xvals,dens7,lty=7)
lines(xvals,dens8,lty=8)
legend("topright",legend=c("log can vol","log weight","age","log bph","svi",
       "log cap pen","gleason","pgg45"),bty="n",lty=1:8)
dev.off()
#
# These are the Bayesian tail areas
#
pt((mean1-0)/se1,df=fullmod2$df.res,lower=F)
pt((mean2-0)/se2,df=fullmod2$df.res,lower=F)
pt((mean3-0)/se3,df=fullmod2$df.res,lower=F)
pt((mean4-0)/se4,df=fullmod2$df.res,lower=F)
pt((mean5-0)/se5,df=fullmod2$df.res,lower=F)
pt((mean6-0)/se6,df=fullmod2$df.res,lower=F)
pt((mean7-0)/se7,df=fullmod2$df.res,lower=F)
pt((mean8-0)/se8,df=fullmod2$df.res,lower=F)
#
xprec <- (1.96/log(10))^2
#
# Informative prior analysis using INLA
#
source("http://www.math.ntnu.no/inla/givemeINLA.R")
inla.upgrade()
library(INLA)
prostdat <- list(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8)
# The default settings in inla are with a very large variance
inlamod1 <- inla(y~x1+x2+x3+x4+x5+x6+x7+x8,data=prostdat)
inlamod2 <- inla(y~x1+x2+x3+x4+x5+x6+x7+x8,control.fixed = list(prec=xprec), 
            data=prostdat)
#, prec.intercept=1))
#
summary(inlamod2)
names(inlamod2)
yval <- seq(1,8)
#
# Fig 5.11
#
pdf("prostatepostint.pdf",width=6.5,height=6.5)
par(mar=c(5,7,1,1)+.1)
plot(yval,yval,xlim=c(-2.2,4.2),type="n",xlab="95% Credible Interval",
     ylab="",axes="F")
box()
axis(1)
axis(2,at=seq(1,8),labels=c("log can vol","log weight","age","log bph","svi",
                   "log cap pen","gleason","pgg45"),las=1)
p1lo <- p1hi <- p2lo <- p2hi <- rep(0,8)
for (i in 1:8){
    p1lo[i] <- inlamod1$summary.fixed[i+1,3]
    p1hi[i] <- inlamod1$summary.fixed[i+1,5]
    p2lo[i] <- inlamod2$summary.fixed[i+1,3]
    p2hi[i] <- inlamod2$summary.fixed[i+1,5]
    lines(y=c(i+.1,i+.1),x=c(p1lo[i],p1hi[i]))
    lines(y=c(i-.1,i-.1),x=c(p2lo[i],p2hi[i]),lty=2)
}
abline(v=0,lty=4)
legend("topright",legend=c("Flat Prior","Informative Prior"),lty=1:2,bty="n")
dev.off()