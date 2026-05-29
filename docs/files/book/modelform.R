#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 4
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Function to implement Selke et al (2001) American Statistician relationship
#
lowerp <- function(pi1,p){
  RPO <- pi1/(1-pi1)
  if ( p < 1/exp(1) ){
    lowerp <- 1/(1-RPO/(exp(1)*p*log(p)))
  }
  else lowerp <- 0
  return(lowerp)
}
#
# Lower bounds on posterior of null given a P-value and as a function of priors
#
pi0 <- c(0.25,0.5,0.75)
logp <- seq(log(0.001),log(.3),.01)
p <- exp(logp)
lobound <- matrix(0,nrow=length(pi0),ncol=length(p))
for (k in 1:length(pi0)){
  for (j in 1:length(p)){
    lobound[k,j] <- lowerp(1-pi0[k],p[j])
  }
}
#
# Fig 4.1
#
pdf("sellk.pdf",width=4.5,height=4.5)
plot(p,lobound[1,],type="n",ylim=c(0,1),ylab="Posterior probability of the null",xlab="p-value")
lines(p,lobound[1,],lty=1)
lines(p,lobound[2,],lty=2)
lines(p,lobound[3,],lty=3)
legend("topleft",xjust=.5,legend=c(expression(pi[0]==0.75),expression(pi[0]==0.50),
    expression(pi[0]==0.25)),lty=c(3,2,1),bty="n",lwd=2)
dev.off()
#
# Examine Lindley's paradox for a normal prior N(0,tau2) on theta.
# Keep the p-value constant but change n.
# p-value = 2*(1-\Phi(overline{y}_n\sqrt{n}/\sigma))
#
#
# Fig 4.2
#
pdf("newlindley1.pdf",height=7.5,width=7.5)
sigma2 <- 1
rectau2 <- 1/(.2^2)
tau2 <- 1/rectau2
n <- 100
pvalue <- 0.05
lambda <- qnorm(1-pvalue/2)
y <- seq(-.5,.5,.01)
fy0 <- dnorm(y,m=0,s=sqrt(sigma2/n))
fy1 <- dnorm(y,m=0,s=sqrt(sigma2/n+tau2)); par(mfrow=c(1,1))
plot(y,fy0,type="n",xlab=expression(bar(y)[n]),ylab="Density")
lines(y,fy0,lty=1)
lines(y,fy1,lty=2)
legend("topleft",legend=c("Under the null","Under the alternative"),bty="n",lty=c(1,2),lwd=2)
ybarstar <- lambda*sqrt(sigma2/n)
abline(v=ybarstar,lty=4,lwd=1.5)
BFlam <- dnorm(ybarstar,m=0,s=sqrt(sigma2/n))/dnorm(ybarstar,m=0,s=sqrt(sigma2/n+tau2))
dev.off()
#
# Check two different Bayes factor calculations for example in Section 4.4
#
cat("Two forms: ",BFlam,sqrt((sigma2+n*tau2)/sigma2)*exp(-.5*lambda^2*tau2/(tau2+sigma2/n)),"\n")
#
#
#
ImpSamp <- function(ybar,n,sigma2,tau2){
  nsim <- 10000
  thetasamp <- rnorm(nsim,m=0,s=sqrt(tau2))
  zl <- (-ybar-thetasamp)*sqrt(n)/sqrt(sigma2)
#  zu <- (ybar-thetasamp)*sqrt(n)/sqrt(sigma2)
#  probu <- mean(pnorm(zu,lower.tail=F)) # these should be equal since the prior is 
  probl <- mean(pnorm(zl,lower.tail=T)) # centered about 0 and the test is 2-sided.
  ImpSamp <- 2*probl #+probu
  list(ImpSamp=ImpSamp,probl=probl,probu=probu)
}
n <- exp(seq(log(1),log(100000),length=1000)) # effect of increasing power
ybar <- lambda*sqrt(sigma2/n)
pi0 <- .5; PO <- pi0/(1-pi0); avepower <- probl <- probu <- NULL
if (rectau2 > 0) {
  for (i in 1:length(n)){
    temp <- ImpSamp(ybar[i],n[i],sigma2,tau2)
    avepower[i] <- temp$ImpSamp
    probl[i] <- temp$probl
    probu[i] <- temp$probu
  }
}
#
# Fig 4.4(a)
#
pdf("lindley3.pdf",height=4.5,width=4.5)
plot(n,avepower,type="n",ylim=c(0,1),log="x",xlab="n",ylab="Average Power")
lines(smooth.spline(n,avepower,nknots=20))
dev.off()
freqBF <- pvalue/avepower
#
# Fig 4.4(b)
#
pdf("lindley4.pdf",height=4.5,width=4.5)
plot(n,freqBF,type="n",ylim=c(0,1),ylab="Tail-area Bayes factor",log="x")
lines(smooth.spline(n,freqBF))
abline(h=0.05,lty=2)
dev.off()
if (rectau2 > 0) {
  tau2 <- 1/rectau2 
#  BF01 <- sqrt((n*tau2+sigma2)/sigma2)*exp(-n*n*ybar*ybar*tau2/(2*sigma2*(n*tau2+sigma2)))
  num <- dnorm(ybar,m=0,s=sqrt(sigma2/n))
  denom <- dnorm(ybar,m=0,s=sqrt(sigma2/n+tau2))
  check <- num/denom
}
if (rectau2==0) { # flat prior on theta to give eqn (2) of Lindley
  check <- exp(-n*ybar^2/(2*sigma2))/sqrt(2*pi*sigma2/n)
  cbar <- pi0*exp(-.5*lambda^2)/(pi0*exp(-.5*lambda^2)+(1-pi0)*sqrt(2*pi*sigma2/n))
}
postodds <- check*PO
#
# Fig 4.3
#
pdf("lindley1.pdf",height=7,width=7.5)
postprob <- postodds/(1+postodds)
plot(n,postprob,type="n",ylim=c(0,1),xlab="n",ylab="Posterior probability of the null",log="x");
lines(n,postprob)
abline(h=0.5,lty=2)
dev.off()
#
# Microarray data example of Section 4.6
#
# Read in data
#
ngenes <- 1000
x1 <- read.table("GeneSamp1.txt",header=T,nrows=ngenes)
x2 <- read.table("GeneSamp2.txt",header=T,nrows=ngenes)
z <- apply(log2(x1[,1:60]),1,mean) - apply(log2(x2[,1:45]),1,mean)
sigma2 <- apply(log2(x1[,1:60]),1,var)/60 + apply(log2(x2[,1:45]),1,var)/45
#
# Compute p-values
#
pvalue <- NULL
for (i in 1:length(z)){
    if (z[i]>=0) pvalue[i] <- 2*(1-pnorm(z[i],0,sd=sqrt(sigma2[i])))
    else if (z[i]<0) pvalue[i] <- 2*pnorm(z[i],0,sd=sqrt(sigma2[i]))
}
#
# Fig 4.5(a)
#
pdf("microzscores.pdf",height=4,width=4)
par(mfrow=c(1,1))
hist(z/sqrt(sigma2),xlab="Z scores",main="",prob=T,ylim=c(0,.4))
x <- seq(min(z/sqrt(sigma2)),max(z/sqrt(sigma2)),.0001)
curve(dnorm(x,mean=0,sd=1),add=T,lty=2)
dev.off()
#
# Fig 4.5(b)
#
pdf("micropvals.pdf",height=4,width=4)
par(mfrow=c(1,1))
hist(pvalue,xlab="p-values",main="")
dev.off()
#
# Multiple testing p-value plot.
#
m <- 100
alpha <- 0.05
p <- runif(m)
expp <- seq(1,m)/(m+1)
p[1:5] <- expp[1:5]/seq(50,5,l=5)
p <- sort(p)
#
# Fig 4.6 -- will not be the same as in book because data is random
#
pdf("multtestfig.pdf",height=4.5,width=4.5)
ymax <- max(-log10(p),-log10(0.05/m),-log10(1/m))
plot(-log10(expp),-log10(p),xlim=c(0,ymax),ylim=c(0,ymax),xlab="-log10(expected)",ylab="-log10(observed)")
abline(h=-log10(0.05/m),lty=2)
abline(h=-log10(1/m),lty=3)
abline(a=0,b=1)
FDRcon <- 0.05
abline(a=-log10(FDRcon),b=1,lty=4)
legend("bottomright",legend=c("Bonferroni 5%","EFD=1","FDR 5%"),bty="n",lwd=1,lty=2:4)
dev.off()
#
# q-values for microarray data in Section 4.6
#
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue) 
run <- qvalue(p=pvalue)
summary(run)
#
# Bonferroni and FDR significant calls
#
bonf <- 0.05/ngenes
cat("Number significant under Bonferroni at 0.05= ",sum(pvalue<bonf),"\n")
# now check
sum(p.adjust(pvalue,method="bonferroni")<0.05)
sum(p.adjust(pvalue,method="BH")<0.05)
#
logcond <- function(z,w,delta,sigma2,sigma21){
	logcond <- 0
	for (i in 1:length(z)){
	    logcond <- logcond + dnorm(z[i],mean=delta*w[i],sd=sqrt(sigma2[i]+sigma21*w[i]),log=T)
	}
logcond
}
#
# Bayes analysis of microarray data treating the data on the m transcripts as independent
#
tau2 <- (log2(1.1)/qnorm(0.975))^2
rat <- tau2/(tau2+sigma2)
zscore <- z/sqrt(sigma2)
BFs <- exp(-0.5*zscore^2*rat)/sqrt(1-rat)
#
# Fig 4.7
#
pdf("microarrayindBF.pdf",height=6.5,width=6.5)
par(mfrow=c(1,1))
plot(sort(-log10(BFs)),xlab="Transcript",ylab="-log10 Bayes factor")
abline(0,0,lty=2)
dev.off()
pi0 <- seq(0.005,.995,.005)
postodds <- matrix(0,nrow=ngenes,ncol=length(pi0))
thresh <- 0.5
for (i in 1:ngenes){
    postodds[i,] <- BFs[i]*pi0/(1-pi0)
}
sig <- NULL
for (j in 1:length(pi0)){
    sig[j] <- sum(postodds[,j] < thresh/(1-thresh))
}
#
# Fig 4.8
#
pdf("microarrayindBayes.pdf",height=6.5,width=6.5)
par(mfrow=c(1,1))
plot(pi0,sig,xlab=expression(pi[0]),ylab="Number of significant transcripts",ylim=c(0,ngenes))
dev.off()
#
# Hierarchical model fitted to all of the microarray data using MCMC 
# with auxiliary variables w.
#
nits <- 10000
delta <- sigma21 <- pi0 <- w <- NULL
count1 <- 0
delta[1] <- 0; sigma21[1] <- var(z); pi0[1] <- .8
for (j in 2:nits){
  for (i in 1:ngenes){
    term1 <- exp(dnorm(z[i],mean=0,sd=sqrt(sigma2[i]),log=T))
    term2 <- exp(dnorm(z[i],mean=delta[j-1],sd=sqrt(sigma2[i]+sigma21[j-1]),
             log=T))
    p <- term2*(1-pi0[j-1])/( term1*pi0[j-1] + term2*(1-pi0[j-1]) )
    w[i] <- rbinom(n=1,size=1,prob=p)
if (is.na(w[i]==T))cat("p = ",p,z[i],delta[j-1],sigma21[j-1],sigma2[i],
      pi0[j-1],"\n")
  }
#cat(w,"\n")
# Conditional for delta
  y1 <- z[w==1]
  sig21 <- sigma2[w==1]
  A <- sum(1/(sig21+sigma21[j-1]))
  B <- sum(y1/(sig21+sigma21[j-1]))
  delta[j] <- rnorm(1,mean=B/A,sd=sqrt(1/A))
#cat(A,B,delta[j],"\n")
# MH for sigma21 
  old <- log(sigma21[j-1])
  new <- old+rnorm(1,mean=0,sd=.1)
  oldcon <- logcond(z,w,delta[j],sigma2,exp(old))
  newcon <- logcond(z,w,delta[j],sigma2,exp(new))
#cat("oldcon newcon: ",oldcon,newcon,"\n")
#cat(old,new,oldcon,newcon,delta[j],sigma21[j-1],"\n")
  if ( log(runif(1)) < newcon-oldcon ) {count1<-count1+1;sigma21[j] <- exp(new)}
  else sigma21[j] <- exp(old)
  a <- ngenes-sum(w)+1
  b <- sum(w)+1
  pi0[j] <- rbeta(1,a,b)
#  cat("Iteration = ",j,a,b,delta[j],sigma21[j],pi0[j],"\n")
}
#
# Summarize and plots
#
burnin <- 1000; p <- NULL
for (i in 1:ngenes){
    term1 <- exp(dnorm(z[i],mean=0,sd=sqrt(sigma2[i]),log=T))
    term2 <- exp(dnorm(z[i],mean=delta[seq(burnin+1,nits)],
             sd=sqrt(sigma2[i]+sigma21[seq(burnin+1,nits)]),log=T))
    p[i] <- mean( term2*(1-pi0[seq(burnin+1,nits)])/( term1*pi0[seq(burnin+1,nits)] + 
            term2*(1-pi0[seq(burnin+1,nits)])) )
}
#
# Fig 4.10
#
cat("Posterior prob > 0.5 = ",sum(p>0.5),"\n")
pdf("microPH.pdf",height=6,width=4.5)
par(mfrow=c(1,1))
plot(sort(p),ylab="Posterior probability of H=1",xlab="Transcript",ylim=c(0,1))
abline(0.5,0,lty=2)
dev.off()
#
sum(sort(p)[1:311])/311
sum(sort(1-p)[1:689])/689
sum(sort(p)[1:397])/397
sum(sort(1-p)[1:603])/603
#
par(mfrow=c(2,2))
plot(delta);plot(sigma21);plot(pi0)
hist(z,prob=T)
x <- seq(min(z),max(z),.001)
curve(dnorm(x,mean=quantile(delta[1001:nits],p=.5),
      sd=sqrt(quantile(sigma21[1001:nits],p=.5))),add=T,lty=2)
cat("accept ratio for sigma21 = ",count1/nits,"\n")
cat("POSTERIOR SUMMARIES:\n")
cat(quantile(delta[seq(burnin+1,nits)],p=c(0.025,0.5,0.975)),"\n")
cat(quantile(sigma21[seq(burnin+1,nits)],p=c(0.025,0.5,0.975)),"\n")
cat(quantile(pi0[seq(burnin+1,nits)],p=c(0.025,0.5,0.975)),"\n")
#
# Fig 4.9(a)
#
pdf("micropost1.pdf",height=4,width=3)
par(mfrow=c(1,1))
hist(delta[seq(burnin+1,nits)],xlab=expression(delta),main="")
dev.off()
#
# Fig 4.9(b)
#
pdf("micropost2.pdf",height=4,width=3)
par(mfrow=c(1,1))
hist(sigma21[seq(burnin+1,nits)],xlab=expression(tau^2),main="")
dev.off()
#
# Fig 4.9(c)
#
pdf("micropost3.pdf",height=4,width=3)
par(mfrow=c(1,1))
hist(pi0[seq(burnin+1,nits)],xlab=expression(pi[0]),main="")
dev.off()
library(MASS)
#
# Fig 4.9(d)
#
pdf("micropost4.pdf",height=4,width=3)
par(mfrow=c(1,1))
image(kde2d(delta[seq(burnin+1,nits)],sigma21[seq(burnin+1,nits)],n=50),col=grey(seq(1,0,length=50)),
   xlab=expression(delta),ylab=expression(tau^2))
box()
dev.off()
#
# Fig 4.9(e)
#
pdf("micropost5.pdf",height=4,width=3)
par(mfrow=c(1,1))
image(kde2d(delta[seq(burnin+1,nits)],pi0[seq(burnin+1,nits)],n=50),
   col=grey(seq(1,0,length=50)),xlab=expression(delta),ylab=expression(pi[0]))
box()
dev.off()
#
# Fig 4.9(f)
#
pdf("micropost6.pdf",height=4,width=3)
par(mfrow=c(1,1))
image(kde2d(sigma21[seq(burnin+1,nits)],pi0[seq(burnin+1,nits)],n=50),
   col=grey(seq(1,0,length=50)),xlab=expression(tau^2),ylab=expression(pi[0]))
box()
dev.off()
#
# Example to demonstrate selection bias in Section 4.8.
#
getdata <- function( nobs, nsims, beta0, beta1, beta2, sigma ){
        error <- rnorm( nsims*nobs, mean=0, sd=sigma )
        covx <- matrix( c(1,0.7,0.7,1), nrow=2, ncol=2 )
        x <- mvrnorm( nobs, mu=c(0,0), Sigma=covx )
        x1 <- x[,1]
        x2 <- x[,2]
        y <- matrix( beta0 + beta1 * x1 + beta2 * x2 + error,
                     nrow=nobs, ncol=nsims )
        list( y=y, x=x, x1=x1, x2=x2, error=error )
}
library(MASS)
beta0 <- 1
beta1 <- 1
beta2 <- 1
nobs <- 10
nsims <- 50000
sigma <- 3.0
simdat <- getdata( nobs, nsims, beta0, beta1, beta2, sigma )
x1 <- simdat$x1
x2 <- simdat$x2
y <- simdat$y
b1hat <- b2hat <- b1hath <- NULL
for (j in 1:nsims) {
       mod <- lm( y[,j] ~ x1+x2 )
       b1hat[j] <- mod$coeff[2]
       b2hat[j] <- mod$coeff[3]
       seb2 <- sqrt(vcov(mod)[3,3])
       if (abs(b2hat[j]) > 1.96*seb2 ) {
          b1hath[j] <- mod$coeff[2] 
       }
       else {
          mod1 <- lm( y[,j] ~ x1 )
          b1hath[j] <- mod1$coeff[2]
       }
}  
#
# Fig 4.11(a)
#
pdf("selbias2Ra.pdf",height=6.5,width=6.5)  
par( mfrow = c(1,1) )
hist(b1hat,xlim=c(min(b1hat,b1hath)-.5,max(b1hat,b1hath)+.5),
      xlab=expression(paste(hat(beta)[1],", without adjustment, ")),
      main="",nclass=30)
dev.off()
#
# Fig 4.11(b)
#
pdf("selbias2Rb.pdf",height=6.5,width=6.5)  
par( mfrow = c(1,1) )
hist( b1hath,xlim=c(min(b1hat,b1hath)-.5,max(b1hat,b1hath)+.5),main="",
       xlab = expression(paste(hat(beta)[1],", with adjustment, ")),nclass=30)
dev.off()
mean(b1hat)
mean(b1hath)
sqrt(var(b1hat))
sqrt(var(b1hath))
#
# Prostate data example
#
library(lasso2)
data(Prostate)
attach(Prostate)
#
# Backwards elimination
#
fullmod <- lm(lpsa~.,data=Prostate)
summary(fullmod) # Remove gleason
back2 <- lm(lpsa~lcavol+lweight+age+lbph+svi+lcp+pgg45)
summary(back2) # remove lcp
back3 <- lm(lpsa~lcavol+lweight+age+lbph+svi+pgg45)
summary(back3) # remove pgg45
back4 <- lm(lpsa~lcavol+lweight+age+lbph+svi)
summary(back4) # remove age
back5 <- lm(lpsa~lcavol+lweight+lbph+svi)
summary(back5) # remove lbph
back6 <- lm(lpsa~lcavol+lweight+svi)
summary(back6) # all p-values now less than 0.05.
#
# forwards selection
#
fall1 <- lm(lpsa~lcavol)
fall2 <- lm(lpsa~lweight)
fall3 <- lm(lpsa~age)
fall4 <- lm(lpsa~lbph)
fall5 <- lm(lpsa~svi)
fall6 <- lm(lpsa~lcp)
fall7 <- lm(lpsa~pgg45)
fall8 <- lm(lpsa~gleason)
#
# add in lcavol
#
fall21 <- lm(lpsa~lcavol+lweight)
fall22 <- lm(lpsa~lcavol+age)
fall23 <- lm(lpsa~lcavol+lbph)
fall24 <- lm(lpsa~lcavol+svi)
fall25 <- lm(lpsa~lcavol+lcp)
fall26 <- lm(lpsa~lcavol+pgg45)
fall27 <- lm(lpsa~lcavol+gleason)
#
# add in lweight
#
fall31 <- lm(lpsa~lcavol+lweight+age)
fall32 <- lm(lpsa~lcavol+lweight+lbph)
fall33 <- lm(lpsa~lcavol+lweight+svi)
fall34 <- lm(lpsa~lcavol+lweight+lcp)
fall35 <- lm(lpsa~lcavol+lweight+pgg45)
fall36 <- lm(lpsa~lcavol+lweight+gleason)
#
# add in svi
#
fall41 <- lm(lpsa~lcavol+lweight+svi+age)
fall42 <- lm(lpsa~lcavol+lweight+svi+lbph)
fall43 <- lm(lpsa~lcavol+lweight+svi+lcp)
fall44 <- lm(lpsa~lcavol+lweight+svi+pgg45)
fall45 <- lm(lpsa~lcavol+lweight+svi+gleason)
# 
###################
# Cp - note p includes the intercept
###################
#
Cpplotprost <- function (cp) 
{
    p <- max(cp$size)
    i <- (cp$Cp < (p + 1.5))
    plot(cp$size[i], cp$Cp[i], xlab = "p", ylab = expression(C[p]), 
         type = "n",xlim=c(3.3,9.7),ylim=c(3.3,9.7))
    labels <- apply(cp$which, 1, function(x) 
              paste(as.character((1:(p - 1))[x]), collapse = ""))
    text(cp$size[i], cp$Cp[i], labels[i])
    abline(0, 1)
}
####
library(leaps)
x <- model.matrix(fullmod)[,-1]
y <- Prostate$lpsa
cp <- leaps(x,y)
library(faraway)
pdf("prostate-cp.pdf",height=6.5,width=6.5)
Cpplotprost(cp)
dev.off()
#########################
adjR2 <- leaps(x,y,method="adjr2")
maxadjr(adjR2,8)
# best model has everything but 7 in
reg <- regsubsets(lpsa~.,data=Prostate)
summary(reg)$bic
#
# BMA
#
library(BMA)
library(MASS)
library(lasso2)
data(Prostate)
fullmod <- lm(lpsa~.,data=Prostate)
x <- model.matrix(fullmod)[,-1]
y <- Prostate$lpsa
glibmod <- glib(x,y,error="gaussian",link="identity")
summary(glibmod)