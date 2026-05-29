#
#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 11
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# Polynomials for the LIDAR data
#
library(SemiPar)
data(lidar)
attach(lidar)
xvals <- seq(min(range),max(range),.1)
lidar2 <- lm(logratio~poly(range,deg=2,raw=T))
lidar3 <- lm(logratio~poly(range,deg=3,raw=T))
lidar4 <- lm(logratio~poly(range,deg=4,raw=T))
lidar8 <- lm(logratio~poly(range,deg=8,raw=T))
#
# Fig 11.1(a)
#
pdf("lidardeg2fig.pdf",height=4.5,width=4.5)
plot(logratio~range,xlab="Range (m)",ylab="Log Ratio",col="grey")
lines(xvals,lidar2$coeff[1]+lidar2$coeff[2]*xvals+lidar2$coeff[3]*xvals*xvals)
dev.off()
#
# Fig 11.1(b)
#
pdf("lidardeg3fig.pdf",height=4.5,width=4.5)
plot(logratio~range,xlab="Range (m)",ylab="Log Ratio",col="grey")
lines(xvals,lidar3$coeff[1]+lidar3$coeff[2]*xvals+lidar3$coeff[3]*xvals*xvals+
      lidar3$coeff[4]*xvals*xvals*xvals,lty=1)
dev.off()
#
# Fig 11.1(c)
#
pdf("lidardeg4fig.pdf",height=4.5,width=4.5)
plot(logratio~range,xlab="Range (m)",ylab="Log Ratio",col="grey")
lines(xvals,lidar4$coeff[1]+lidar4$coeff[2]*xvals+lidar4$coeff[3]*xvals*xvals+
      lidar4$coeff[4]*xvals*xvals*xvals+lidar4$coeff[5]*xvals*xvals*xvals*xvals,
      lty=1)
dev.off()
#
# Fig 11.1(d)
#
pdf("lidardeg8fig.pdf",height=4.5,width=4.5)
plot(logratio~range,xlab="Range (m)",ylab="Log Ratio",col="grey")
lines(xvals,lidar8$coeff[1]+lidar8$coeff[2]*xvals+lidar8$coeff[3]*xvals^2+
      lidar8$coeff[4]*xvals^3+lidar8$coeff[5]*xvals^4+lidar8$coeff[6]*xvals^5+
      lidar8$coeff[7]*xvals^6+lidar8$coeff[8]*xvals^7+lidar8$coeff[9]*xvals^8,
      lty=1)
dev.off()
#
# Illustration of splines with 2 knots for the lidar data
#
y <- logratio
x <- range
xi1 <- min(x) + (range(x)[2]-range(x)[1])/3
xi2 <- min(x) + 2*(range(x)[2]-range(x)[1])/3
# Horizontal lines
xtruncflat1 <- ifelse(x>xi1,1,0)
xtruncflat2 <- ifelse(x>xi2,1,0)
pieceflat <- lm(y~1+xtruncflat1+xtruncflat2)
#
# Fig 11.2(a)
#
pdf("lidarpiece0.pdf",height=4.5,width=4.5)
plot(x,y,ylim=c(min(y),max(y)),ylab="y",type="n")
points(x,y,col="grey")
segments(x0=min(x),x1=xi1,y0=pieceflat$coeff[1],y1=pieceflat$coeff[1])
segments(x0=xi1,x1=xi2,y0=pieceflat$coeff[1]+pieceflat$coeff[2],
         y1=pieceflat$coeff[1]+pieceflat$coeff[2])
segments(x0=xi2,x1=max(x),y0=pieceflat$coeff[1]+pieceflat$coeff[2]+
         pieceflat$coeff[3],y1=pieceflat$coeff[1]+pieceflat$coeff[2]+
         pieceflat$coeff[3])
abline(v=xi1,lty=2)
abline(v=xi2,lty=2)
dev.off()
#
xtrunclin1 <- ifelse(x>xi1,x-xi1,0)
xtrunclin2 <- ifelse(x>xi2,x-xi2,0)
piecelin <- lm(y~1+x+xtrunclin1+xtrunclin2)
#
# Fig 11.2(b)
#
pdf("lidarpiece1.pdf",height=4.5,width=4.5)
plot(x,predict(piecelin),type="l",ylim=c(min(y),max(y)),ylab="y")
points(x,y,col="grey")
abline(v=xi1,lty=2)
abline(v=xi2,lty=2)
dev.off()
# Now for quadratic
x2 <- x*x
xtruncquad1 <- ifelse(x>xi1,(x-xi1)^2,0)
xtruncquad2 <- ifelse(x>xi2,(x-xi2)^2,0)
piecequad <- lm(y~1+x+x2+xtruncquad1+xtruncquad2)
#
# Fig 11.2(c)
#
pdf("lidarpiece2.pdf",height=4.5,width=4.5)
plot(x,predict(piecequad),type="l",ylim=c(min(y),max(y)),ylab="y")
points(x,y,col="grey")
abline(v=xi1,lty=2)
abline(v=xi2,lty=2)
dev.off()
# Now for cubic
x3 <- x2*x
xtrunccub1 <- ifelse(x>xi1,(x-xi1)^3,0)
xtrunccub2 <- ifelse(x>xi2,(x-xi2)^3,0)
piececub <- lm(y~1+x+x2+x3+xtrunccub1+xtrunccub2)
#
# Fig 11.2(d)
#
pdf("lidarpiece3.pdf",height=4.5,width=4.5)
plot(x,predict(piececub),type="l",ylim=c(min(y),max(y)),ylab="y")
points(x,y,col="grey")
abline(v=xi1,lty=2)
abline(v=xi2,lty=2)
dev.off()
#
# Piecewise splines
#
xi1 <- 1/3
xi2 <- 2/3
x <- y <- seq(0,1,.01)
#
# Fig 11.3
#
pdf("piecewiselin.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1),pty="s")
plot(x,y,type="l",lty=1,ylab="f(x)")
lines(c(0,1),c(1,1),lty=1)
lines(c(xi1,1),c(0,2*xi1),lty=2)
lines(c(xi2,1),c(0,xi1),lty=2)
text(.31,.07,expression(xi[1]))
lines(x=c(.33,.33),y=c(-.01,.01))
text(.64,.07,expression(xi[2]))
lines(x=c(.66,.66),y=c(-.01,.01))
dev.off()
#
# Basis functions for piecewise cubic
#
# Fig 11.4(a)
#
pdf("piecewisecuball1.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1),pty="s")
plot(x,y,type="n",lty=1,ylab="f(x)")
lines(c(0,1),c(1,1),lty=1)
lines(x,y,lty=1)
lines(seq(0,1,.01),seq(0,1,.01)^2,lty=1)
lines(seq(0,1,.01),seq(0,1,.01)^3,lty=1)
dev.off()
beta3 <- 3.3
beta4 <- 12
x3 <- ifelse(x>xi1,(x-xi1)^3,0)
x4 <- ifelse(x>xi2,(x-xi2)^3,0)
#
# Fig 11.4(b)
#
pdf("piecewisecuball2.pdf",height=4.5,width=4.5)
par(mfrow=c(1,1),pty="s")
plot(x,y,type="n",lty=1,ylab="f(x)")
lines(x,x3*beta3,lty=2)
text(.31,.07,expression(xi[1]))
lines(x=c(.33,.33),y=c(-.01,.01))
text(.69,.07,expression(xi[2]))
lines(x=c(.66,.66),y=c(-.01,.01))
lines(x,x4*beta4,lty=2)
text(.31,.07,expression(xi[1]))
lines(x=c(.33,.33),y=c(-.01,.01))
text(.69,.07,expression(xi[2]))
lines(x=c(.66,.66),y=c(-.01,.01))
dev.off()
#
# Lidar with natural smoothing spline and CV
#
library(pspline)
library(SemiPar)
data(lidar)
attach(lidar)
dfval <- seq(4,20,.1)
ocvval <- gcvval <- rocvval <- rgcvval <- NULL
for (i in 1:length(dfval)){	
    origmod <- smooth.Pspline(x=range,y=logratio,method=2,df=dfval[i])
    ocvval[i] <- origmod$cv
    gcvval[i] <- origmod$gcv
}
#
# Fig 11.5
#
pdf("lidarnatsmoothcvplot.pdf",height=4.5,width=4.5)
plot(ocvval~dfval,type="l",lty=1,ylab="Cross-Validation Score",
     xlab="Effective Degrees of Freedom")
lines(gcvval~dfval,lty=2)
legend("topright",legend=c("Ordinary Cross-Validation",
       "Generalized Cross-Validation"),lty=1:2,bty="n")
dev.off()
dfminocvorig <- dfval[which.min(ocvval)]
dfmingcvorig <- dfval[which.min(gcvval)]
cat("Min DF for original and OCV = ",dfminocvorig,"\n") # 9.3
cat("Min DF for original and GCV = ",dfmingcvorig,"\n") # 9.4
#
# Lidar with smoothing splines and mixed model
#
library(SemiPar)
data(lidar)
attach(lidar)
knots.pos <- seq(400,700,length=20)
fit <- spm(lidar$logratio~f(lidar$range,basis="trunc.poly",
       degree=3,knots=knots.pos))          
summary(fit)
par(mfrow=c(1,1))
#
# Fig 11.6
#
pdf("lidarmixed.pdf",height=6,width=7)
plot(logratio~range,col="grey",xlab="Range (m)",ylab="Log Ratio")
lines(sm.spline(x=range,y=logratio,df=dfmingcvorig),lty=1)
lines(x=lidar$range,y=fit$fit$fitted,lty=2)
legend("bottomleft",legend=c("GCV natural cubic spline","Mixed model cubic spline"),lty=1:2,bty="n")
dev.off()
#
# B-splines figure
#
library(splines)
x <- seq(0,1,.01)
K <- 9 # K=no of knots, so basis is of dimension K+4 (inc intercept).
knots <- seq(1,K)/(K+1)
y <- bs(x,knots=knots,intercept=T)
nbasis <- K+4
#
# Fig 11.7
#
pdf("Bsplinefig.pdf",h=6,w=4.5)
plot(y[,1]~x,type="l",xlab="x",ylab="B-spline",ylim=c(0,max(y)))
for (j in 2:nbasis){
    lines(y[,j]~x,lty=j)
}
points(x=knots,y=rep(0,K))
dev.off()
#
# LIDAR example at the end of Section 11.2.7
#
#
# Tatiana Maravina's code for carrying out confidence band calculation.
#
# First get function that computes knots locations
#
source("http://matt-wand.utsacademics.info/webspr/default.knots.r")
#
#    Fit using mixed models representations
#    for later comparison with my 'by hand' results
#
library(SemiPar)
data(lidar)
attach(lidar)
fit=spm(logratio~f(range, basis="trunc.poly", degree=3))   
plot(fit)
points(range, logratio, pch=20)
#
# 'By hand' fitting of penalized splines and computation of CIs
# Define some functions first
#
Tprimenorm <- function (x, knots, W, p) {
# computes ||T'(x)|| - function to be integrated to get kappa0
# x can be a vector
	K <- length(knots)
	m <- length(x)
	S <- Sprime <- matrix(nrow=m, ncol=1+p+K)
	for (i in 1:m) {
		S[i,] <- c(x[i]^(0:p), ((x[i]-knots)^p)*(x[i]>knots))
		Sprime[i,] <- c(0, 1:p*(x[i]^(0:(p-1))), 
                           (3*(x[i]-knots)^2)*(x[i]>knots))
	}
	S <- S%*%W		 # m-by-n 
	Sprime <- Sprime%*%W	 # m-by-n
	
	Snorm2 <- apply(S^2, 1, sum) 
	numer <- sqrt(Snorm2*apply(Sprime^2, 1, sum)-
                 apply((S*Sprime)^2,1,sum))  # m-by-1
	return(numer/Snorm2)
}
#
f <- function (x, kappa0, alpha) {
	# solution of f(x)=0 gives c
	2*(1-pnorm(x))+kappa0*exp(-x^2/2)/pi-alpha
}
#
## data
#
x <- lidar$range
y <- lidar$logratio
## Standardize x 
# Otherwise, I kept running into "system is computationally singular" problem
#    when inverting (t(X)%*%X + lambda*D)
mu <- mean(x)
sig <- sd(x)
xs <- (x-mu)/sig
#
n <- length(y)
K <- 35 # number of knots
knots <- default.knots(xs, K)   # knots
p <- 3  # order (3 for cubic)
X <- (matrix(nrow=n, ncol=K, data=xs)-matrix(nrow=n, ncol=K, data=knots, 
     byrow=TRUE))
X <- (X^p)*(X>0)  # use truncated power basis of order p
Xmat <- cbind(model.matrix(~xs+I(xs^2)+I(xs^3)), X)  # EDIT TO MATCH ORDER!
D <- diag(c(rep(0,p+1), rep(1,K)))
#
## Find lambda that minimizes GCV
#
lambda <- seq(0,1, by=0.01)
N <- length(lambda)
GCV <- df <- numeric(N)
for (i in 1:N) {
# following the Semiparametric Regression book and 
# http://www.uow.edu.au/~mwand/SPmanu.pdf
	#  the multiplier is lambda^(2*p)
	W <- solve(t(Xmat)%*%Xmat+(lambda[i]^(2*p))*D)%*%t(Xmat)
	beta.hat <- W%*%y
	y.hat <- Xmat%*%beta.hat
	S <- Xmat%*%W
	Sii <- diag(S)
	df[i] <- sum(Sii)
	GCV[i] <- sum((y-y.hat)^2)*n/(n-df[i])^2
}
plot(lambda, GCV, type="l", xlab=expression(lambda))
lambdaopt <- lambda[which.min(GCV)]  # optimal lambda is 0.41
abline(v=lambdaopt, lty=2, col="darkgrey")
#
## Fit model with optimal value of lambda
#
W <- solve(t(Xmat)%*%Xmat+(lambdaopt^(2*p))*D)%*%t(Xmat)
beta.hat <- W%*%y
y.hat <- Xmat%*%beta.hat
## Find kappa0 and c
a <- min(xs)
b <- max(xs)
kappa0 <- integrate(Tprimenorm, a, b, knots=knots, W=W, p=p)$value
kappa0	# 15.44
cc <- uniroot(f, c(0, 10), kappa0=19, alpha=0.05)$root # c=3.11 as in text
## Variance estimation
# Clearly, the data are heteroscedastic
resids <- y-y.hat
z <- log(resids^2)
plot(x, z, pch=20)  
# looks linear so I will simply use lm() instead of nonparametric methods
# q=spm(z~f(x))$fit$fitted
q <- lm(z~x)
# evaluate approximate standard errors at the observed range values
sigma.e <- sqrt(exp(predict(q)))
S <- Xmat%*%W
Snorm <- sqrt(apply(S^2, 1, sum))
se <- sigma.e*Snorm
# what if assumed constant variance?
v <- sum(diag(S))  		# trace(S)
vv <- sum(diag(t(S)%*%S))  # trace(t(S)%*%S)
sigma.e.homosc <- sqrt(sum(resids^2)/(n-2*v-vv))
sigma.e.homosc   # compare with fit$fit$sigma   0.083 vs 0.080
se.homosc <- sigma.e.homosc*Snorm  
# compare with predict(fit, newdata=lidar, se=TRUE)$se
#
# Confidence bands
#
# Fig 11.8(a)
#
pdf("lidarsim2.pdf",h=6,w=6)
par(mfrow=c(1,1))
plot(x, y, pch=20, xlab="Range (m)", ylab="Log Ratio",col="grey") 
lines(x, y.hat, lwd=2) 
lines(x, y.hat-cc*se, lty=2)
lines(x, y.hat+cc*se, lty=2)
lines(x, y.hat-1.96*se, lty=3)
lines(x, y.hat+1.96*se, lty=3)
legend("bottomleft", legend=c("Pointwise", "Simultaneous"),lty=3:2, ,bty="n")
dev.off()
#
# Fig 11.8(b)
#
pdf("lidarsim1.pdf",h=6,w=6)
plot(x, y, pch=20, xlab="Range (m)", ylab="Log Ratio", col="grey")
lines(x, y.hat,  lwd=2) 
lines(x, y.hat-cc*se.homosc, lty=2)
lines(x, y.hat+cc*se.homosc, lty=2)
lines(x, y.hat-1.96*se.homosc, lty=3)
lines(x, y.hat+1.96*se.homosc, lty=3)
legend("bottomleft", legend=c("Pointwise", "Simultaneous"), lty=3:2, bty="n")
dev.off()
#
# Mixed model LIDAR spline example at the end of Section 11.2.7
#
knots.pos <- seq(400,700,length=20)
fit <- spm(lidar$logratio~f(lidar$range,basis="trunc.poly",
       degree=3,knots=knots.pos))
xseq <- seq(400,700,1)
cubfit <- fit$fit$coef$fixed[1] + fit$fit$coef$fixed[2]*xseq + 
          fit$fit$coef$fixed[3]*xseq^2 +fit$fit$coef$fixed[4]*xseq^3
xbas <- matrix(0,nrow=20,ncol=301)
# totfit <- cubfit+randomfit
randomfit <- 0
rand <- matrix(0,nrow=20,ncol=301)
figname <- c("lidars1.pdf","lidars2.pdf","lidars3.pdf","lidars4.pdf","lidars5.pdf","lidars6.pdf","lidars7.pdf","lidars8.pdf","lidars9.pdf","lidars10.pdf","lidars11.pdf","lidars12.pdf","lidars13.pdf","lidars14.pdf","lidars15.pdf","lidars16.pdf","lidars17.pdf","lidars18.pdf","lidars19.pdf","lidars20.pdf")
caps <- c("1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th","11th","12th","13th","14th","15th","16th","17th","18th","19th","20th")
for (i in 1:20){
    xbas[i,] <- ifelse(xseq > knots.pos[i],(xseq-knots.pos[i])^3,0)
    rand[i,] <- xbas[i,]*fit$fit$coef$random[i]
# Fig 11.9 -- create 20 different files
    pdf(figname[i],height=5,width=5)
    par(mfrow=c(1,1),pty="s")
    plot(logratio~range,col="grey",xlab="",ylab="",ylim=c(-1.5,1.5),
    main=caps[i],cex.main=2)
    lines(xseq,cubfit,lty=1)
    lines(xseq,rand[i,],lty=2)
    randomfit <- randomfit + xbas[i,]*fit$fit$coef$random[i]
    abline(v=knots.pos[i],lty=3)
#    totfit <- cubfit+randomfit
#    lines(xseq,totfit,lty=3)
    dev.off()
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Putting prior on the effective degrees of freedom in Section 11.2.9
# and then analyzing the SBMD data -- Code written by Youyi Fong
#--------------------------------------------------------------------
# Code to evaluate induced prior distribution on the degrees of freedom
#
DFfun <- function(m,n,Cmat,sigmae2,sigmab2){
      Lambda <- diag(c(0,rep(sigmae2/sigmab2,m)))
      Totdf <- tr( solve(crossprod(Cmat) + Lambda) %*% crossprod(Cmat) )
}
#
create.matrix=function (n, m){
    out = matrix(0,n*m,m)
    for (i in 1:m)
        out[(i-1)*n+1:n,i]=1
    out
}
#
# Simple one-way ANOVA, example at end of Section 11.2.9
#
m <- 10 # Number of groups
n <- 5  # Number of observations per group
#
Xmat <- matrix(1,nrow=m*n,ncol=1)
Zmat <- create.matrix(n,m)
tr <- function(M) sum(diag(M))
Cmat <- cbind(Xmat,Zmat)
sigmae2 <- 1
nsim <- 10000
sigmab2 <- 1/rgamma(nsim,.5,.005)
par(mfrow=c(1,1))
cuto <- 2.5
#
# Fig 11.10(a)
#
pdf("DFPrior1.pdf",w=4,h=4) 
hist(sqrt(sigmab2[sqrt(sigmab2)<cuto]),nclass=30,xlim=c(0,cuto),main="",
     xlab=expression(sigma[b]))
dev.off()
DFvals <- DFvalsred <- NULL
count <- 0
big <- 10^10
for (i in 1:length(sigmab2)){
    if (sigmab2[i]<big) {
       count <- count+1
       DFvalsred[i] <- DFfun(m,n,Cmat,sigmae2,sigmab2[i])
    }
    else DFvalsred[i] <- m
}
cat("Count = ",count,"\n")
cat("Degrees of freedom above 9.9 = ",length(DFvalsred[DFvalsred>9.9])/nsim,
     "\n")
cat("Proportion below cut = ",length(sigmab2[sqrt(sigmab2)<cuto])/nsim,"\n")
#
# Fig 11.10(b)
#
pdf("DFPrior2.pdf",h=4,w=4) 
hist(DFvalsred,main="",xlab="Effective Degrees of Freedom")
dev.off()
#
# Fig 11.10(c)
#
pdf("DFprior3.pdf",h=5,w=5) 
maxs <- max(sqrt(sigmab2[sqrt(sigmab2)<cuto]))
plot(sqrt(sigmab2[sqrt(sigmab2)<cuto]),
     DFvalsred[sqrt(sigmab2)<cuto],ylim=c(1,m),xlim=c(0,maxs),
ylab="Effective Degrees of Freedom",xlab=expression(sigma[b]),main="",type="n")
lines(sort(sqrt(sigmab2[sqrt(sigmab2)<cuto])),
      sort(DFvalsred[sqrt(sigmab2)<cuto]))
abline(h=m,lty=2,col="black")
abline(h=1,lty=2,col="black")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spinal bone minreral density example at the end of Section 11.2.9
#
GammaPriorCh <- function(theta,prob,d){
             a <- d/2
             b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
             cat("Gamma Parameters: ",a,b,"\n")
             list(a=a,b=b)
}
# INLA download
source("http://www.math.ntnu.no/inla/givemeINLA.R")
library(INLA)
library(mvtnorm)
library(pixmap)
require(splines)
require(nlme)
bone <- read.table("spinalbonedata.txt", header=T)
bone <- subset(bone, sex=="fem")
bone$ethnic <- as.factor(bone$ethnic)
bone$groupVec <- 1
bone$seqno <- 1:nrow(bone)
#
# Matt Wand code
#
formOmega <- function(a,b,intKnots)
{
   allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- t(Bdd*wts)%*%Bdd     
   return(Omega)
}
# Obtain the spline component of the Z matrix:
a <- 8 ; b <- 28; 
numIntKnots <- 15 
intKnots <- quantile(unique(bone$age),
            seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
Omega <- formOmega(a,b,intKnots)
eigOmega <- eigen(Omega)
indsZ <- 1:(numIntKnots+2)
UZ <- eigOmega$vectors[,indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
B <- bs(bone$age,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
ZSpline <- B%*%LZ   
ZBlock <- list(list(groupVec=pdIdent(~ZSpline-1)),list(idnum=pdIdent(~1)))
ZBlock <- unlist(ZBlock,recursive=FALSE)
bone.inla.fit <- inla(spnbmd ~ ethnic + age + f(idnum,model="iid",
                param=c(.5, 5e-6),diagonal=0) + f(seqno,model="z",Z=ZSpline,
                initial=3,param=c(.5, 0.00113)), data=bone, family="gaussian", 
                control.predictor=list(compute=TRUE) )
bone.hyperpar <- inla.hyperpar(bone.inla.fit)
summary(bone.inla.fit)
#
# C matrix
X  <- model.matrix( ~ age + ethnic, data=bone)
X1 <- model.matrix( ~ -1 + as.factor(idnum), data=bone)
C <- cbind (X, ZSpline, X1)
#
# return a vector of 4 dofs: total, fixed effects, spline, and
# individual random effects
#
get.dof = function(sigma2.epsilon, sigma2.1, sigma2.2) {
    # colInd specifies the columns for each group
    colInd1 = 1:ncol(X)
    colInd2 = ncol(X)+1:ncol(ZSpline)
    colInd3 = ncol(X)+ncol(ZSpline)+1:ncol(X1)
# cat(length(colInd3),"\n")
    # Lambda matrix
    Lambda = matrix(0, ncol(C), ncol(C))
    diag(Lambda)[colInd2] = sigma2.epsilon/sigma2.1
    diag(Lambda)[colInd3] = sigma2.epsilon/sigma2.2
#  cat(diag(Lambda)[colInd3],"\n\n")

#  cat("Dim of Lambda = ",dim(Lambda),"\n")
# return degrees of freedom
    tmp = diag ( solve(t(C) %*% C + Lambda) %*% t(C) %*% C )
    c(sum(tmp), sum(tmp[colInd1]), sum(tmp[colInd2]), sum(tmp[colInd3]))
}
#
nsim <- 1000
sigma2.1 <- 1/rgamma(nsim,.5,5e-06)
DFvals <- DFvalsred <- NULL
count <- 0
big <- 10^10
sig.eps <- 0.03276967
for (i in 1:length(sigma2.1)){
    if (sigma2.1[i]<big) {
       count <- count+1
       DFvalsred[i] <- get.dof(sig.eps^2,sigma2.1[i],.00001)[3]
    }
    else DFvalsred[i] <- 17
}
cat("Count = ",count,"\n")
sigma2ch <- GammaPriorCh(.3,.95,1)
nsim <- 1000
sigma2.2 <- 1/rgamma(nsim,sigma2ch$a,sigma2ch$b)
cuto <- 2
DFvals2 <- NULL
count <- 0
for (i in 1:length(sigma2.2)){
    if (sigma2.2[i]<big) {
       count <- count+1
       DFvals2[i] <- get.dof(sig.eps^2,.00001,sigma2.2[i])[4]
    }
    else DFvals2[i] <- 226
}
cat("Count = ",count,"\n")
#
# Construct the model fit figure
#
library(lattice)
lattice.options(default.theme = "col.whitebg")
#
# Fig 11.11
#
pdf("spinalbonev2.pdf") 
bone$fitted1 <- bone.inla.fit$summary.fitted.values$mean
bone$spline <- ZSpline %*% sapply(bone.inla.fit$summary.random[-1],
       function(obj) obj$mean)  + X %*%
bone.inla.fit$summary.fixed[colnames(X),1]
xyplot(spnbmd~age|ethnic, data = bone, groups=idnum, xlab="Age (years)",
  ylab="Spinal BMD",layout=c(2,2),aspect=1, panel=function(x,y,subscripts,groups) {
        lpoints (x,y, pch=20, col=1, cex=.3)
        ids=unique(groups[subscripts])
        for (id in ids) {
            subset=groups[subscripts]==id
            llines(x[subset], y[subset], col="darkgrey", lty=1)
            llines(x[subset], bone$fitted1[subscripts[subset]], col=1,lty=2)
        }
# spline curve
        llines(x[order(x)], bone$spline[subscripts][order(x)], col=1,lwd=4)
    }
      )
dev.off()
# 
# Prior and posterior picture
#
# Fig 11.12(a)
#
pdf("splinepriorfig1.pdf")
par(mfrow=c(1,1))
scuto <- .1
plot1 <- sqrt(sigma2.1)[sqrt(sigma2.1)<scuto]
hist( plot1, nclass=30,xlim=c(0,max(plot1)),main="",xlab=expression(sigma[1]))
abline(v=0.018,lty=2,col="black") # Post est
dev.off()
#
# Fig 11.12(b)
#
pdf("splinepriorfig2.pdf")
hist(DFvalsred,main="",xlab="Effective Degrees of Freedom")
abline(v=get.dof(0.033^2,0.018^2,0.109^2)[3],lty=2,col="black")
dev.off()
#
# Fig 11.12(c)
#
pdf("splinepriorfig3.pdf")
plot(plot1,DFvalsred[sqrt(sigma2.1)<scuto],ylim=c(1,17),ylab="Effective Degrees of Freedom",xlab=expression(sigma[1]),main="",type="n")
lines(sort(plot1),sort(DFvalsred[sqrt(sigma2.1)<scuto]))
abline(h=17,lty=2,col="black")
abline(h=1,lty=2,col="black")
scuto <- sqrt(.8)
plot2 <- sqrt(sigma2.2)[sqrt(sigma2.2)<scuto]
dev.off()
#
# Fig 11.12(d)
#
pdf("splinepriorfig4.pdf")
hist(plot2,nclass=30,xlim=c(0,max(plot2)),main="",xlab=expression(sigma[2]))
abline(v=0.109,lty=2,col="black") # Post est
dev.off()
#
# Fig 11.12(e)
#
pdf("splinepriorfig5.pdf")
hist(DFvals2,main="",xlab="Effective Degrees of Freedom")
abline(v=get.dof(0.033^2,0.018^2,0.109^2)[4],lty=2,col="black")
dev.off()
#
# Fig 11.12(f)
#
pdf("splinepriorfig6.pdf")
plot(plot2,DFvals2[sqrt(sigma2.2)<scuto],ylim=c(1,226),
  ylab="Effective Degrees of Freedom",xlab=expression(sigma[2]),main="",type="n")
lines(sort(plot2),sort(DFvals2[sqrt(sigma2.2)<scuto]))
abline(h=226,lty=2,col="black")
abline(h=1,lty=2,col="black")
dev.off()
#
# REML fit for Table 11.1
#
lme.model = spnbmd ~ age + ethnic
remlfit <- lme( lme.model, random=ZBlock, data=bone )
summary(remlfit)
#
# INLA fit for Table 11.1
#
summary(bone.inla.fit)
m1 <- bone.inla.fit$marginals.hyperpar$`Precision for seqno.1`
m2 <- bone.inla.fit$marginals.hyperpar$`Precision for idnum`
m3 <- bone.inla.fit$marginals.hyperpar$`Precision for the Gaussian observations`
postmeansigma1 <- inla.emarginal(function(x) 1/sqrt(x), m1)
temp1 <- inla.emarginal(function(x) 1/x, m1)
temp2 <- inla.emarginal(function(x) 1/x, m2)
temp3 <- inla.emarginal(function(x) 1/x, m3)
postmeansigma2 <- inla.emarginal(function(x) 1/sqrt(x), m2)
postmeansigmae <- inla.emarginal(function(x) 1/sqrt(x), m3)
postsdsigma1 <- sqrt(temp1-postmeansigma1^2)
postsdsigma2 <- sqrt(temp2-postmeansigma2^2)
postsdsigmae <- sqrt(temp3-postmeansigmae^2)
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Kernel Section 11.3
#
#  Fig 11.13
#
pdf("kernelfig.pdf",height=5.5,width=7.5)
par(mfrow=c(1,1))
x <- seq(-3,3,.01)
xs <- seq(-1,1,.01)
plot(x,dnorm(x,0,1),type="l",lty=1,ylim=c(0,.9),ylab="Kernel Density")
lines(xs,.75*(1-xs**2),lty=2,lwd=1)
lines(xs,(70/81)*(1-abs(xs)**3)**3,lty=3,lwd=1)
lines(x=c(-1,-1),y=c(0,.5),lty=4,lwd=1)
lines(x=c(-1,1),y=c(0.5,.5),lty=4,lwd=1)
lines(x=c(1,1),y=c(0.5,0),lty=4,lwd=1)
legend("topleft",legend=c("Gaussian","Epanechnikov","Tricube","Boxcar"),
       lty=c(1,2,3,4),lwd=1,bty="n")
dev.off()
#
# Local polynomials
#
library(locfit)
gcv1 <- gcvplot(logratio~range,kern="gauss",alpha=seq(.1,.9,by=.01),deg=1)
plot(gcv1)
gcv2 <- gcvplot(logratio~range,kern="tcub",alpha=seq(.1,.9,by=.01),deg=1)
plot(gcv2)
gcv3 <- gcvplot(logratio~range,kern="epan",alpha=seq(.1,.9,by=.01),deg=1)
plot(gcv3)
locfit1 <- locfit(logratio~range,kern="gauss",deg=1,
           alpha=gcv1$alpha[which(gcv1$values==min(gcv1$values))])
locfit2 <- locfit(logratio~range,kern="tcub",deg=1,
           alpha=gcv2$alpha[which(gcv2$values==min(gcv2$values))])
locfit3 <- locfit(logratio~range,kern="epan",deg=1,
           alpha=gcv3$alpha[which(gcv3$values==min(gcv3$values))])
#
# Fig 11.14
#
pdf("lidarlocpolyfig1.pdf",height=4.5,width=4.5)
plot(locfit2,band="none",xlab="Range (m)",ylab="Log Ratio",
     ylim=c(min(logratio),max(logratio)))
points(logratio~range,col="grey")
plot(locfit2,band="none",add=T,lty=2)
plot(locfit3,band="none",add=T,lty=3)
legend("bottomleft",legend=c("Gaussian","Tricube","Epanechnikov"),
       lty=c(1,2,3),bty="n")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Lidar variance estimation -- end of Section 11.4
#
library(SemiPar)
data(lidar)
attach(lidar)
n <- length(logratio)
firstdiffs <- (logratio[2]-logratio[1])^2
ai <- bi <- ci <- epsi <- epsicheck <- NULL
ci[1] <- epsi[1] <- epsicheck[1] <- 0
for (i in 2:(n-1)){
    firstdiffs <- firstdiffs + (logratio[i+1]-logratio[i])^2
    ai[i] <- (range[i+1]-range[i])/(range[i+1]-range[i-1])
    bi[i] <- (range[i]-range[i-1])/(range[i+1]-range[i-1])
    ci[i] <- sqrt(1/(ai[i]^2+bi[i]^2+1))
    epsi[i] <- ai[i]*logratio[i-1]+bi[i]*logratio[i+1]-logratio[i]
}
#
# Naive estimate uses 
naivesd <- sqrt(sum((logratio-fitted(fit))^2)/(n-fit$aux$df.fit))
firstdiffsd <- sqrt(0.5*firstdiffs/(n-1))
secondiffsd <- sqrt(sum(ci^2*epsi^2)/(n-2))
cat("Naive sd = ",naivesd,"\n") 
cat("First diff error sd = ",firstdiffsd,"\n") 
cat("Second diff error sd = ",secondiffsd,"\n") 
#
# Heterogeneous
#
library(pspline)
optmod <- smooth.Pspline(x=range,y=logratio,method=2,df=dfmingcvorig)
resids <- logratio - predict(optmod)$ysmth
#
# Fig 11.15(a)
#
pdf("lidarresid.pdf",height=5,width=5)
plot(resids~range,ylim=c(-.3,.3),xlab="Range (m)",ylab="Residuals")
abline(0,0,lty=2)
dev.off()
z <- log(resids^2)
#
# Fig 11.15(b)
#
pdf("lidarlogresid.pdf",h=5,w=5)
plot(z~range,xlab="Range (m)",ylab="Log Squared Residuals")
zlinmod <- lm(z~range)
abline(zlinmod$coeff)
dev.off()
standresids <- resids/sqrt(exp(predict(zlinmod)))
#
# Fig 11.15(d)
#
pdf("lidarstandresid.pdf",h=5,w=5)
plot(standresids~range,ylim=c(-max(abs(standresids)),max(abs(standresids))),
    xlab="Range (m)",ylab="Standardized Residuals")
abline(0,0,lty=2)
dev.off()
rvals <- seq(min(range),max(range),1)
sdvals <- sqrt(exp(zlinmod$coeff[1]+zlinmod$coeff[2]*rvals))
#
# Fig 11.15(c)
#
pdf("lidarsdest.pdf",h=5,w=5)
plot(sdvals~rvals,xlab="Range (m)",ylab="Estimated Standard Deviation",type="n")
lines(sdvals~rvals)
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BPD example in Section 11.5
#
library(locfit)
myexpit <- function(x,b0,b1){expit <- exp(b0+b1*x)/( 1+exp(b0+b1*x) )} 
bw <- read.table("birthweight.txt",header=T)
attach(bw)
lrmod1 <- glm(BPD~birthweight,family=binomial)
x <- seq(min(birthweight),max(birthweight))
y <- bw$BPD
meanbw <- mean(bw$birthweight)
sdbw <- sd(bw$birthweight)
x <- (bw$birthweight-meanbw)/sdbw
#
# My own AIC function
#
aicbinfun <- function(y,n,p,df){
	  ndat <- length(y)
	  dfterm <- 2*df
	  loglik <- 0
	  for (i in 1:ndat){
	      loglik <- loglik + lchoose(n[i],y[i]) + y[i]*log(p[i]) + 
                        (n[i]-y[i])*log(1-p[i])
	  }
	  aicbin <- (-2*loglik  + dfterm)
	  list(aicbin=aicbin,dfterm=dfterm,loglik=loglik)
}
n <- rep(1,223)
hval <- seq(1,5,.1)
aicval <- myaicval <- df <- loglik <- NULL
for (i in 1:length(hval)){
    locfitbw <- locfit(y~lp(x,deg=1,nn=0,h=hval[i]),data=bw,family="binomial")	
    preds <- predict(locfitbw,newdata=x)
    df[i] <- locfitbw[9]$dp[6]
    myaicval[i] <- aicbinfun(y,n,preds,df[i])$aicbin
    loglik[i] <- aicbinfun(y,n,preds,df[i])$loglik
}
par(mfrow=c(1,2))
plot(myaicval~hval)
plot(df~hval)
hopt <- hval[which.min(myaicval)]
dfopt <- df[which.min(myaicval)]
cat("Minimizing values h and df =",hopt,dfopt,"\n") # df=4.1
#
# Now obtain the fit with the minimizing value of h
#
locgcv <- locfit(y~lp(x,deg=1,nn=0,h=hopt),data=bw,family="binomial")
xseq <- seq(min(x),max(x),.01)
preds <- predict(locgcv,newdata=xseq)
newx <- xseq*sdbw + meanbw
#
# Penalized regression splines for the BPD data
#
library(mgcv)
attach(bw)
# 
gammod <- gam(bw[,2]~s(x=bw$birthweight,k=10,fx=F,bs="cr",m=2),family=binomial)
orderbw <- cbind(BPD,birthweight)
orderbw <- orderbw[order(orderbw[,2]),]
orderfit <- cbind(birthweight,gammod$fitted)
orderfit <- orderfit[order(orderfit[,1]),]
#
# Fig 11.16
#
pdf("BPDcubic.pdf",width=4.5,height=4.5)  
par(mfrow=c(1,1))
plot(birthweight,BPD,pch="|",xlab="Birthweight (grams)")
lines(orderfit[,1],orderfit[,2])
lrmod1 <- glm(BPD~birthweight,family=binomial)
xseq2 <- seq(min(birthweight),max(birthweight))
myexpit <- function(x,b0,b1){expit <- exp(b0+b1*x)/( 1+exp(b0+b1*x) )} 
lines(xseq2,myexpit(xseq2,b0=lrmod1$coeff[1],b1=lrmod1$coeff[2]),lty=3)
lines(newx,preds,lty=2)
legend(x=1000,y=.95,legend=c("Cubic Spline","Local Likelihood","Linear Logistic"),
       bty="n",lty=1:3)
dev.off()
