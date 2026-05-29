#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 9
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Contraception data
#
library(xtable)
alldat <- read.table("contracept.txt",header=F,na.strings=".",
          col.names=c("id","dose","time","y"))
#
# first plot averages over time
y <- dosemat <- idmat <- timemat <- matrix(0,nrow=1151,ncol=4)
dose <- NULL
y0 <- matrix(0,nrow=576,ncol=4)
time0 <- matrix(0,nrow=576,ncol=4)
dose0 <- matrix(0,nrow=576,ncol=4)
id0 <- matrix(0,nrow=576,ncol=4)
y1 <- matrix(0,nrow=575,ncol=4)
time1 <- matrix(0,nrow=575,ncol=4)
dose1 <- matrix(0,nrow=576,ncol=4)
id1 <- matrix(0,nrow=575,ncol=4)
count0 <- count1 <- 1
p00 <- p01 <- p10 <- p11 <- rep(0,4)
counter <- 1
n0 <- n1 <- rep(0,4)
for (i in 1:1151){
    dose[i] <- alldat[counter,2]
# cat("dose counter: ",dose[i],counter,"\n")
    for (j in 1:4){
    	y[i,j] <- alldat[counter,4]
        if (dose[i]==0){
 	   y0[count0,j] <- y[i,j]
	   if (is.na(y[i,j]) == F) {n0[j] <- n0[j]+1}
           if(j==4){count0 <-count0+1}
	}
	if (dose[i]==1){
	   y1[count1,j] <- y[i,j]; if(j==4){count1 <- count1+1}
	   if (is.na(y[i,j]) == F) {n1[j] <- n1[j]+1}
        }
        if (is.na(y[i,j])==F ){
        if (dose[i]==0){
	   if (y[i,j]==0) {p00[j]<-p00[j]+1}
           if (y[i,j]==1) {p01[j]<-p01[j]+1}
        }
	if (dose[i]==1){
           if (y[i,j]==0) {p10[j]<-p10[j]+1}
           if (y[i,j]==1) {p11[j]<-p11[j]+1}
           
        }
        }
    counter <- counter+1
    }
}
alldat0 <- as.table(id0,dose0,y0,time0)
cor0 <- cor(y0,use="pairwise")
print(xtable(cor0))
cor1 <- cor(y1,use="pairwise")
# 
# Creates Tab 9.1
#
print(xtable(cor1))
cov(y0,use="pairwise")
cov(y1,use="pairwise")
for (j in 1:4){
    p00[j] <- p00[j]/(p00[j]+p01[j])
    p01[j] <- 1-p00[j]
    p10[j] <- p10[j]/(p10[j]+p11[j])
    p11[j] <- 1-p10[j]
}
# variance of logits:
vlogit0 <- 1/(n0*p01*(1-p01))
vlogit1 <- 1/(n1*p11*(1-p11))
cat("SEs of logits for low dose: ",sqrt(vlogit0),"\n")
cat("SEs of logits for high dose: ",sqrt(vlogit1),"\n")
#
# Fig 9.1
#
pdf("contracept-prob.pdf",h=4,w=4)
plot(seq(1,4),p01,pch=19,ylim=c(min(p01,p11),max(p01,p11)),
   ylab="Probability of Amenorrhea",xlab="Measurement Occasion",xaxp=c(1,4,3))
points(seq(1,4),p11,pch=21)
legend("topleft",legend=c("High Dose","Low Dose"),bty="n",pch=c(21,19))
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Seizure data
#
# THE FOLLOWING DESCRIBES THE COLUMNS OF  thall5.dat:
#       ID
#       Seizure counts
#       Visit           (visit==0, through visit==4)
#       TX              1 if progabide; 0 otherwise
#       Age             
#       Weeks           # weeks over which seizures were recorded
# REFERENCE:   Thall and Vail (1991) Biometrics
#
seiz <- read.table("seizure.dat",col.names=
        c("ID","y","redun","x1","age","time"))
x2 <- rep(c(0,1,1,1,1),59) # Time effect
x3 <- seiz$x1*x2 
# x3 is the Treatment by period (pre-/post) 
# interaction-covariate of interest
#
newtime <- rep(c(4.5,9,11,13,15))
seiz2 <- cbind(seiz,seiz$x1,x2,x3,newtime)
#
denom <- c(8,2,2,2,2)
visit <- seq(0,4)
#
# Fig 9.2
#
pdf("seizinit.pdf")
plot(rep(visit,59),(seiz2$y+.5)/8,type="n",ylim=c(0.2,90),log="y",
     xlab="Period",ylab="Seizure Rate",cex.lab=1.5)
for (i in 1:59){
    index <- (i-1)*5+seq(1,5)    
    if (i<=30 ) {
       meas0 <- ((seiz2$y[index]+.5)/denom)
       lines(visit,meas0,lty=1)
    }
    else {
       meas1 <- ((seiz2$y[index]+.5)/denom)
       lines(visit,meas1,lty=2)
    }
}
legend("topright",bty="n",legend=c("Placebo","Progabide"),lty=1:2,cex=1.5)
dev.off()
#
ly0mean <- ly1mean <- y0mean <- y1mean <- y0var <- y1var <- rep(0,5)
plactot <- 28*5
progtot <- 31*5
for (j in 1:5){
        ind0 <- seq(j,plactot,5)
        ind1 <- plactot+seq(j,progtot,5)
	y0data <- seiz2$y[ind0]
	y1data <- seiz2$y[ind1]
        ly0data <- log((seiz2$y[ind0]+.5)/denom[j])
	ly1data <- log((seiz2$y[ind1]+.5)/denom[j])
        y0mean[j] <- mean(y0data)
        y1mean[j] <- mean(y1data)
        ly0mean[j] <- mean(ly0data)
        ly1mean[j] <- mean(ly1data)
        y0var[j] <- var(y0data)
        y1var[j] <- var(y1data)
}
#
# Fig 9.3
#
pdf("seizsum.pdf")
plot(visit,y0mean/denom,type="n",xlab="Period",log="y",ylim=c(3.3,4.7),
ylab="Average Seizure Rate",cex.lab=1.5)
# title("Means versus time, 0=placebo,1=progabide")
lines(visit,y0mean/denom,lty=1)
lines(visit,y1mean/denom,lty=2)
legend("bottomleft",legend=c("Placebo","Progabide"),lty=1:2,bty="n",cex=1.5)
dev.off()
#
# Examine the ratio of the variance to the mean, as in Table 9.2
#
y0var/y0mean
y1var/y1mean
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Theophylline
#
library(nlme)
data(Theoph)
newTheoph <- Theoph[Theoph$Time>0,]
#  
# PLOT THE THEOPHYLLINE DATA
# The object Theoph is a groupedData object
#
# Fig 9.4
#
pdf("theoph_fig1.pdf",h=9,w=7)
plot(newTheoph)
dev.off()
#
# Seizure example in Section 9.4
#
attach(seiz2)
library(lme4)
library(MASS)
#
# Use Gauss-Hermite rule with 50 points
#
lmermod3 <- glmer(y ~ x1+x2+x3+(1|ID)+offset(log(time)),
            family=poisson,data=seiz2,nAGQ=50)
#
# PQL is inaccurate here
#
lmermodPQL <- glmmPQL(y ~ x1+x2+x3+offset(log(time)),random=~1|ID,
            family=poisson,data=seiz2)
#
# Laplace approximation
#
lmermod2 <- lmer(y ~ x1+x2+x3+(1|ID)+offset(log(time)),family=poisson,
            data=seiz2,nAGQ=1)
#
# Function to work out the model-based estimate of the correlation for two
# means of sizes mu1 and mu2
#
poiscor <- function(mu1,mu2,kappa){
	kappa*sqrt(mu1*mu2)/sqrt((1+kappa*mu1)*(1+kappa*mu2))
}
sigma2 <- 0.61; kappa <- exp(sigma2)-1
#
# Correlations in seizure example in Section 9.5
#
mu1 <- 1; mu2 <- 1; poiscor(mu1,mu2,kappa)
mu1 <- 2; mu2 <- 2; poiscor(mu1,mu2,kappa)
mu1 <- 5; mu2 <- 5; poiscor(mu1,mu2,kappa)
#
# Conditional likelihood seizure example in Section 9.5
#
seizy <- c(963,987)
seizn <- c(1825,1967)
x <- c(0,1)
binmod <- glm(cbind(seizy,seizn-seizy)~x,family="binomial")
summary(binmod)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bayesian analysis of seizure data (Section 9.6): see WinBUGS code file
# SeizureWinBUGS.txt
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s# Lung cancer and radon example with spatial model (Section 9.7), see
# seperate code
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Seizure data analyzed with GEE at the end of Section 9.11
#
library(geepack)
seiz <- read.table("seizure.dat",col.names=c("ID","y","redun","x1","age","time"))
x2 <- rep(c(0,1,1,1,1),59) # Time effect
x3 <- seiz$x1*x2 # Treatment by period (pre-/post) interaction which
                 # is the covariate of interest
newtime <- rep(c(4.5,9,11,13,15))
seiz2 <- cbind(seiz,x1,x2,x3,newtime)
attach(seiz2)
#
# Table 9.7, GEE independence results for seizure data
#
modgee1 <- geese(y ~ x1+x2+x3+offset(log(time)),id=ID,family=poisson,data=seiz,
           corstr="independence")
#
# Table 9.7, GEE exchangeable results for seizure data
#
modgee2 <- geese(y ~ x1+x2+x3+offset(log(time)), id=ID,family=poisson,data=seiz,
           corstr="exchangeable")
#
# Table 9.7, GEE AR(1) results for seizure data
#
modgee3 <- geese(y ~ x1+x2+x3+offset(log(time)), id=ID,family=poisson,data=seiz,
           corstr="ar1")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Illustrate attenuation in Section 9.13
#
beta0 <- 0; beta1 <- 1; sigma0 <- 2; x <- seq(-4,4,.1)
K <- 40 # Number of curves
prob <- rep(0,length(x)); b <- rnorm(K,0,sigma0)
#
# Fig 9.6
#
pdf("logistint.pdf")
plot(x,prob,type="n",ylim=c(0,1),ylab="Probability")
for (i in 1:K){
  prob <- exp(beta0+b[i]+beta1*x)/(1+exp(beta0+b[i]+beta1*x))
  lines(x,prob,lty=3)
}
const <- 16*sqrt(3)/(15*pi); const2 <- sqrt(const^2*sigma0^2+1)
margprob <- exp((beta0+beta1*x)/const2)/(1+exp((beta0+beta1*x)/const2))
lines(x,margprob)
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Contraception Data example in Section 9.13
#
library(lme4)
alldat <- read.table("contracept.txt",header=F,na.strings=".",
          col.names=c("id","dose","time","y"))
time2 <- alldat$time*alldat$time
alldat$time2 <- time2
#
# Fit model using Laplace to integrate out the random effects
#
mod1 <- lmer(y~time+time2+dose:time+dose:time2+(1|id),data=alldat,
        family="binomial")
alldat0 <- subset(alldat,dose==0)
alldat1 <- subset(alldat,dose==1)
library(glmmML)
#
# Fit model using Gauss-Hermite to integrate out the random effects
#
modGH <- glmmML(y~time+time2+dose:time+dose:time2,cluster=id,data=alldat,
         family="binomial",method="ghq",n.points=50) # adaptive GH
#
# Now fit the model to each of the dose groups separately to assess the 
# common random effects variance assumption
#
modGH0 <- glmmML(y~time+time2,cluster=id,data=alldat0,family="binomial",
          method="ghq",n.points=50,start.sigma=2.25) # adaptive GH
modGH1 <- glmmML(y~time+time2,cluster=id,data=alldat1,family="binomial",
          method="ghq",n.points=50,start.sigma=2.25) # adaptive GH
#
# For Bayes analysis via MCMC see the WinBUGS file ContraceptWinBUGS.txt
#
# Marginal var-cov structure and also marginal probabilities
#
beta0 <- -3.80568
beta1 <- 1.13320
beta2 <- -0.04192
beta3 <- 0.56444
beta4 <- -0.10955
#
nsim <- 100000
bsamp <- rnorm(nsim,mean=0,sd=sigmaest)
mcov0 <- mcov1 <- matrix(0,nrow=4,ncol=4)
m0 <- m1 <- mvar0 <- mvar1 <- rep(0,4)
expitj0 <- expitprod0 <- expitj1 <- expitprod1 <- NULL
for (j in 1:4){
    expitj0 <- exp(beta0+beta1*j+beta2*j^2+bsamp)/
    	  (1+exp(beta0+beta1*j+beta2*j^2+bsamp))
    expitj1 <- exp(beta0+(beta1+beta3)*j+(beta2+beta4)*j^2+bsamp)/
    	  (1+exp(beta0+(beta1+beta3)*j+(beta2+beta4)*j^2+bsamp))
    m0[j] <- mean(expitj0)
    m1[j] <- mean(expitj1)
    mvar0[j] <- m0[j]*(1-m0[j])
    mvar1[j] <- m1[j]*(1-m1[j])
    mcov0[j,j] <- mvar0[j]
    mcov1[j,j] <- mvar1[j]
# cat(j,m0[j],mvar0[j],mean(expitj0)-m0[j]^2,"\n")
    for (k in 1:j-1){
          expitprod0 <- (exp(beta0+beta1*k+beta2*k^2+bsamp)/
    	        (1+exp(beta0+beta1*k+beta2*k^2+bsamp)))*
                         (exp(beta0+beta1*j+beta2*j^2+bsamp)/
                (1+exp(beta0+beta1*j+beta2*j^2+bsamp)))
          expitprod1 <- (exp(beta0+(beta1+beta3)*k+(beta2+beta4)*k^2+bsamp)/
    	        (1+exp(beta0+(beta1+beta3)*k+(beta2+beta4)*k^2+bsamp)))*
                         (exp(beta0+(beta1+beta3)*j+(beta2+beta4)*j^2+bsamp)/
                (1+exp(beta0+(beta1+beta3)*j+(beta2+beta4)*j^2+bsamp)))
          mcov0[j,k] <- mean(expitprod0) - m0[j]*m0[k]
          mcov0[k,j] <- mcov0[j,k]/sqrt(mvar0[j]*mvar0[k])
	  mcov1[j,k] <- mean(expitprod1) - m1[j]*m1[k]
          mcov1[k,j] <- mcov1[j,k]/sqrt(mvar1[j]*mvar1[k])
# cat(j,k,m0[j],m0[k],mcov0[j,k],"\n")
   }
}
# 
# Plot the marginal curves with data points added.
#
# Fig 9.7
#
pdf("contracept-prob-fitted.pdf",h=4,w=4)
plot(seq(1,4),p01,pch=19,ylim=c(min(p01,p11,m0,m1),max(p01,p11,m0,m1)),
   ylab="Probability of Amenorrhea",xlab="Measurement Occasion",xaxp=c(1,4,3))
points(seq(1,4),p11,pch=21)
points(seq(1,4),m0,pch=17)
points(seq(1,4),m1,pch=2)
legend("topleft",legend=c("High Dose Obs","Low Dose Obs","High Dose Fit",
       "Low Dose Fit"),bty="n",pch=c(21,19,2,17))
dev.off()
#
#
# Table 9.9
#
library(xtable)
print(xtable(mcov0))
print(xtable(mcov1))



#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Contraception data example in Section 9.14, analyzed using GEE
#
library(geepack)
time2 <- alldat$time*alldat$time
alldat$time2 <- time2
modge <- geese(y~time+time2+dose:time+dose:time2,id=id,corstr="exch",
         data=alldat,family="binomial")
modgi <- geese(y~time+time2+dose:time+dose:time2,id=id,corstr="independence",
         data=alldat,family="binomial")
modgu <- geese(y~time+time2+dose:time+dose:time2,id=id,corstr="unstructured",
         data=alldat,family="binomial")
#
# ALR analysis performed using R-forge link from here:
# http://www.biostat.harvard.edu/~carey/
# 
install.packages("alr", repos= c("http://R-Forge.R-project.org", 
                 getOption("repos")))
library(alr)
modalr <- alr(y~time+time2+dose:time+dose:time2,id=id,data=alldat,ainit=.01,
         depm="exchangeable")
#
# Fitted logits versus empirical
#
gamma0 <- modge$beta[1]
gamma1 <- modge$beta[2]
gamma2 <- modge$beta[3]
gamma3 <- modge$beta[4]
gamma4 <- modge$beta[5]
t4 <- seq(1,4,1)
exp(gamma0+gamma1*t4+gamma2*t4^2)/(1+exp(gamma0+gamma1*t4+gamma2*t4^2))
exp(gamma0+(gamma1+gamma3)*t4+(gamma2+gamma4)*t4^2)/
    (1+exp(gamma0+(gamma1+gamma3)*t4+(gamma2+gamma4)*t4^2))
tseq <- seq(1,4,.01)
#
# Fig 9.8
#
pdf("FittedLogitAmer.pdf",h=4,w=4)
plot(gamma0+gamma1*tseq+gamma2*tseq^2~tseq,type="l",xlab="Measurement Occasion",
      ylab="Logit of Amerorrhea",ylim=c(-1.7,0.3),xaxp=c(1,4,3))
lines(gamma0+(gamma1+gamma3)*tseq+(gamma2+gamma4)*tseq^2~tseq,type="l",
      col="grey")
logit0 <- log(p01/(1-p01))
logit1 <- log(p11/(1-p11))
shift <- .03
points(logit0~seq(1+shift,4+shift,1),pch=19)
points(logit1~seq(1-shift,4-shift,1),pch=21)
for(j in 1:4){
      lines(x=c(j+shift,j+shift),y=c(logit0[j]-1.96*sqrt(vlogit0[j]),
            logit0[j]+1.96*sqrt(vlogit0[j])),lwd=2)
      lines(x=c(j-shift,j-shift),y=c(logit1[j]-1.96*sqrt(vlogit1[j]),
            logit1[j]+1.96*sqrt(vlogit1[j])),lwd=2,col="grey")
}
legend("topleft",legend=c("High Dose","Low Dose"),bty="n",pch=c(21,19))
dev.off()
#
# Compare the marginal and conditional estimates:
#
const<- 16*sqrt(3)/(15*pi); 
sigmaest <- 2.25
# Attenuation factor
atten <- sqrt(const^2*sigmaest^2+1)
# These are the GLMM estimates from Table 9.8
beta0 <- -3.80568
beta1 <- 1.13320
beta2 <- -0.04192
beta3 <- 0.56444
beta4 <- -0.10955
# These values can be compared with the GEE estimates in Table 9.12
cat(c(beta0,beta1,beta2,beta3,beta4)/atten,"\n") 
# These are the GLMM standard errors divided by the attenuation factor
cat(c(0.274,.250,0.0524,0.183,0.0576)/atten,"\n")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Section 9.17: Theophylline example with nonlinear LS fitting.
#
TheoNLS <- nlsList(conc~SSfol(Dose, Time, lKe, lKa, lCl),
               data=Theoph) # This is the quick way 
intervals(TheoNLS)
#
# 95% CI on the coefficients constructed "by hand" for plotting
#
j <- 10*(seq(1,12)-1)+1
Dose <- newTheoph$Dose[j] 
sigma <- NULL
est <- lower <- upper <- matrix(0,nrow=12,ncol=3)
resids <- fitted <- time <-  matrix(0,nrow=12,ncol=10)
for (i in 1:12){
    k <- (i-1)*10+seq(1,10)
    t <- newTheoph$Time[k]
    y <- newTheoph$conc[k]
    fit <- nls(y~SSfol(Dose[i],t,lKe,lKa,lCl))
    est[i,] <- coef(fit)
    for (j in 1:3){
    	lower[i,j] <- est[i,j] - 1.96*sqrt(vcov(fit)[j,j])
	upper[i,j] <- est[i,j] + 1.96*sqrt(vcov(fit)[j,j])
    }
    time[i,] <- t
    fitted[i,] <- fitted(fit)
    resids[i,] <- y-fitted(fit)
    sigma[i] <- sqrt(sum(resids[i,]^2)/(10-3))
    resids[i,] <- resids[i,]/sigma[i]
}
#
# Fig 9.9(a)
#
pdf("theoph-indiv1.pdf")
plot(seq(1,12)~est[c(1:12),1],xlim=c(min(lower[,1],upper[,1]),
     max(lower[,1],upper[,1])),xlab=expression(paste("log ",k[e])),
     ylab="Individual Number",yaxp=c(1,12,11))
for(i in 1:12){lines(y=c(i,i),x=c(lower[i,1],upper[i,1]))}
dev.off()
#
# Fig 9.9(b)
#
pdf("theoph-indiv2.pdf")
plot(seq(1,12)~est[c(1:12),2],xlim=c(min(lower[,2],upper[,2]),
     max(lower[,2],upper[,2])),xlab=expression(paste("log ",k[a])),
     ylab="Individual Number",yaxp=c(1,12,11))
for(i in 1:12){lines(y=c(i,i),x=c(lower[i,2],upper[i,2]))}
dev.off()
#
# Fig 9.9(c)
#
pdf("theoph-indiv3.pdf")
plot(seq(1,12)~est[c(1:12),3],xlim=c(min(lower[,3],upper[,3]),
     max(lower[,3],upper[,3])),xlab=expression(paste("log Cl")),
     ylab="Individual Number",yaxp=c(1,12,11))
for(i in 1:12){lines(y=c(i,i),x=c(lower[i,3],upper[i,3]))}
dev.off()
# 
# Fig 9.10
#
pdf("Fitted-nls-Theoph.pdf")
plot(augPred(TheoSTS.nls))
dev.off()
#
# NLMM fitting
#
TheoNLME <- nlme(TheoNLS,
        fixed = lKe + lKa + lCl ~ 1,
        random = lKe + lKa + lCl ~ 1,
        data=Theoph)
summary(TheoNLME) # To get summaries in Likelihood column of Table 9.13
Theo-nlme <- nlme(TheoNLS,fixed = lKe + lKa + lCl ~ 1,
        random = lKe + lKa + lCl ~ 1,data=Theoph)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# For Bayesian analysis see WinBUGS file
#
#
# Create flip flop figure
#
chain1 <- read.table("betachain1.txt")    
chain2 <- read.table("betachain2.txt")    	
nits <- dim(chain1)[1]/3
t1seq <- seq(1,nits,20)
t2seq <- t1seq + nits
t3seq <- t1seq + 2*nits
upper <- 2000
chain1par1 <- chain1[t1seq[1:upper],2]
chain2par1 <- chain2[t1seq[1:upper],2]
par(mfrow=c(1,1))
#
# Fig 9.11(a)
#
pdf("FlipFlopk1.pdf")
plot(chain1par1,typ="l",ylim=c(min(chain1par1,chain2par1),
     max(chain1par1,chain2par1)),xlab="Thinned Iteration Number",
     ylab=expression(beta[1]))
lines(chain2par1,lty=3)
dev.off()
#
chain1par2 <- chain1[t2seq[1:upper],2]
chain2par2 <- chain2[t2seq[1:upper],2]
#
# Fig 9.11(b)
#
pdf("FlipFlopk2.pdf")
plot(chain1par2,typ="l",ylim=c(min(chain1par2,chain2par2),
     max(chain1par2,chain2par2)),xlab="Thinned Iteration Number",
     ylab=expression(beta[2]))
lines(chain2par2,lty=3)
dev.off()
#
chain1par3 <- chain1[t3seq[1:upper],2]
chain2par3 <- chain2[t3seq[1:upper],2]
#
# Fig 9.11(c)
#
pdf("FlipFlopCl.pdf")
plot(chain1par3,typ="l",ylim=c(min(chain1par3,chain2par3),
     max(chain1par3,chain2par3)),xlab="Thinned Iteration Number",
     ylab=expression(beta[3]))
lines(chain2par3,lty=3)
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GEE in Section 9.19
#
#  
# GEE for nonlinear model
#
TheoGEENLS <- nls(conc~SSfol(Dose, Time, lKe, lKa, lCl),data=Theoph)
#
# These are the GEE point estimates in Table 9.14
#
gammaGEE <- coef(TheoGEENLS)
#
#--------------------------------------------------------------
# First form the n_i x 3 matrix D_i using the following functions
#
gamma1deriv <- function(dose,time,gamma){
   u <- dose*exp(gamma[1]+gamma[2])*( exp(-exp(gamma[2])*time) - 
        exp(-exp(gamma[1])*time) )
   v <- exp(gamma[3])*(exp(gamma[2])-exp(gamma[1]))
   deltau <- u + dose*exp(gamma[1]+gamma[2])*exp(-exp(gamma[1])*time)*
             exp(gamma[1])*time
   deltav <- -exp(gamma[1]+gamma[3])
   (u*deltav - v*deltau)/v^2
}
gamma2deriv <- function(dose,time,gamma){
   u <- dose*exp(gamma[1]+gamma[2])*( exp(-exp(gamma[2])*time) - 
        exp(-exp(gamma[1])*time) )
   v <- exp(gamma[3])*(exp(gamma[2])-exp(gamma[1]))
   deltau <- u - dose*exp(gamma[1]+gamma[2])*exp(-exp(gamma[2])*time)*
             exp(gamma[2])*time
   deltav <- exp(gamma[2]+gamma[3])
   (u*deltav - v*deltau)/v^2
}
gamma3deriv <- function(dose,time,gamma){
   u <- dose*exp(gamma[1]+gamma[2])*( exp(-exp(gamma[2])*time) - 
        exp(-exp(gamma[1])*time) )
   v <- exp(gamma[3])*(exp(gamma[2])-exp(gamma[1]))
   mu <- u/v
   -mu
}
Darray <- array(0,dim=c(12,11,3))
j <- 11*(seq(1,12)-1)+1
Dose <- Theoph$Dose[j] 
counter <- 0
for (i in 1:12){
    for (j in 1:11){
    	counter <- counter + 1
        Darray[i,j,1] <- gamma1deriv(Dose[i],Theoph$Time[counter],gammaGEE)
	Darray[i,j,2] <- gamma2deriv(Dose[i],Theoph$Time[counter],gammaGEE)
	Darray[i,j,3] <- gamma3deriv(Dose[i],Theoph$Time[counter],gammaGEE)
    }
}
# Need to form 3 matrices for the sandwich var-cav matrix
#
Breadmat <- Fillingmat <- matrix(0,3,3)
residi <- matrix(0,11,1)
covYi <- matrix(0,11,11)
for (i in 1:12){
    Breadmat <- Breadmat + t(Darray[i,,]) %*% Darray[i,,]	
    k <- (i-1)*11 + seq(1,11)
    residi <- resid(TheoGEENLS)[k]
    covYi <- residi %*% t(residi)
    Fillingmat <- Fillingmat + t(Darray[i,,]) %*% covYi %*% Darray[i,,]
}
vargammaGEE <- solve(Breadmat) %*% Fillingmat %*% solve(Breadmat)
#
# Table 9.14 standard errors
#
print(sqrt(diag(vargammaGEE)))
#
# Theophylline example in Section 9.20
#
#
# Fig 9.12(a)
#
par(mfrow=c(1,1))
pdf("Theoph-qqnormke.pdf",h=4,w=4)
qqnorm(est[,1],main="",ylab=expression(paste("Sample quantiles log ",k[e])))
dev.off()
#
# Fig 9.12(c)
#
pdf("Theoph-qqnormka.pdf",h=4,w=4)
qqnorm(est[,2],main="",ylab=expression(paste("Sample quantiles log ",k[a])))
dev.off()
#
# Fig 9.12(c)
#
pdf("Theoph-qqnormCl.pdf",h=4,w=4)
qqnorm(est[,3],main="",ylab=expression(paste("Sample quantiles log ",Cl)))
dev.off()
#
# Fig 9.12(b)
#
pdf("Theoph-keka.pdf",h=4,w=4)
plot(est[,1],est[,2],ylab=expression(paste("log ",k[a])),
     xlab=expression(paste("log ",k[e])))
abline(0,1) # To check if flop-flop likely, ie are ka and ke of similar size?
dev.off()
#
# Fig 9.12(d)
#
pdf("Theoph-keCl.pdf",h=4,w=4)
plot(est[,1],est[,3],ylab=expression(paste("log ",Cl)),
     xlab=expression(paste("log ",k[e])))
dev.off()
#
# Fig 9.12(f)
#
pdf("Theoph-kaCl.pdf",h=4,w=4)
plot(est[,2],est[,3],ylab=expression(paste("log ",Cl)),
     xlab=expression(paste("log ",k[a])))
dev.off()
#
# Now examine residuals
#
# Fig 9.13(a)
#
pdf("theoph-nls-resids1.pdf",h=5,w=5)
plot(resids~fitted,ylab="Residuals",xlab="Fitted")
abline(0,0,lty=2)
dev.off()
#
# Fig 9.13(b)
#
pdf("theoph-nls-resids2.pdf",h=5,w=5)
plot(resids~time,ylab="Residuals",xlab="Time (hours)")
abline(0,0,lty=2)
dev.off()
#
# Read in lognormal model output
#
temp <- read.table("TheophWinBUGSLognormalOut.txt",nrows=330)
bayesresids <- temp[171:290,6]/temp[294,6]
bayesfitted <-  temp[26:145,6]
#
# Fig 9.13(c)
#
pdf("BayesResids1.pdf",h=5,w=5)
plot(bayesresids~bayesfitted,ylab="Residuals",xlab="Fitted")
abline(0,0,lty=2)
dev.off()
#
# Fig 9.13(d)
#
pdf("BayesResids2.pdf",h=5,w=5)
plot(bayesresids~as.vector(time),ylab="Residuals",xlab="Time (hours)")
abline(0,0,lty=2)
dev.off()
#
# Read in power model output
#
tempPow <- read.table("TheophWinBUGSPowerOut.txt",nrows=297)
bayesfittedPow <-  tempPow[26:145,6] # no fitted at zero
bayesresidsPow <- newTheoph$conc-bayesfittedPow
sigma0 <- tempPow[296,6]
sigma1 <- tempPow[297,6]
gamma <- tempPow[146,6]
#
# estimate of sigma under the power model
#
sigmaPow <- sqrt(sigma0*sigma0+sigma1*sigma1*bayesfittedPow^gamma)
#
# Form residuals
#
bayesresidsPow <- bayesresidsPow/sigmaPow
#
# Fig 9.13(e)
#
pdf("BayesResidsPow1.pdf",h=5,w=5)
par(mfrow=c(1,1))
plot(bayesresidsPow~bayesfittedPow,ylab="Residuals",xlab="Fitted")
abline(0,0,lty=2)
dev.off()
#
# Fig 9.13(f)
#
pdf("BayesResidsPow2.pdf",h=5,w=5)
plot(bayesresidsPow~as.vector(time),ylab="Residuals",xlab="Time (hours)")
abline(0,0,lty=2)
dev.off()
#
# Posteriors for power model
#
par(mfrow=c(1,1))
gammat <- read.table("WinBUGSgamma.txt")
gamma <- gammat[,2]
#
# Fig 9.14(e)
#
pdf("Powgamma.pdf",h=4,w=4)
hist(gamma,main="",xlab=expression(gamma))
dev.off()
#
# Fig 9.14(a)
#
sigma0t <- read.table("WinBUGSsigma0.txt")
sigma0 <- sigma0t[,2]
pdf("Powsigma0.pdf",h=4,w=4)
hist(sigma0,main="",xlab=expression(sigma[0]))
dev.off()
sigma1t <- read.table("WinBUGSsigma1.txt")
sigma1 <- sigma1t[,2]
#
# Fig 9.14(c)
#
pdf("Powsigma1.pdf",h=4,w=4)
hist(sigma1,main="",xlab=expression(sigma[1]))
dev.off()
#
# Fig 9.14(b)
#
pdf("Powsigma1sigma0.pdf",h=4,w=4)
plot(sigma1~sigma0,xlab=expression(sigma[0]),ylab=expression(sigma[1]))
dev.off()
#
# Fig 9.14(d)
#
pdf("Powgammasigma0.pdf",h=4,w=4)
plot(gamma~sigma0,ylab=expression(gamma),xlab=expression(sigma[0]))
dev.off()
#
# Fig 9.14(f)
#
pdf("Powgammasigma1.pdf",h=4,w=4)
plot(gamma~sigma1,ylab=expression(gamma),xlab=expression(sigma[1]))
dev.off()
#
#  Now posteriors on pars of interest
#
Clmedt <- read.table("WinBUGSClmed.txt")
Clmed <- Clmedt[,2]
#
# Fig 9.15(e)
#
pdf("PowClmed.pdf",h=4,w=4)
hist(Clmed,main="",xlab=expression(Cl))
dev.off()
kamedt <- read.table("WinBUGSkamed.txt")
kamed <- kamedt[,2]
#
# Fig 9.15(a)
#
pdf("Powkamed.pdf",h=4,w=4)
hist(kamed,main="",xlab=expression(k[a]),nclass=20)
dev.off()
kemedt <- read.table("WinBUGSkemed.txt")
kemed <- kemedt[,2]
#
# Fig 9.15(c)
#
pdf("Powkemed.pdf",h=4,w=4)
hist(kemed,main="",xlab=expression(k[e]))
dev.off()
#
# Fig 9.15(b)
#
pdf("Powkamedkemed.pdf",h=4,w=4)
plot(kamed~kemed,xlab=expression(k[e]),ylab=expression(k[a]))
dev.off()
#
# Fig 9.15(d)
#
pdf("PowClmedkemed.pdf",h=4,w=4)
plot(Clmed~kemed,ylab=expression(Cl),xlab=expression(k[e]))
dev.off()
#
# Fig 9.15(f)
#
pdf("PowClmedkamed.pdf",h=4,w=4)
plot(Clmed~kamed,ylab=expression(Cl),xlab=expression(k[a]))
dev.off()
