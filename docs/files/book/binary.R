#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 7
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Head injury data EDA
#
par(mfrow=c(1,1))
age <- seq(1,3)
p1 <- c(9/56,19/34,7/8)
p2 <- c(58/87,45/56,20/20)
p3 <- c(5/16,21/39,19/21)
p4 <- c(32/45,61/72,42/49)
agevar <- seq(1,3)
#
# Fig 7.1(a)
#
pdf("headinjury-EDA1.pdf",h=5,w=5)
plot(p1~agevar,type="n",ylim=c(0,1),xlab="Age",ylab="Probability of Death",ax=F)
axis(2)
box() 
axis(1,at=c(1,2,3))
lines(agevar,p1,lty=1)
lines(agevar,p2,lty=2)
lines(agevar,p3,lty=3)
lines(agevar,p4,lty=4)
legend("bottomright",lty=c(1,2,3,4),legend=c("Haem=No,Pup=Good","Haem=No,Pup=Poor",
       "Haem=Yes,Pup=Good","Haem=Yes,Pup=Poor"),bty="n")
dev.off()
#
p5 <- c(5/82,6/50,12/18)
p6 <- c(11/35,7/23,7/9)
p7 <- c(7/31,14/52,25/40)
p8 <- c(12/28,15/36,17/24)
agevar <- seq(1,3)
#
# Fig 7.1(b)
#
pdf("headinjury-EDA2.pdf",h=5,w=5)
plot(p1~agevar,type="n",ylim=c(0,1),xlab="Age",ylab="Probability of Death",ax=F)
axis(2)
box()
axis(1,at=c(1,2,3))
lines(agevar,p5,lty=1)
lines(agevar,p6,lty=2)
#
# BPD and birthweight example
#
bw <- read.table("birthweight.txt",header=T)
attach(bw)
myexpit <- function(x,b0,b1){expit <- exp(b0+b1*x)/( 1+exp(b0+b1*x) )}
mycllit <- function(x,b0,b1){cllit <- 1-exp(-exp(b0+b1*x))}
lrmod1 <- glm(BPD~birthweight,family=binomial)
clmod1 <- glm(BPD~birthweight,family=binomial(link="cloglog"))
par(mfrow=c(1,1))
pdf("bwfig1.pdf",height=4.5,width=4.5)
plot(birthweight,BPD,pch="|",xlab="Birthweight (grams)")
x <- seq(min(birthweight),max(birthweight))
lines(x,myexpit(x,b0=lrmod1$coeff[1],b1=lrmod1$coeff[2]),lty=2)
lines(x,mycllit(x,b0=clmod1$coeff[1],b1=clmod1$coeff[2]),lty=3)
legend("topright",legend=c("Logistic","Comp Log Log"),lty=2:3,bty="n")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outcome after head injury example in Section 7.6
#
yagg <- c(9,5,5,7,58,11,32,12,19,6,21,14,45,7,61,15,7,12,19,25,20,7,42,17)
zagg <- c(47,77,11,24,29,24,13,16,15,44,18,38,11,16,11,21,1,6,2,15,0,2,7,7)
haem <- as.factor(c(0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1))
pup <- as.factor(c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1))
coma <- as.factor(c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1))
agec <- as.factor(c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2))
head <- data.frame(yagg,zagg,haem,pup,coma,agec)
#
# Forward selection and backwards elimination
#
modnull <- glm(cbind(yagg,zagg)~1,family="binomial")
modfull <- glm(cbind(yagg,zagg)~haem*pup*coma*agec,family="binomial",data=head)
# Backwards and forwards
backstep <- step(modfull,scope="cbind(yagg,zagg)~haem*pup*coma*agec")
summary(backstep) # To find deviance of chosen model
forstep <- step(modnull,scope="cbind(yagg,zagg)~haem*pup*coma*agec")
summary(forstep) # To find deviance of chosen model
#
# Now carry out exhaustivs search over all models that obey the hierarchy principle
#
# Create matrix of 0's and 1's to keep track of interactions and model elements
#
design.it <- function(k){
  design <- matrix(nrow=2^k,ncol=k)
  for(i in 1:k){
    design[,i]=rep(rep(c(0,1),each=2^(k-i)),2^(i-1))
  }
  return(design)
}
#
# is a subset of b? 
#
is.subset <- function(a,b,k){
  yes <- 1
  for(i in 1:k)
    {
      yes=yes*(a[i]<=b[i])
    }
  return(yes)
}
#
# is this a hierarchical  model? 
#
is.hier <- function(int.mat,model,p,k){
  yes=1
  for(i in 1:p){
    if(model[i]){
      for(j in 1:p){
        if(is.subset(int.mat[j,],int.mat[i,],k)){
          yes= yes * model[j]
        }
      }
    }
  }
  return(yes)
}
#
# matrix to help identify all possible interactions
#
interactions<-design.it(4)[-1,]

##matrix to identify which interactions to include####
possible.model<-design.it(15)

check.hier<-function(c){
  return(is.hier(interactions,c,15,4))
}
#
# keep models that obey hierarchy principle
#
check.it <- apply(possible.model,1,check.hier)
keep <- which(check.it==1)
hier.models <- possible.model[keep,]
#
vars <- names(head)[-c(1,2)]
colnames(interactions)=vars
varints <- apply(interactions, 1, function(x){
       paste(colnames(interactions)[as.logical(x)], collapse="*")})
#
library("MASS")
get.aic<-function(model,var,hdata,k){
  my.formula <- as.formula(paste(c("cbind(yagg,zagg) ~ 1", var[as.logical(model)]), 
                collapse=" + "))
  my.model<-glm(my.formula, family = binomial, data=hdata)
#
  if (k==1) return(extractAIC(my.model)[2]) 
  if (k==2) return(extractAIC(my.model,k=log(sum(yagg+zagg)))[2])
}

# calculate aic for each hierarchical model
get.my.aic<-function(c){return(get.aic(c,varints,head,k=1))}
get.my.bic<-function(c){return(get.aic(c,varints,head,k=2))}
#
aic.vals <- apply(hier.models,1,get.my.aic)
bic.vals <- apply(hier.models,1,get.my.bic)
#
# which is the best?
#
min(aic.vals)
hier.models[which.min(aic.vals),]
best.model <- varints[as.logical(hier.models[which.min(aic.vals),])]
best.model
#
# Which is the best model uing BIC as crition?
#
min(bic.vals)
hier.models[which.min(bic.vals),]
best.model<-varints[as.logical(hier.models[which.min(bic.vals),])]
best.model
#
modreport <- glm(cbind(yagg,zagg)~haem+pup+coma+agec+haem*pup+haem*agec+pup*agec,
             family="binomial")
# MLEs and ses in Table 7.2
#
summary(modreport)
#
# Now for Bayesian analysis with INLA
#
headdat <- list(yagg=yagg,zagg=zagg,haem=haem,pup=pup,coma=coma,agec=agec)
source("http://www.math.ntnu.no/inla/givemeINLA.R")
inla.upgrade()
library(INLA)
upper <- 10000
W <- (log(upper)/1.96)^2
inlamod <- inla(yagg~haem+pup+coma+agec+haem*pup+haem*agec+pup*agec,family="binomial",
           Ntrials=yagg+zagg,data=headdat,control.fixed=list(mean.intercept=c(0),
           prec.intercept=c(0),mean=c(0,0,0,0,0,0,0),prec=c(seq(1/W,7))))
summary(inlamod)
#
#
fitreport <- fitted(modreport)
fitodds <- fitreport/(1-fitreport)
obsodds <- (yagg+.5)/(zagg+.5)
denoms <- yagg+zagg
p1 <- c(9/56,19/34,7/8)
seq1 <- c(1,9,17)
fit1 <- fitreport[seq1]
p2 <- c(58/87,45/56,20/20)
seq2 <- c(5,13,21)
fit2 <- fitreport[seq2]
p3 <- c(5/16,21/39,19/21)
seq3 <- c(3,11,19)
fit3 <- fitreport[seq3]
p4 <- c(32/45,61/72,42/49)
seq4 <- c(7,15,23)
fit4 <- fitreport[seq4]
agevar <- seq(1,3)
#
# Fig 7.3(a)
#
pdf("headinjury-FIT1.pdf",h=5,w=5)
plot(p1~agevar,type="n",ylim=c(0,1),xlab="Age",ylab="Probability of Death",ax=F)
axis(2)
box() 
axis(1,at=c(1,2,3))
lines(agevar,p1,lty=1)
lines(agevar,p2,lty=2)
lines(agevar,p3,lty=3)
lines(agevar,p4,lty=4)
legend("bottomright",lty=c(1,2,3,4),legend=c("Haem=No,Pup=Good","Haem=No,Pup=Poor",
       "Haem=Yes,Pup=Good","Haem=Yes,Pup=Poor"),bty="n")
points(agevar,fit1,pch=21)
points(agevar,fit2,pch=21)
points(agevar,fit3,pch=21)
points(agevar,fit4,pch=21)
lines(x=c(1,1),c(p1[1],fit1[1]),lty=1)
lines(x=c(1,1),c(p2[1],fit2[1]),lty=2)
lines(x=c(1,1),c(p3[1],fit3[1]),lty=3)
lines(x=c(1,1),c(p4[1],fit4[1]),lty=4)
lines(x=c(2,2),c(p1[2],fit1[2]),lty=1)
lines(x=c(2,2),c(p2[2],fit2[2]),lty=2)
lines(x=c(2,2),c(p3[2],fit3[2]),lty=3)
lines(x=c(2,2),c(p4[2],fit4[2]),lty=4)
lines(x=c(3,3),c(p1[3],fit1[3]),lty=1)
lines(x=c(3,3),c(p2[3],fit2[3]),lty=2)
lines(x=c(3,3),c(p3[3],fit3[3]),lty=3)
lines(x=c(3,3),c(p4[3],fit4[3]),lty=4)
dev.off()
#
#
p5 <- c(5/82,6/50,12/18)
seq5 <- c(2,10,18)
fit5 <- fitreport[seq5]
p6 <- c(11/35,7/23,7/9)
seq6 <- c(6,14,22)
fit6 <- fitreport[seq6]
p7 <- c(7/31,14/52,25/40)
seq7 <- c(4,12,20)
fit7 <- fitreport[seq7]
p8 <- c(12/28,15/36,17/24)
seq8 <- c(8,16,24)
fit8 <- fitreport[seq8]
agevar <- seq(1,3)
#
# Fig 7.3(b)
#
pdf("headinjury-FIT2.pdf",h=5,w=5)
plot(p1~agevar,type="n",ylim=c(0,1),xlab="Age",ylab="Probability of Death",ax=F)
axis(2)
lines(agevar,p5,lty=1)
lines(agevar,p6,lty=2)
lines(agevar,p7,lty=3)
lines(agevar,p8,lty=4)
box()
axis(1,at=c(1,2,3))
points(agevar,fit5,pch=21)
points(agevar,fit6,pch=21)
points(agevar,fit7,pch=21)
points(agevar,fit8,pch=21)
lines(x=c(1,1),c(p5[1],fit5[1]),lty=1)
lines(x=c(1,1),c(p6[1],fit6[1]),lty=2)
lines(x=c(1,1),c(p7[1],fit7[1]),lty=3)
lines(x=c(1,1),c(p8[1],fit8[1]),lty=4)
lines(x=c(2,2),c(p5[2],fit5[2]),lty=1)
lines(x=c(2,2),c(p6[2],fit6[2]),lty=2)
lines(x=c(2,2),c(p7[2],fit7[2]),lty=3)
lines(x=c(2,2),c(p8[2],fit8[2]),lty=4)
lines(x=c(3,3),c(p5[3],fit5[3]),lty=1)
lines(x=c(3,3),c(p6[3],fit6[3]),lty=2)
lines(x=c(3,3),c(p7[3],fit7[3]),lty=3)
lines(x=c(3,3),c(p8[3],fit8[3]),lty=4)
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Aircract fastener example in Section 7.6
#
x <- c(2500,2700,2900,3100,3300,3500,3700,3900,4100,4300)
n <- c(50,70,100,60,40,85,90,50,80,65)
y <- c(10,17,30,21,18,43,54,33,60,51)
#
# Fig 7.4
#
pdf("aircraft1.pdf")
par(mfrow=c(1,1))
plot(y/n~x,ylim=c(0,1),ylab="Probability of Failure",xlab="Pressure Load (psi)",
     xlim=c(2000,4800))
logistmod <- glm(cbind(y,n-y)~x,family="binomial")
cloglogmod <- glm(cbind(y,n-y)~x,family=binomial(link="cloglog"))
loglogmod <- glm(cbind(n-y,y)~x,family=binomial(link="cloglog"))
xseq <- seq(2000,6000,1)
fitted <- exp(logistmod$coef[1]+xseq*logistmod$coef[2])/(1+exp(logistmod$coef[1]+
              xseq*logistmod$coef[2]))
fittedcloglog <- 1-exp(-exp(cloglogmod$coef[1]+xseq*cloglogmod$coef[2]))
fittedloglog <- exp(-exp(loglogmod$coef[1]+xseq*loglogmod$coef[2]))
lines(fitted~xseq)
lines(fittedcloglog~xseq,lty=2)
lines(fittedloglog~xseq,lty=3)
abline(v=3000,lty=2)
legend("topleft",lty=c(1:3),bty="n",legend=c("Logistic","Comp log-log","Log-log"))
dev.off()
#
# Numerical summaries 
#
cat("exp(beta1) 95 CI: ",exp(logistmod$coef[2]-1.96*sqrt(vcov(logistmod)[2,2])),"\n")
cat("exp(beta1) 95 CI: ",exp(logistmod$coef[2]+1.96*sqrt(vcov(logistmod)[2,2])),"\n")
cat("exp(500 beta1) 95 CI: ",exp(500*logistmod$coef[2]-1.96*sqrt(500*vcov(logistmod)[2,2])),"\n")
cat("exp(500 beta1) 95 CI: ",exp(500*logistmod$coef[2]+1.96*sqrt(500*vcov(logistmod)[2,2])),"\n")
#
# Now for prediction at 3000 psi
#
xtilde <- 3000
theta <- logistmod$coef[1]+xtilde*logistmod$coef[2]
vtheta <- vcov(logistmod)[1,1] + 2*xtilde*vcov(logistmod)[1,2] + 
               xtilde^2*vcov(logistmod)[2,2]
loci <- exp(theta-1.96*sqrt(vtheta))/(1+exp(theta-1.96*sqrt(vtheta)))
hici <- exp(theta+1.96*sqrt(vtheta))/(1+exp(theta+1.96*sqrt(vtheta)))
cat("95 CI for prob at 3000: ",exp(logistmod$coef[2]-1.96*sqrt(vcov(logistmod)[2,2])),"\n")
#
# Residual analysis of Section 7.8
#
logistres <- sqrt(n)*(y/n - fitted(logistmod))/sqrt(fitted(logistmod)*(1-fitted(logistmod)))
clogres <- sqrt(n)*(y/n - fitted(cloglogmod))/sqrt(fitted(cloglogmod)*(1-fitted(cloglogmod)))
loglogres <- sqrt(n)*(y/n - (1-fitted(loglogmod)))/sqrt(1-fitted(loglogmod)*(fitted(loglogmod)))
rlo <- min(logistres,clogres,loglogres)
rhi <- max(logistres,clogres,loglogres)
#
# Fig 7.6(a)
#
par(mfrow=c(1,1))
pdf("air-res1.pdf",height=4.5,width=4.5)
plot(logistres~x,ylab="Pearson Residuals",xlab="Pressure Load (psi)",ylim=c(rlo,rhi))
abline(0,0,lty=2)
dev.off()
#
# Fig 7.6(b)
#
par(mfrow=c(1,1))
pdf("air-res2.pdf",height=4.5,width=4.5)
plot(clogres~x,ylab="Pearson Residuals",xlab="Pressure Load (psi)",ylim=c(rlo,rhi))
abline(0,0,lty=2)
dev.off()
#
# Fig 7.6(c)
#
par(mfrow=c(1,1))
pdf("air-res3.pdf",height=4.5,width=4.5)
plot(loglogres~x,ylab="Pearson Residuals",xlab="Pressure Load (psi)",ylim=c(rlo,rhi))
abline(0,0,lty=2)
dev.off()
#
# Bayesian analysis using MCMC
#
#----------------------------------------------------
# Function to evaluate log-likelihood
#
bloglik <- function(y,n,x,beta){
	nexpts <- length(y)
	bloglik <- 0
	for (i in 1:nexpts){
	    bloglik <- bloglik + y[i]*(beta[1]+beta[2]*x[i]) - 
                       n[i]*log(1+exp(beta[1]+beta[2]*x[i]))
	}
	bloglik
}
library(MASS)
burnin <- 1000
nsims <- burnin + 10000
logistmod <- glm(cbind(y,n-y)~x,family="binomial")
beta <- matrix(0,nrow=nsims,ncol=2)
beta[1,1] <- logistmod$coeff[1]
beta[1,2] <- logistmod$coeff[2]
vmat <- vcov(logistmod)
c <- 2; count <- 0
prevloglik <- bloglik(y,n,x,beta[1,])
for(i in 2:nsims){
        beta[i,] <- mvrnorm(1,mu=c(beta[i-1,]),Sigma=c*c*vmat)
        newloglik <- bloglik(y,n,x,beta[i,])
        u <- runif(1)
        if (log(u) < newloglik-prevloglik){
                prevloglik <- newloglik
                count <- count + 1
        }
        else {
                beta[i,] <- beta[i-1,]
        }
}
cat("Accept = ",count/nsims,"\n")
cat("beta0 2.5 50 97.5: ",quantile(beta[seq(burnin+1,nsims),1],p=c(0.025,.5,.975)),"\n")
cat("beta1 2.5 50 97.5: ",quantile(beta[seq(burnin+1,nsims),2],p=c(0.025,.5,.975)),"\n")
cat("exp(beta1) 2.5 50 97.5: ",exp(quantile(beta[seq(burnin+1,nsims),2],
                p=c(0.025,.5,.975))),"\n")
#
# Examination of trace plots to assess convergence (not in book)
#
par(mfrow=c(2,1))
plot(beta[,1])
plot(beta[,2])
par(mfrow=c(1,1))
#
# Fig 7.5(a)
#
pdf("aircraftpost1.pdf")
hist(beta[seq(burnin+1,nsims),1],xlab=expression(beta[0]),main="")
dev.off()
#
# Fig 7.5(b)
#
pdf("aircraftpost2.pdf")
hist(beta[seq(burnin+1,nsims),2],xlab=expression(beta[1]),main="")
dev.off()
#
# Fig 7.5(c)
#
pdf("aircraftpost3.pdf")
plot(beta[seq(burnin+1,nsims),1],beta[seq(burnin+1,nsims),2],
     xlab=expression(beta[0]),ylab=expression(beta[1]))
dev.off()
#
# Fig 7.5(d)
#
xtilde <- 3000
theta <- beta[,1] + xtilde*beta[,2]
prob <- exp(theta)/(1+exp(theta))
cat("theta prob 95%: ",quantile(prob,p=c(0.025,0.975)),"\n")
pdf("aircraftpost4.pdf")
hist(prob,main="",xlab="Failure Probability at 3000 psi")
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Tumor appearance within mice example in Section 7.7
#
# This function returns the mean, variance and probabilities for an extended
# hypergeometric distribution
#
myhyper <- function(n0dot,n1dot,ndot1,alpha){
	const <- 0
	lo <- max(0,ndot1-n1dot)
	hi <- min(n0dot,ndot1)
	for (u in lo:hi){
	    const <- const + exp(lchoose(n0dot,u) + lchoose(n1dot,ndot1-u) + alpha*u)
	}
	hmean <- hvar <- h2 <- 0
	probs <- NULL
	counter <- 1
	for (u in lo:hi){
	    hmean <- hmean + exp(log(u) + lchoose(n0dot,u) + lchoose(n1dot,ndot1-u) + 
                     alpha*u)
	    h2 <- h2 + exp(2*log(u) + lchoose(n0dot,u) + lchoose(n1dot,ndot1-u) + alpha*u)
	    probs[counter] <- exp(lchoose(n0dot,u) + lchoose(n1dot,ndot1-u) + alpha*u)/const
	    counter <- counter + 1
	}
	hmean <- hmean/const
	hvar <- h2/const - hmean^2
	list(hmean=hmean,hvar=hvar,probs=probs)
}
#
# Tumor data in Table 7.4
#
n01 <- 21
n0dot <- 23
n1dot <- 32
ndot1 <- 40
#
# Unconditional likelihood analysis
#
x <- c(0,1)
logitmod <- glm(cbind(c(2,13),c(21,19))~x,family="binomial")
summary(logitmod)
#
# Conditional likelihood analysis: equate the expected value with the empirical 
# mean, as in equation (7.23) in the book
#
f <- function(alpha,n0dot,n1dot,ndot1,n01){
     myhyper(n0dot,n1dot,ndot1,alpha)$hmean - n01}
MCLE <- uniroot(f,interval <- c(-5,5),n0dot,n1dot,ndot1,n01)
calpha <- MCLE$root
cat("Unconditional estimate = ",MCLE$root,"\n")
con <- myhyper(n0dot,n1dot,ndot1,calpha)
cat("Unconditional se: ",sqrt(con$hvar),"\n")
#
# Fisher's exact probs:
#
fish <- myhyper(n0dot,n1dot,ndot1,0)
#
# the last 3 of these probabilities give the one-sided p-value
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Further Section 7.8 examples on assessment of assumptions
#
# Outcome after head injury
#
modreport <- glm(cbind(yagg,zagg)~haem+pup+coma+agec+haem*pup+haem*agec+pup*agec,
     family="binomial")
fitreport <- fitted(modreport)
headresids <- (yagg-(yagg+zagg)*fitreport)/sqrt((yagg+zagg)*fitreport*(1-fitreport))
#
# Fig 7.7
#
pdf("headqqnorm.pdf",h=5,w=5)
qqnorm(headresids,main="")
dev.off()
#
# BPD and birthweight
#
bwlogit <- glm(cbind(1-BPD,BPD)~birthweight,family="binomial")
bwcloglog <- glm(cbind(1-BPD,BPD)~birthweight,family=binomial(link="cloglog"))
bwlogistres <- (BPD - fitted(bwlogit))/(fitted(bwlogit)*(1-fitted(bwlogit)))
clogres <- (BPD - fitted(bwcloglog))/(fitted(bwcloglog)*(1-fitted(bwcloglog)))
rlo <- min(clogres,bwlogistres)
rhi <- max(clogres,bwlogistres)
#
# Fig 7.8(a)
#
par(mfrow=c(1,1))
pdf("bwfig2.pdf",height=4.5,width=4.5)
plot(bwlogistres~birthweight,ylab="Pearson Residuals",xlab="Birthweight (grams)",
     ylim=c(rlo,rhi))
abline(0,0,lty=2)
dev.off()
#
# Fig 7.8(b)
#
pdf("bwfig3.pdf",height=4.5,width=4.5)
plot(clogres~birthweight,ylab="Pearson Residuals",xlab="Birthweight (grams)",
     ylim=c(rlo,rhi))
abline(0,0,lty=2) 
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# Collapsibility example in Section 7.9
#
y0 <- c(95,10,90,5)
y1 <- c(5,90,10,95)
x <- c(0,0,1,1)
z <- c(0,1,0,1)
modxz <- glm(cbind(y1,y0)~x+z,family="binomial")
summary(modxz)
modx <- glm(cbind(y1,y0)~x,family="binomial")
summary(modx)











