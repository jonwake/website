#
#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 10
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Lidar initial plot
# 
# Fig 10.1
#
library(SemiPar)
data(lidar)
attach(lidar)
pdf("lidarfig1.pdf",height=4.5,width=4.5)
plot(logratio~range,xlab="Range (m)",ylab="Log Ratio")
dev.off()
#
# Ethanol initial plot
#
# Fig 10.2
#
library(SemiPar)
data(ethanol)
attach(ethanol)
library(scatterplot3d)
pdf("ethanol-3d-data.pdf")
scatterplot3d(E,C,NOx,angle=60)
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Ridge and lasso explanation plot
#
#
library(ellipse)
#
# Plot an ellipse corresponding to a 95% probability region for a 
# bivariate normal distribution with mean scale, variances scale and 
# correlation 0.8. 
#
x <- y <- seq(-12,12,.1)
covmat <- matrix(c(1,0.9,0.9,1),nrow=2,ncol=2)
#
# Fig 10.4
#
pdf("ridgelasso1.pdf",height=4,width=4)
par(pty="s")
plot(ellipse(x=covmat,scale=c(1.56,1.56),centre=c(6,7),level=0.95),
type = 'l',xlim=c(-12,12),ylim=c(-12,12),xlab=expression(beta[1]),
ylab=expression(beta[2]),xaxt="n",yaxt="n") 
abline(a=0,b=0,col="grey")
abline(v=0,col="grey")
lines(ellipse(x=covmat,scale=c(1,1),centre=c(6,7),level=0.95))
lines(ellipse(x=covmat,scale=c(.4,.4),centre=c(6,7),level=0.95))
axis(1,at=0)
axis(2,at=0)
covmatc <- matrix(c(1,0,0,1),nrow=2,ncol=2)
lines(ellipse(x=covmatc,scale=c(2,2),level=0.85))
dev.off()
#
# Fig 10.6
#
pdf("ridgelasso2.pdf",height=4,width=4) # Fig 10.6
par(pty="s")
plot(ellipse(x=covmat,scale=c(1.58,1.58),centre=c(6,7),level=0.95),
type = 'l',xlim=c(-12,12),ylim=c(-12,12),xlab=expression(beta[1]),
ylab=expression(beta[2]),xaxt="n",yaxt="n") 
abline(a=0,b=0,col="grey")
abline(v=0,col="grey")
lines(ellipse(x=covmat,scale=c(1,1),centre=c(6,7),level=0.95))
lines(ellipse(x=covmat,scale=c(.4,.4),centre=c(6,7),level=0.95))
axis(1,at=0)
axis(2,at=0)
edge <- 5.4
lines(x=c(-edge,0),y=c(0,edge))
lines(x=c(-edge,0),y=c(0,-edge))
lines(x=c(0,edge),y=c(-edge,0))
lines(x=c(0,edge),y=c(edge,0))
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Hard/soft plot
#
x <- y <- seq(-5,5,.1)
# 
# Fig 10.4(a)
#
pdf("hardsoft1.pdf",height=4.5,width=4.5)  
plot(x,y,type="n",main="",xlab=expression(paste("Least Squares Estimate")),
     ylab="Ridge Regression Estimate",xaxt="n",yaxt="n")
abline(a=0,b=1,lty=2)
lines(x=c(-5,5),y=c(-3,3),lwd=2)
abline(a=0,b=0,col="grey")
abline(v=0,col="grey")
dev.off()
# 
# Fig 10.4(b)
#
pdf("hardsoft2.pdf",height=4.5,width=4.5)  # Fig 10.4(b)
plot(x,y,type="n",main="",xlab=expression(paste("Least Squares Estimate")),
     ylab="Soft Thresholding Estimate",xaxt="n",yaxt="n")
abline(a=0,b=1,lty=2)
lines(x=c(-5,-2),y=c(-3,0),lwd=2)
lines(x=c(5,2),y=c(3,0),lwd=2)
lines(x=c(-2,2),y=c(0,0),lwd=2)
abline(a=0,b=0,col="grey")
abline(v=0,col="grey")
dev.off()
# 
# Fig 10.4(c)
#
pdf("hardsoft3.pdf",height=4.5,width=4.5) # Fig 10.4(c)
x <- y <- seq(-5,5,.1)
plot(x,y,type="n",main="",xlab=expression(paste("Least Squares Estimate")),
     ylab="Hard Thresholding Estimate",xaxt="n",yaxt="n")
abline(a=0,b=1,lty=2)
lines(x=c(-5,-2),y=c(-5,-2),lwd=2)
lines(x=c(5,2),y=c(5,2),lwd=2)
lines(x=c(-2,2),y=c(0,0),lwd=2)
abline(a=0,b=0,col="grey")
abline(v=0,col="grey")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Prostate cancer -- full data shrinkage with ridge regression and the lasso.
#
library(lasso2)
data(Prostate)
lambda <- seq(1,10000,10)
betaridge <- matrix(0,nrow=length(lambda),ncol=8)
df <- NULL
x <- as.matrix(Prostate[,1:8])
y <- Prostate[,9] - mean(Prostate[,9])
means <- apply(Prostate[,1:8],2,mean)
sds <- apply(Prostate[,1:8],2,sd)
for (j in 1:8) {x[,j] <- (x[,j]-means[j])/sds[j]}
for(i in 1:length(lambda)){  
      usemat <-  solve(t(x) %*% x + lambda[i]*diag(x=1,nrow=8)) %*% t(x)  
      betaridge[i,] <- usemat %*% y
      df[i] <- tr(x %*% usemat)
}
# 
# Fig 10.5
#
pdf("prostridge.pdf",height=7,width=7.5)  
matplot(df,betaridge[,1:8],type="l",ylab=expression(hat(beta)),
        xlab="Degrees of Freedom",lty=c(1:6,8,9),col="black",
        xlim=c(0,8),xaxt="n",yaxt="n")
axis(side=1,at=c(0,2,4,6,8))
axis(side=2)
legend("topleft",bty="n",legend=c("lcavol","lweight","age","lbph","svi",
       "lcp","gleason","pgg45"),lty=c(1:6,8,9))
abline(0,0,lty=1,col="grey")
dev.off()
#
# Now the lasso
#
scale <- seq(0.01,1,.005)
x <- as.matrix(Prostate[,1:8])
y <- Prostate[,9]
means <- apply(Prostate[,1:8],2,mean)
sds <- apply(Prostate[,1:8],2,sd)
for (j in 1:8) {x[,j] <- (x[,j]-means[j])/sds[j]}
standProst <- data.frame(y,x)
res <-  l1ce(y~.,data=standProst,bound=scale)
coefmat <- cbind(coef(res)[,2],coef(res)[,3],coef(res)[,4],coef(res)[,5],
           coef(res)[,6],coef(res)[,7],coef(res)[,8],coef(res)[,9])
# 
# Fig 10.7
#
pdf("prostlasso.pdf",height=7,width=7.5)  
matplot(scale,coefmat[,1:8],type="l",xlim=c(0,1),
  ylab=expression(hat(beta)),xlab="Shrinkage Factor",
  lty=c(1:6,8,9),col="black")
legend("topleft",bty="n",legend=c("lcavol","lweight","age","lbph",
  "svi","lcp","gleason","pgg45"),lty=c(1:6,8,9))
abline(0,0,lty=1,col="grey")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Comparison of methods on splits of data. Example after Section 10.6.
#
library(lasso2)
data(Prostate)
library(leaps)
library(rpart)
library(BMA)  
library(lars)  
#
# Randomly sample 67 observations and then form training and test datasets.
#
# PLACE TO START!
#
set.seed(1294)
nsim <- 500 # Number of repeat split datasets to look over
testerrbase <- testerrfull <- testerrridge <- testerrridge5CV <- testerrmodave <- testerrlasso <- testerrCp <- testerrlasso5 <- testerrtree <- dfminOCV <- dfminGCV <- dfminMal <- dfmin5CV <- pesummin <- smin <- smin5 <- NULL
whichcoeff <- matrix(0,nrow=nsim,ncol=8)
#
for(si in 1:nsim){
  scramble <- sample(c(1:97),size=97,replace=F)
  traini <- scramble[1:67]
  testi <- scramble[68:97]
  traindat <- Prostate[traini,]
  means <- apply(traindat[,1:8],2,mean)
  sds <- apply(traindat[,1:8],2,sd)
  testdat <- Prostate[testi,]
  for (j in 1:8){traindat[,j] <- (traindat[,j]-means[j])/sds[j]}
  for (j in 1:8){testdat[,j] <- (testdat[,j]-means[j])/sds[j]}
  testerrbase[si] <- mean((testdat[,9]-mean(testdat[,9]))^2)
  trainfullmod <- lm(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45,
             data=traindat)
  testx <- as.matrix(cbind(1,testdat[,1:8]))
  betahat <- as.matrix(trainfullmod$coeff,nrow=1,ncol=9)
  testfit <- testx %*% betahat
#  testfit2 <- predict.lm(trainfullmod,testdat)# the same as the pevious 3 lines
  testerrfull[si] <- mean((testfit-testdat[,9])^2)
#
# Now let's look at ridge regression
#
  sigma2full <- sum((traindat[,9]-trainfullmod$fitted)^2)/trainfullmod$df.resid
  lambda <- exp(seq(log(.1),log(1000),.1))
  df <- sumsq <- gsumsq <- aveRSS <- AMSE <- penalty <- pesum <- NULL
  betaridge <- matrix(0,nrow=length(lambda),ncol=9)
  x <- as.matrix(cbind(rep(1,dim(traindat)[1]),traindat[,1:8]))
  y <- traindat[,9]
# Set up 5-fold cross-validation datasets
  cvbetaridge <- NULL
  cvscramble <- sample(c(1:67),size=67,replace=F)
  set1 <- cvscramble[1:13]
  set2 <- cvscramble[14:26]
  set3 <- cvscramble[27:39]
  set4 <- cvscramble[40:53]
  set5 <- cvscramble[54:67]
  cvtraindat1 <- traindat[c(set1,set2,set3,set4),]
  cvtraindat2 <- traindat[c(set1,set2,set3,set5),]
  cvtraindat3 <- traindat[c(set1,set2,set4,set5),]
  cvtraindat4 <- traindat[c(set1,set3,set4,set5),]
  cvtraindat5 <- traindat[c(set2,set3,set4,set5),]
  means1 <- apply(cvtraindat1[,1:8],2,mean)
  sds1 <- apply(cvtraindat1[,1:8],2,sd)
  means2 <- apply(cvtraindat2[,1:8],2,mean)
  sds2 <- apply(cvtraindat2[,1:8],2,sd)
  means3 <- apply(cvtraindat3[,1:8],2,mean)
  sds3 <- apply(cvtraindat3[,1:8],2,sd)
  means4 <- apply(cvtraindat4[,1:8],2,mean)
  sds4 <- apply(cvtraindat4[,1:8],2,sd)
  means5 <- apply(cvtraindat5[,1:8],2,mean)
  sds5 <- apply(cvtraindat5[,1:8],2,sd)
  cvtestdat1 <- traindat[set5,1:8]
  cvtestdat2 <- traindat[set4,1:8]
  cvtestdat3 <- traindat[set3,1:8]
  cvtestdat4 <- traindat[set2,1:8]
  cvtestdat5 <- traindat[set1,1:8]
  for (j in 1:8){cvtraindat1[,j] <- (cvtraindat1[,j]-means1[j])/sds1[j]}
  for (j in 1:8){cvtestdat1[,j] <- (cvtestdat1[,j]-means1[j])/sds1[j]}
  for (j in 1:8){cvtraindat2[,j] <- (cvtraindat2[,j]-means2[j])/sds2[j]}
  for (j in 1:8){cvtestdat2[,j] <- (cvtestdat2[,j]-means2[j])/sds2[j]}
  for (j in 1:8){cvtraindat3[,j] <- (cvtraindat3[,j]-means3[j])/sds3[j]}
  for (j in 1:8){cvtestdat3[,j] <- (cvtestdat3[,j]-means3[j])/sds3[j]}
  for (j in 1:8){cvtraindat4[,j] <- (cvtraindat4[,j]-means4[j])/sds4[j]}
  for (j in 1:8){cvtestdat4[,j] <- (cvtestdat4[,j]-means4[j])/sds4[j]}
  for (j in 1:8){cvtraindat5[,j] <- (cvtraindat5[,j]-means5[j])/sds5[j]}
  for (j in 1:8){cvtestdat5[,j] <- (cvtestdat5[,j]-means5[j])/sds5[j]}
  cvx1 <- as.matrix(cbind(rep(1,dim(cvtraindat1)[1]),cvtraindat1[,1:8]))
  cvy1 <- cvtraindat1[,9]
  cvx2 <- as.matrix(cbind(rep(1,dim(cvtraindat2)[1]),cvtraindat2[,1:8]))
  cvy2 <- cvtraindat2[,9]
  cvx3 <- as.matrix(cbind(rep(1,dim(cvtraindat3)[1]),cvtraindat3[,1:8]))
  cvy3 <- cvtraindat3[,9]
  cvx4 <- as.matrix(cbind(rep(1,dim(cvtraindat4)[1]),cvtraindat4[,1:8]))
  cvy4 <- cvtraindat4[,9]
  cvx5 <- as.matrix(cbind(rep(1,dim(cvtraindat5)[1]),cvtraindat5[,1:8]))
  cvy5 <- cvtraindat5[,9]
# Start the lambda loop
  for(i in 1:length(lambda)){ 
# First: leave-one-out cross validation 
      usemat <-  solve(t(x) %*% x + lambda[i]*diag(c(0,rep(1,8)))) %*% t(x)  
      betaridge[i,] <- usemat %*% y
      df[i] <- tr(x %*% usemat)
      Smat <- x %*% usemat
      fitted <- Smat %*% y
      Sii <- diag(Smat)
      sumsq[i] <- sum((y-fitted)^2/(1-Sii)^2)/length(y)
      gsumsq[i] <- sum( (y-fitted)^2/(1-tr(Smat)/length(y))^2 )/length(y)
      aveRSS[i] <- sum((y-fitted)^2)/length(y) 
      penalty[i] <- (length(y)-2*df[i])*sigma2full/length(y)
      AMSE[i] <- aveRSS[i] - penalty[i]
# Second: K=5 cross validation
      pesum[i] <- 0
      cvusemat1<-solve(t(cvx1)%*%cvx1+lambda[i]*diag(c(0,rep(1,8))))%*%t(cvx1) 
      cvusemat2<-solve(t(cvx2)%*%cvx2+lambda[i]*diag(c(0,rep(1,8))))%*%t(cvx2) 
      cvusemat3<-solve(t(cvx3)%*%cvx3+lambda[i]*diag(c(0,rep(1,8))))%*%t(cvx3) 
      cvusemat4<-solve(t(cvx4)%*%cvx4+lambda[i]*diag(c(0,rep(1,8))))%*%t(cvx4) 
      cvusemat5<-solve(t(cvx5)%*%cvx5+lambda[i]*diag(c(0,rep(1,8))))%*%t(cvx5)  
      cvbetaridge1 <- cvusemat1 %*% cvy1
      cvbetaridge2 <- cvusemat2 %*% cvy2
      cvbetaridge3 <- cvusemat3 %*% cvy3
      cvbetaridge4 <- cvusemat4 %*% cvy4
      cvbetaridge5 <- cvusemat5 %*% cvy5
      cvfitted1 <- as.matrix(cbind(rep(1,length(set5)),cvtestdat1)) %*% 
       cvbetaridge1; pesum[i] <- pesum[i] + sum((traindat[set5,9]-cvfitted1)^2)
      cvfitted2 <- as.matrix(cbind(rep(1,length(set4)),cvtestdat2)) %*% 
       cvbetaridge2; pesum[i] <- pesum[i] + sum((traindat[set4,9]-cvfitted2)^2)
      cvfitted3 <- as.matrix(cbind(rep(1,length(set3)),cvtestdat3)) %*% 
       cvbetaridge3; pesum[i] <- pesum[i] + sum((traindat[set3,9]-cvfitted3)^2)
      cvfitted4 <- as.matrix(cbind(rep(1,length(set2)),cvtestdat4)) %*% 
       cvbetaridge4; pesum[i] <- pesum[i] + sum((traindat[set2,9]-cvfitted4)^2)
      cvfitted5 <- as.matrix(cbind(rep(1,length(set1)),cvtestdat5)) %*% 
       cvbetaridge5; pesum[i] <- pesum[i] + sum((traindat[set1,9]-cvfitted5)^2)
	pesum[i] <- pesum[i]/67
  }
  lamminOCV <- lambda[which.min(sumsq)]
  lamminGCV <- lambda[which.min(gsumsq)]
  lamminMal <- lambda[which.min(AMSE)]
  lammin5CV <- lambda[which.min(pesum)]
  dfminOCV[si] <- df[which.min(sumsq)]
  dfminGCV[si] <- df[which.min(gsumsq)]
  dfminMal[si] <- df[which.min(AMSE)]
  dfmin5CV[si] <- df[which.min(pesum)]  
  pesummin[si] <- min(pesum)
#
  betahatridge <- as.matrix(betaridge[which.min(AMSE),1:9])
  testfitridge <- testx %*% betahatridge
  testerrridge[si] <- mean((testfitridge-testdat[,9])^2)
  betahatridge5CV <- as.matrix(betaridge[which.min(pesum),1:9])
  testfitridge5CV <- testx %*% betahatridge5CV
  testerrridge5CV[si] <- mean((testfitridge5CV-testdat[,9])^2)
# Now model averaging
  fullmod <- lm(lpsa~.,data=traindat)
  x <- model.matrix(fullmod)[,-1]
  bicmod <- bic.glm(x,y,glm.family="gaussian")
  testfitmodave <- testx %*% bicmod$postmean
  testerrmodave[si] <- mean((testfitmodave-testdat[,9])^2)
#
# Exhaustive search using Cp (equivalent to AIC)
#
  cp <- leaps(x,y,method="Cp")
  whichcoeff[si,] <- cp$which[which.min(cp$Cp),]
  xred <- x[,whichcoeff[si,]==T]
  redmod <- lm(y~xred)
  testfitmodCp <- testx[,whichcoeff[si,]==T] %*% redmod$coeff
  testerrCp[si] <- mean((testfitmodCp-testdat[,9])^2)
#
# Lasso with both K=n and K=5 cross-validation
#
  lasscv <- cv.lars(x=x,y=y,type="lasso",trace=F,K=length(y),plot.it=F)
  lasscv5 <- cv.lars(x=x,y=y,type="lasso",trace=F,K=5,plot.it=F)
  smin[si] <- lasscv$index[which.min(lasscv$cv)]
  smin5[si] <- lasscv5$index[which.min(lasscv5$cv)]
  ldatframe <- data.frame(y,x)
  lassomod <-  l1ce(ldatframe$y~ldatframe$lcavol+ldatframe$lweight+
               ldatframe$age+ldatframe$lbph+ldatframe$svi+ldatframe$lcp+
               ldatframe$gleason+ldatframe$pgg45,bound=smin[si],absolute.t=F)
  lassomod5 <-  l1ce(ldatframe$y~ldatframe$lcavol+ldatframe$lweight+
                ldatframe$age+ldatframe$lbph+ldatframe$svi+ldatframe$lcp+
                ldatframe$gleason+ldatframe$pgg45,bound=smin5[si],absolute.t=F)
  testfitlasso <- testx %*% lassomod$coeff
  testfitlasso5 <- testx %*% lassomod5$coeff
  testerrlasso[si] <- mean((testfitlasso-testdat[,9])^2)
  testerrlasso5[si] <- mean((testfitlasso5-testdat[,9])^2)
#
# Regression tree 
#
  treefit <- rpart(lpsa~.,data=traindat,method="anova",
             control=list(minsplit=5,minbucket=3,cp=.001))
#
  prunedtreefit <- prune(treefit,cp=treefit$cptable[
                   which.min(treefit$cptable[,"xerror"]),"CP"])
  dftestx <- as.data.frame(testx)
  testfittree <- predict(prunedtreefit,newdata=dftestx)
  testerrtree[si] <- mean((testfittree-testdat[,9])^2)
#
#  cat("Simulation Test error for baseline, full model, ridge regression, 5CV ridge regression, BMA Exhaustive Cp lasso tree= \n",si,testerrbase[si],testerrfull[si],testerrridge[si],testerrridge5CV[si],testerrmodave[si],testerrCp[si],testerrlasso[si],testerrlasso5[si],testerrtree[si],"\n")
  if ((floor(si/50)==ceiling(si/50))) cat("Simulation no = ",si,"\n")
}
#
# Summaries in Table 10.3
#
cat("Test base mean and sd = ",mean(testerrbase),sd(testerrbase),"\n")
cat("Test full mean and sd = ",mean(testerrfull),sd(testerrfull),"\n")
cat("Test ridge mean and sd = ",mean(testerrridge),sd(testerrridge),"\n")
cat("Test ridge mean and sd 5-fold= ",mean(testerrridge5CV),
    sd(testerrridge5CV),"\n")
cat("Test BMA mean and sd = ",mean(testerrmodave),sd(testerrmodave),"\n")
cat("Test Cp mean and sd = ",mean(testerrCp),sd(testerrCp),"\n")
cat("Test LASSO mean and sd = ",mean(testerrlasso),sd(testerrlasso),"\n")
cat("Test LASSO 5-fold mean and sd = ",mean(testerrlasso5),sd(testerrlasso5),
    "\n")
cat("Test TREE mean and sd = ",mean(testerrtree),sd(testerrtree),"\n")
#
# Exhaustive search summary
#
psum <- 0
sequ=unique(whichcoeff, margin=1)
for (i in 1:nrow(sequ)){
   tempsum <- 0
   for (j in 1:nsim){
       temp <- sequ[i,]==whichcoeff[j,]
       tempsum <- tempsum + ifelse(sum(temp)==8,1,0)
   }
   psum <- psum + tempsum/nsim
   cat("Proportion of ",sequ[i,]," is ",tempsum/nsim,"\n")
}
#
# Table 10.4 percentages of different models selected
#
cat("Total Prob = ",psum,"\n")
#
# Fig 10.8(a)
#
pdf("RidgedfminOCV2.pdf") 
hist(dfminOCV-1,nclass=5,xlab="Optimal DF",main="",xlim=c(min(dfminOCV-1,
       dfmin5CV-1,dfminGCV-1,dfminMal-1),max(dfminOCV-1,dfmin5CV-1,dfminGCV-1,
       dfminMal-1)))
dev.off()
# Minimizing K=5 CV parameters -- not in book
hist(dfmin5CV-1,nclass=5,xlab="Optimal DF",main="",xlim=c(min(dfminOCV-1,
       dfmin5CV-1,dfminGCV-1,dfminMal-1),max(dfminOCV-1,dfmin5CV-1,dfminGCV-1,
       dfminMal-1)))
# Minimizing GCV parameters -- not in book
hist(dfminGCV-1,nclass=5,xlab="Optimal DF",main="",
       xlim=c(min(dfminOCV-1,dfmin5CV-1,dfminGCV-1,dfminMal-1),
       max(dfminOCV-1,dfmin5CV-1,dfminGCV-1,dfminMal-1)))
# 
# Fig 10.8(b)
#
pdf("RidgedfminMal2.pdf") 
hist(dfminMal-1,nclass=5,xlab="Optimal DF",main="",
     xlim=c(min(dfminOCV-1,dfmin5CV-1,dfminGCV-1,dfminMal-1),
     max(dfminOCV-1,dfmin5CV-1,dfminGCV-1,dfminMal-1)))
dev.off()
#
# Fig 10.8(c)
#
pdf("RidgedfminMalOCV2.pdf") 
plot(dfminOCV-1,dfminMal-1,xlim=c(4.7,8),ylim=c(4.7,8),
     xlab="OCV Optimal DF",ylab="Cp Optimal DF")
abline(0,1)
dev.off()
#
# Plot of the predictive risk curves versus deg of freedom
#
# Fig 10.9
#
pdf("ridgecv2.pdf",height=4.5,width=5.5) 
par(mfrow=c(1,1))
plot(df-1,sumsq,type="l",ylab="Error",xlab="Degrees of Freedom",
       ylim=c(min(aveRSS,sumsq),max(aveRSS,sumsq)),xaxt="n",yaxt="n",
       xlim=c(-.02,8.02))
axis(side=1,at=c(0,2,4,6,8))
axis(side=2)
lines(df-1,gsumsq,lty=2)
lines(df-1,AMSE+sigma2full,lty=3)
lines(df-1,aveRSS,lty=4)
lines(df-1,pesum,lty=6)
legend("topright",bty="n",legend=c("Ordinary Cross-Validation",
       "Generalized Cross-Validation","Mallows Cp","Residual Sum of Squares",
       "5-Fold Cross-Validation"),lty=c(1:4,6))
vlim <- .85
lines(x=c(dfminOCV[si]-1,dfminOCV[si]-1),y=c(min(aveRSS),vlim),lty=1)
lines(x=c(dfminGCV[si]-1,dfminGCV[si]-1),y=c(min(aveRSS),vlim),lty=2)
lines(x=c(dfminMal[si]-1,dfminMal[si]-1),y=c(min(aveRSS),vlim),lty=3)
lines(x=c(dfmin5CV[si]-1,dfmin5CV[si]-1),y=c(min(aveRSS),vlim),lty=6)
dev.off()
#
# Lasso plots
lmax <- max(lasscv$cv+lasscv$cv.error,lasscv5$cv+lasscv5$cv.error)
lmin <- min(lasscv$cv-lasscv$cv.error,lasscv5$cv-lasscv5$cv.error)
#
# Fig 10.10(a)
#
pdf("lassoSEOCV2.pdf") 
plot(lasscv$cv~lasscv$index,type="n",ylab="Error",xlab="Scaling Factor",
     ylim=c(lmin,lmax))
lines(lasscv$cv~lasscv$index)
abline(v=min(lasscv$index[which.min(lasscv$cv)]))
for (i in 1:length(lasscv$cv)){
    lines(x=c(lasscv$index[i],lasscv$index[i]),
     y=c(lasscv$cv[i]+lasscv$cv.error[i],lasscv$cv[i]-lasscv$cv.error[i]))
}
dev.off()
# 
# Fig 10.10(b)
#
pdf("lassoSE5CV2.pdf") 
plot(lasscv5$cv~lasscv5$index,type="n",ylab="Error",xlab="Scaling Factor",
     ylim=c(lmin,lmax))
lines(lasscv5$cv~lasscv5$index)
abline(v=min(lasscv5$index[which.min(lasscv5$cv)]))
for (i in 1:length(lasscv5$cv)){
    lines(x=c(lasscv5$index[i],lasscv5$index[i]),
     y=c(lasscv5$cv[i]+lasscv5$cv.error[i],lasscv5$cv[i]-lasscv5$cv.error[i]))
}
dev.off()
#
# BMA plot now
#
#---------------------------------------------------------------------------
# I changed the original imageplot.bma.R a bit to give:
#
myimageplot.bma <- function (bma.out, color = c("red", "blue", "#FFFFD5"), 
    order = c("input","probne0", "mds"), ...) 
{
    clr <- color
    if (length(color) == 1) {
        if (color == "default") 
            clr <- c("#FF0000", "#FF0000", "#FFFFD5")
        if (color == "blackandwhite") 
            clr <- c("black", "black", "white")
    }
    keep.mar <- par(mar = c(5, 6, 4, 2) + 0.1)
    nmodel <- nrow(bma.out$which)
    which <- bma.out$which
    probne0 <- bma.out$probne0
    if (class(bma.out) == "bic.surv") 
        mle <- bma.out$mle
    else mle <- bma.out$mle[, -1, drop = FALSE]
    nvar <- ncol(mle)
    rownms <- bma.out$namesx
    if (ifelse(!is.null(bma.out$factor.type), bma.out$factor.type, 
        FALSE)) {
        which <- matrix(NA, ncol = nvar, nrow = nmodel)
        probne0 <- rep(NA, times = nvar)
        rownms <- rep(NA, times = nvar)
        assign <- bma.out$assign
        offset <- 1
        if (class(bma.out) == "bic.surv") 
            offset <- 0
        assign[[1]] <- NULL
        for (i in 1:length(assign)) {
            probne0[assign[[i]] - offset] <- bma.out$probne0[i]
            which[, assign[[i]] - offset] <- bma.out$which[, 
                i]
            nm <- names(bma.out$output.names)[i]
            if (!is.na(bma.out$output.names[[i]][1])) 
                nm <- paste(nm, bma.out$output.names[[i]][-1], 
                  sep = ".")
            rownms[assign[[i]] - offset] <- nm
        }
    }
    ordr.type <- match.arg(order)
    if (ordr.type == "probne0") 
        ordr <- order(-probne0)
    else if (ordr.type == "mds") {
        postprob.rep <- matrix(bma.out$postprob, ncol = nvar, 
            nrow = nmodel)
        k11 <- t(which + 0) %*% ((which + 0) * postprob.rep)
        k00 <- t(1 - which) %*% ((1 - which) * postprob.rep)
        k01 <- t(which + 0) %*% ((1 - which) * postprob.rep)
        k10 <- t(1 - which) %*% ((0 + which) * postprob.rep)
        ktau <- 4 * (k00 * k11 - k01 * k10)
        dissm <- 1 - abs(ktau)
        diag(dissm) <- 0
        ordr <- order(as.vector(cmdscale(dissm, k = 1)))
    }
    else ordr <- 1:nvar
    ordr <- rev(ordr)
    postprob <- bma.out$postprob
    which <- which[, ordr, drop = FALSE]
    mle <- mle[, ordr, drop = FALSE]
    rownms <- rownms[ordr]
    color.matrix <- (which) * (2 - (mle > 0)) + 3 * (!which)
    par(las = 1)
    image(c(0, cumsum(postprob)), 1:nvar, color.matrix, col = clr, 
        xlab = "Model Number", ylab = "", xaxt = "n", yaxt = "n", 
        xlim = c(0, 1), main = "", ...)
    xat <- (cumsum(postprob) + c(0, cumsum(postprob[-nmodel])))/2
    axis(1, at = xat, labels = 1:nmodel, ...)
    axis(2, at = 1:nvar, labels = rownms, ...)
    par(mar = keep.mar)
}
# 
# Fig 10.11
#
pdf("BMAimage2.pdf") 
# +ve coef, -ve coef, not in mod
myimageplot.bma(bicmod,color=c("black","blue","grey"),yaxt="n") 
axis(side=2,at=seq(1,8),labels=c("pgs45","gleason","svi","lcp",
     "lbph","age","lweight","lcanvol"))
dev.off()
