library(HWEBayes)
# Four allele example
bvec0 <- c(1,1,1,1)
nvec0 <- rep(0,10)
# First sample from the prior under the null
priorsampH0 <- DirichSampHWE(nvec0,bvec0,nsim=1000)
par(mfrow=c(2,2))
hist(priorsampH0$pvec[,1],xlab=expression(p[1]),main="")
hist(priorsampH0$pvec[,2],xlab=expression(p[2]),main="")
hist(priorsampH0$pvec[,3],xlab=expression(p[3]),main="")
hist(priorsampH0$pvec[,4],xlab=expression(p[4]),main="")
#
data(DiabRecess)
nvec <- DiabRecess
postsampH0 <- DirichSampHWE(nvec,bvec0,nsim=1000)
MLE4 <- HWEmodelsMLE(nvec)
par(mfrow=c(2,2))
hist(postsampH0$pvec[,1],xlab=expression(p[1]),main="")
abline(v=MLE4$qhat[1],col="red")
hist(postsampH0$pvec[,2],xlab=expression(p[2]),main="")
abline(v=MLE4$qhat[2],col="red")
hist(postsampH0$pvec[,3],xlab=expression(p[3]),main="")
abline(v=MLE4$qhat[3],col="red")
hist(postsampH0$pvec[,4],xlab=expression(p[4]),main="")
abline(v=MLE4$qhat[4],col="red")
#
bvec1 <- rep(1,10)
nvec1 <- rep(0,10)
priorsampH1sat <- DirichSampSat(nvec=nvec1,bvec1,nsim=1000)
par(mfrow=c(2,3))
hist(priorsampH1sat$pvec[,1],xlab=expression(p[11]),main="")
hist(priorsampH1sat$pvec[,2],xlab=expression(p[12]),main="")
hist(priorsampH1sat$pvec[,3],xlab=expression(p[22]),main="")
hist(priorsampH1sat$pmarg[,1],xlab=expression(p[1]),main="")
hist(priorsampH1sat$pmarg[,2],xlab=expression(p[2]),main="")
hist(priorsampH1sat$fixind[,2,1],xlab=expression(f[12]),main="")
#
# Sample from the saturated posterior for the 4 allele data
postsampH1sat <- DirichSampSat(nvec,bvec1,nsim=1000)
par(mfrow=c(2,3))
hist(postsampH1sat$pvec[,1],xlab=expression(p[11]),main="")
abline(v=MLE4$phat[1,1],col="red")
hist(postsampH1sat$pvec[,2],xlab=expression(p[12]),main="")
abline(v=MLE4$phat[1,2],col="red")
hist(postsampH1sat$pvec[,3],xlab=expression(p[22]),main="",xlim=c(0,.3))
abline(v=MLE4$phat[2,2],col="red")
hist(postsampH1sat$pmarg[,1],xlab=expression(p[1]),main="")
abline(v=MLE4$qhat[1],col="red")
hist(postsampH1sat$pmarg[,2],xlab=expression(p[2]),main="")
abline(v=MLE4$qhat[2],col="red")
hist(postsampH1sat$fixind[,2,1],xlab=expression(f[12]),main="")
abline(v=MLE4$fixind[1,2],col="red")
#
# Single f example
bvec <- c(1,1,1,1)
# Find the parameters for the prior for f
init <- c(-3,log(1.1)) # Good starting values needed 
lampr <- LambdaOptim(nsim=10000,bvec=bvec,f1=0,f2=0.26,p1=0.5,p2=0.95,init)
nsim <- 100
postsampf1 <- SinglefReject(nsim,bvec,lambdamu=lampr$lambdamu,
                            lambdasd=lampr$lambdasd,nvec)
par(mfrow=c(2,3))
hist(postsampf1$psamp[,1],xlab=expression(p[1]),main="")
abline(v=MLE4$fqhat[1],col="red")
hist(postsampf1$psamp[,2],xlab=expression(p[2]),main="")
abline(v=MLE4$fqhat[2],col="red")
hist(postsampf1$psamp[,3],xlab=expression(p[3]),main="")
abline(v=MLE4$fqhat[3],col="red")
hist(postsampf1$psamp[,4],xlab=expression(p[4]),main="")
abline(v=MLE4$fqhat[4],col="red")
hist(postsampf1$fsamp,xlab="f",main="")
abline(v=MLE4$fsingle,col="red")
#
PrnH0 <- DirichNormHWE(nvec,bvec0)
PrnH1sat <- DirichNormSat(nvec,bvec1)
BFH0H1sat <- PrnH0/PrnH1sat
#
alpha <- rep(1,4) 
# First simulate from a normal proposal using mean vector and covariance
# matrix from a WinBUGS run
gmu <- c(-0.4633092,0.3391625,0.3397936,-3.5438008)
gsigma <- matrix(c(
0.07937341,0.02819656,0.02766583,0.04607996,
0.02819656,0.07091320,0.04023827,0.01657028,
0.02766583,0.04023827,0.07042278,0.01752266,
0.04607996,0.01657028,0.01752266,0.57273683),nrow=4,ncol=4)
est1 <- HWEImportSamp(nsim=5000,nvec,ischoice=1,lambdamu=lampr$lambdamu,
             lambdasd=lampr$lambdasd,alpha=alpha,gmu,gsigma)
# Now let's evaluate using the prior
est2 <- HWEImportSamp(nsim=20000,nvec,ischoice=2,lambdamu=lampr$lambdamu,
             lambdasd=lampr$lambdasd,alpha=alpha,gmu,gsigma)