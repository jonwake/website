#
#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 12
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Example at the end of Section 12.2
#
library(lasso2)
data(Prostate)
attach(Prostate)
library(mgcv)
gammod <- gam(lpsa~s(x=lcavol,k=7,fx=F,bs="cr",m=2)+
      s(x=lweight,k=7,fx=F,bs="cr",m=2)+s(x=age,k=7,fx=F,bs="cr",m=2)+
      s(x=lbph,k=7,fx=F,bs="cr",m=2)+s(x=lcp,k=7,fx=F,bs="cr",m=2)+
      svi+gleason+s(x=pgg45,k=7,fx=F,bs="cr",m=2),method="GCV.Cp")
#
# Fig 12.1(a)
#
pdf("ProstateSmoothFig1.pdf",h=4,w=4)
plot(gammod,select=1,ylab="Fitted Cubic Spline",xlab="log(can vol)",shade=T,seWithMean=T)
dev.off()
#
# Fig 12.1(b)
#
pdf("ProstateSmoothFig2.pdf",h=4,w=4)
plot(gammod,select=2,ylab="Fitted Cubic Spline",xlab="log(weight)",shade=T,seWithMean=T)
dev.off()
#
# Fig 12.1(c)
#
pdf("ProstateSmoothFig3.pdf",h=4,w=4)
plot(gammod,select=3,ylab="Fitted Cubic Spline",xlab="Age",shade=T,seWithMean=T)
dev.off()
#
# Fig 12.1(d)
#
pdf("ProstateSmoothFig4.pdf",h=4,w=4)
plot(gammod,select=4,ylab="Fitted Cubic Spline",xlab="log(BPH)",shade=T,seWithMean=T)
dev.off()
#
# Fig 12.1(e)
#
pdf("ProstateSmoothFig5.pdf",h=4,w=4)
plot(gammod,select=5,ylab="Fitted Cubic Spline",xlab="log(cap pen)",shade=T,seWithMean=T)
dev.off()
#
# Fig 12.1(f)
#
pdf("ProstateSmoothFig6.pdf",h=4,w=4)
plot(gammod,select=6,ylab="Fitted Cubic Spline",xlab="PGS45",shade=T,seWithMean=T)
dev.off()
# 
# Prostate Cancer Example in Section 12.3.2
# 
# Linear model first
#
linmod <- lm(lpsa~.,data=Prostate)
wgrid <- seq(min(Prostate$lweight),max(Prostate$lweight),length=25)
cgrid <- seq(min(Prostate$lcavol),max(Prostate$lcavol),length=25)
zvals <- matrix(0,nrow=length(cgrid),ncol=length(wgrid))
for (i in 1:length(cgrid)){
    for (j in 1:length(wgrid)){
    	zvals[i,j] <- coef(linmod)[3]*wgrid[i] + coef(linmod)[2]*cgrid[j]
    }
}
#
# Fig 12.2(a)
#
pdf("Prostate-2dpersp.pdf",h=8,w=8)
persp(x=wgrid,y=cgrid,z=zvals,theta=35,phi=25,xlab="log(can vol)",
      ylab="log(weight)",zlab="Fitted Linear")
dev.off()
#
# Now TPRS
#
gammod <- gam(lpsa~s(lcavol,lweight,bs="tp")+s(x=age,k=7,fx=F,bs="cr",m=2)+
          s(x=lbph,k=7,fx=F,bs="cr",m=2)+s(x=lcp,k=7,fx=F,bs="cr",m=2)+
          svi+gleason+s(x=pgg45,k=7,fx=F,bs="cr",m=2),method="GCV.Cp")
#
# Fig 12.2(b)
#
pdf("TPRSprost1.pdf",h=8,w=8)
plot(gammod,select=1,pers=T,xlab="log(can vol)",ylab="log(weight)",
     main="Fitted TPRS",theta=35,phi=25)
dev.off()
#
# Now tensor products
#
gammod2 <- gam(lpsa~te(lcavol,lweight,k=c(6,6))+s(x=age,k=7,fx=F,bs="cr",m=2)+
    s(x=lbph,k=7,fx=F,bs="cr",m=2)+s(x=lcp,k=7,fx=F,bs="cr",m=2)+svi+gleason+
    s(x=pgg45,k=7,fx=F,bs="cr",m=2),method="GCV.Cp")
#
# Fig 12.2(c)
#
pdf("TensProdprost1.pdf",h=8,w=8)
plot(gammod2,select=1,pers=T,xlab="log(can vol)",ylab="log(weight)",
     main="Fitted Tensor Product",theta=35,phi=25)
dev.off()
summary(gammod2) # For effective degrees of freedom of tensor product = 12.4
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Varying coefficient model for the ethanol data -- example at the
# end of Section 12.6
#
library(SemiPar)
data(ethanol)
attach(ethanol)
nint <- 9
quantE <- midE <- NULL
quantE[1] <- min(E)-.01
quantE[2:nint] <- quantile(E,p=seq(1:(nint-1))/nint)
quantE[nint+1] <- max(E)+.01
#
# Fig 12.3
#
pdf("ethfig1.pdf")
par(mfrow=c(3,3))
betas <- matrix(0,nrow=nint,ncol=2)
for (i in 1:nint){
    y <- NOx[E>quantE[i]&E<=quantE[i+1]]
    x <- C[E>quantE[i]&E<=quantE[i+1]]
    midE[i] <- quantE[i] +(quantE[i+1]-quantE[i])/2
    plot(y~x,ylim=c(min(NOx),max(NOx)),xlim=c(min(C),max(C)),xlab="C",ylab="NOx")
    mod <- lm(y~x)
    abline(coef(mod))
    betas[i,] <- coef(mod)
    cat("Quintile n beta: ",i,length(y),coef(mod),"\n")
}
dev.off()
inter <- 1.29
Elims <- c(.51,1.25)
ylims1 <- c(-2.5,2)
ylims2 <- c(-.03,.15)
#
modgam2 <- gam(NOx~s(E,bs="cr")+s(E,by=C,bs="cr"))
summary(modgam2) # effective degrees of freedom of 6.4 and 4.7 on the intercept and slope smooths
#
# Fig 12.4(a)
#
pdf("ethvarycoeffig1.pdf",h=4,w=4)
par(mfrow=c(1,1))
plot.gam(modgam2,select=1,se=F,ylab="Intercept(E)",ylim=ylims1,xlim=Elims)
points(midE,adjust)
dev.off()
#
# Fig 12.4(b)
#
pdf("ethvarycoeffig2.pdf",h=4,w=4)
par(mfrow=c(1,1))
plot.gam(modgam2,select=2,se=F,ylab="Slope(E)",ylim=ylims2,xlim=Elims)
points(midE,betas[,2])
dev.off()
#
# Fitted surface
#
Eval <- seq(min(E),max(E),length=25)
Cval <- seq(min(C),max(C),length=25)
fitf <- matrix(0,nrow=length(Eval),ncol=length(Cval))
newd <- data.frame(C=seq(min(C),max(C),length=100),E=seq(min(E),max(E),length=100))
for (i in 1:length(Eval)){
    for (j in 1:length(Cval)){
    	newd <- data.frame(C=Cval[j],E=Eval[i])
    	fitf[i,j] <- predict.gam(modgam2,newdata=newd)
    }
}
#
# Fig 12.5
#
pdf("ethvaryingimage.pdf",h=5,w=5)
par(mfrow=c(1,1))
image(Eval,Cval,fitf,xlab="E",ylab="C",col=gray((0:32)/32))
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Regression trees
#
x1 <- seq(0,1)
x2 <- seq(0,1)
#
# Fig 12.6(a)
#
pdf("regtreefig0.pdf",h=5,w=5) 
par(pty="s")
plot(x1,x2,type="n",xlab="",ylab="",axes=F)
lines(x=c(0,1),y=c(0,0))
lines(x=c(0,0),y=c(0,1))
lines(x=c(0,1),y=c(1,1))
lines(x=c(1,1),y=c(0,1))
axis(side=1,at=c(0.5),labels=expression(x[1]),tick=F,cex.axis=2)
axis(side=2,at=c(0.5),labels=expression(x[2]),tick=F,cex.axis=2)
lines(x=c(0,.25),y=c(0,.6))
lines(x=c(0,.5),y=c(.7,.5))
lines(x=c(0.4,.6),y=c(1,0))
lines(x=c(.54,1),y=c(.3,.7))
dev.off()
#
# Fig 12.6(b)
#
pdf("regtreefig1.pdf",h=5,w=5) 
par(pty="s")
plot(x1,x2,type="n",xlab="",ylab="",axes=F)
axis(side=1,at=c(0.5),labels=expression(x[1]),tick=F,cex.axis=2.0)
 axis(side=2,at=c(0.5),labels=expression(x[2]),tick=F,cex.axis=2.0)
lines(x=c(0,1),y=c(0,0))
lines(x=c(0,0),y=c(0,1))
lines(x=c(0,1),y=c(1,1))
lines(x=c(1,1),y=c(0,1))
lines(x=c(0,.3),y=c(.6,.6))
lines(x=c(0.3,.9),y=c(.4,.4))
lines(x=c(0.3,.9),y=c(.75,.75))
lines(x=c(0.3,.3),y=c(.4,.75))
lines(x=c(0.9,.9),y=c(0,.9))
lines(x=c(0.5,1),y=c(.9,.9))
lines(x=c(0.5,.5),y=c(0.9,1))
dev.off()
#
# Fig 12.7
#
pdf("regtreefig2.pdf")
par(pty="s")
x1 <- seq(0,1)
x2 <- seq(0,1)
plot(x1,x2,type="n",xlab=expression(x[1]),ylab=expression(x[2]),axes=F)
lines(x=c(0,1),y=c(0,0))
lines(x=c(0,0),y=c(0,1))
lines(x=c(0,1),y=c(1,1))
lines(x=c(1,1),y=c(0,1))
t1 <- .3
lines(x=c(0,1),y=c(t1,t1))
t2 <- .2
lines(x=c(t2,t2),y=c(0,t1))
t3 <- .8
lines(x=c(t3,t3),y=c(t1,1))
t4 <- .7
lines(x=c(0,t3),y=c(t4,t4))
axis(2,at=c(t1,t4),labels=c(expression(t[1]),expression(t[4])),las=1,tick=F)
axis(1,at=c(t2,t3),labels=c(expression(t[2]),expression(t[3])),las=1,tick=F)
text(x=t2/2,y=t1/2,label=expression(R[1]))
text(x=t2+(1-t2)/2,y=t1/2,label=expression(R[2]))
text(x=t3+(1-t3)/2,y=t1+(1-t1)/2,label=expression(R[3]))
text(x=t3/2,y=t1+(t4-t1)/2,label=expression(R[4]))
text(x=t3/2,y=t4+(1-t4)/2,label=expression(R[5]))
dev.off()
#
# Fig 12.8
#
pdf("regtreefig3.pdf",h=6,w=6) 
x1 <- seq(0,1)
x2 <- seq(0,1)
par(pty="s")
plot(x1,x2,type="n",xlab="",ylab="",axes=F)
height <- .25
lev1height <- .95
lev2height <- lev1height - height
lev3height <- lev1height - 2*height
lev4height <- lev1height - 3*height
vline1 <- 0
vline2 <- .2
vline3 <- .4
vline4 <- .5
vline5 <- .6
vline6 <- .65
vline7 <- .8
vline8 <- 1
lines(x=c(vline2,vline7),y=c(lev1height,lev1height))
lines(x=c(vline7,vline7),y=c(lev1height,lev2height))
lines(x=c(vline2,vline2),y=c(lev1height,lev2height))
lines(x=c(vline1,vline3),y=c(lev2height,lev2height))
lines(x=c(vline5,vline8),y=c(lev2height,lev2height))
lines(x=c(vline1,vline1),y=c(lev2height,lev3height))
lines(x=c(vline3,vline3),y=c(lev2height,lev3height))
lines(x=c(vline5,vline5),y=c(lev2height,lev3height))
lines(x=c(vline8,vline8),y=c(lev2height,lev3height))
lines(x=c(vline4,vline4),y=c(lev4height,lev3height))
lines(x=c(vline7,vline7),y=c(lev4height,lev3height))
lines(x=c(vline4,vline7),y=c(lev3height,lev3height))
lines(x=c(vline2+(vline7-vline2)/2,vline2+(vline7-vline2)/2),y=c(lev1height-.01,lev1height+.01))
eps <- .04
eps2 <- .09
cexno <- .7
text(x=vline2+(vline7-vline2)/2,y=lev1height+eps,label=expression(x[2]<=t[1]),
     family="Arial",cex=cexno)
text(x=vline2-eps2,y=lev2height+eps,label=expression(x[1]<=t[2]),family="Arial",cex=cexno)
text(x=vline8-eps2,y=lev2height+eps,label=expression(x[1]<=t[3]),family="Arial",cex=cexno)
text(x=vline5-eps2,y=lev3height+eps,label=expression(x[2]<=t[4]),family="Arial",cex=cexno)
text(x=vline1,y=lev3height-eps,label=expression(R[1]),cex=cexno,font=3)
text(x=vline3,y=lev3height-eps,label=expression(R[2]),cex=cexno,font=3)
text(x=vline4,y=lev4height-eps,label=expression(R[4]),cex=cexno,font=3)
text(x=vline7,y=lev4height-eps,label=expression(R[5]),cex=cexno,font=3)
text(x=vline8,y=lev3height-eps,label=expression(R[3]),cex=cexno,font=3)
dev.off()
#
#------------------------------------------------------------------
# Function to create fitted surface for Fig 12.9
#
predict.example <- function(x1,x2){
	t1 <- .3; t2 <- .2; t3 <- .8; t4 <- .7
        R1 <- -3; R2 <- -1; R3 <- 0; R4 <- 1; R5 <- 3
	if (x2 <= t1){
	   if (x1 <= t2) predict.example <- R1
	   else if (x1 > t2) predict.example <- R2
	}
	else if (x2 > t1){
	   if (x1 > t3) predict.example <- R3
	   else if (x1 <= t3){
		if (x2 <= t4) predict.example <- R4
		else if (x2 > t4) predict.example <- R5
           }
        }
}
x1val <- seq(0,1,length=100)
x2val <- seq(0,1,length=100)
fitf <- matrix(0,nrow=length(x1val),ncol=length(x2val))
for (i in 1:length(x1val)){
    for (j in 1:length(x2val)){
        fitf[i,j] <- predict.example(x1val[i],x2val[j])
    }
}
library(lattice)
library(reshape)
M <- melt(data.frame(x1val,fitf,check.names=FALSE),id=1,variable="x2val")
M2 <- M
M2$x2val <- as.numeric(M$x2val)/100
#
# Fig 12.9
#
pdf("regtreefig4.pdf",h=6,w=6) 
par(mfrow=c(1,1))
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(fitf~x1val*x2val,data=M2,xlab=expression(x[1]),ylab=expression(x[2]),zlab="y",
screen=list(z=30,x=-60),colkey=F,col="grey",scales=list(arrows=F))
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Regreesion tree prostate data example in Section 12.7
#
#-------------------------------------------------------------------------------
# My version of function to plot CV versus complexity in Fig 12.10
#
myplotcp <- function (x, minline = TRUE, lty = 3, col = 1, upper = c("size", 
    "splits", "none"), ...) 
{
    dots <- list(...)
    if (!inherits(x, "rpart")) 
        stop("Not legitimate rpart object")
    upper <- match.arg(upper)
    p.rpart <- x$cptable
    if (ncol(p.rpart) < 5L) 
        stop("cptable does not contain cross-validation results")
    xstd <- p.rpart[, 5L]
    xerror <- p.rpart[, 4L]
    nsplit <- p.rpart[, 2L]
    ns <- seq_along(nsplit)
    cp0 <- p.rpart[, 1L]
    cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
    if (!"ylim" %in% names(dots)) 
        dots$ylim <- c(min(xerror - xstd) - 0.1, max(xerror + 
            xstd) + 0.1)
    do.call(plot, c(list(ns, xerror, axes = FALSE, xlab = "Complexity", 
        ylab = "CV Estimate", type = "o"), dots))
    box()
    axis(2, ...)
    segments(ns, xerror - xstd, ns, xerror + xstd)
    axis(1L, at = ns, lab = as.character(signif(cp, 2L)), ...)
    switch(upper, size = {
        axis(3L, at = ns, lab = as.character(nsplit + 1), ...)
        mtext("Size of Tree", side = 3, line = 3)
    }, splits = {
        axis(3L, at = ns, lab = as.character(nsplit), ...)
        mtext("number of splits", side = 3, line = 3)
    }, )
    minpos <- min(seq_along(xerror)[xerror == min(xerror)])
    if (minline) 
        abline(h = (xerror + xstd)[minpos], lty = lty, col = col)
    invisible()
}
#
library(rpart)
library(lasso2)
data(Prostate)
attach(Prostate)
treefit <- rpart(lpsa~.,data=Prostate,method="anova",control=list(minsplit=5,minbucket=3,cp=0.0))
#
# Plot of tree -- very cluttered and not included in book
#
plot(treefit,margin=.07,uniform=T,compress=T)
text(treefit,use.n=TRUE,cex=1)
# Details of fit
printcp(treefit)
# Now prune the tree
prunedtreefit <- prune(treefit,cp=treefit$cptable[which.min(treefit$cptable[,"xerror"]),"CP"])
#
# Fig 12.10
#
pdf("TreeProstateFig1.pdf",h=6,w=6) 
par(pty="s")
myplotcp(treefit,minline=F)
dev.off()
#
# Fig 12.11
#
pdf("TreeProstateFig2.pdf",h=6,w=6) 
par(pty="s")
plot(prunedtreefit,margin=.07,uniform=T)
text(prunedtreefit,use.n=TRUE,cex=1)
dev.off()
#
# Ethanol with regression trees
#
library(SemiPar)
data(ethanol)
attach(ethanol)
library(rpart)
TREEmod <- rpart(NOx~E+C,method="anova")
#
# Hierarchical tree etc initial plots: not in book
#
par(mfrow=c(2,2))
plot(TREEmod,margin=.1)
text(TREEmod)
plot(NOx~E)
lines(loess.smooth(y=NOx,x=E,span=.3))
# Lines are from tree model
abline(v=1.095,col="red",lty=2)
abline(v=0.796,col="red",lty=2)
abline(v=0.646,col="red",lty=2)
abline(v=1.023,col="red",lty=2)
plot(NOx~C)
lines(loess.smooth(y=NOx,x=C))
abline(v=8.25,col="red",lty=2)
# Cross-validation plot
plotcp(TREEmod)
#
# Now prune the tree
#
prunedtreefit <- prune(TREEmod,cp=TREEmod$cptable[which.min(TREEmod$cptable[,"xerror"]),"CP"])
#
Eval <- seq(min(E),max(E),length=100)
Cval <- seq(min(C),max(C),length=100)
fitf <- matrix(0,nrow=length(Eval),ncol=length(Cval))
#
# Create fitted values using a sledgehammer
#
for (i in 1:length(Eval)){
    for (j in 1:length(Cval)){
        if (Eval[i] >= 1.095) {fitf[i,j] <- 0.8753}
        else {
          if (Eval[i] < 0.796){ fitf[i,j]<-ifelse (Eval[i]<0.646,0.6402,1.651)}
          else {
             if (Eval[i] >= 1.023) {fitf[i,j] <- 2.16}
             else{ fitf[i,j] <- ifelse(Cval[j] > 8.25,3.023,3.544) }
       }}}
}
#
# Image plot
#
# Fig 12.12(a)
#
pdf("EthanolTree.pdf",h=5,w=5) 
par(mfrow=c(1,1))
image(Eval,Cval,fitf,xlab="E",ylab="C",col=gray((0:32)/32))
dev.off()
#
#
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Ethanol MARS
#
# First illustrate the use of the earth package
#
library(SemiPar)
data(ethanol)
attach(ethanol)
library(earth)
Emod <- earth(NOx~E+C,data=ethanol)
Eval <- seq(min(E),max(E),length=25)
Cval <- seq(min(C),max(C),length=25)
fitf <- matrix(0,nrow=length(Eval),ncol=length(Cval))
for (i in 1:length(Eval)){
    for (j in 1:length(Cval)){
    	fitf[i,j] <- predict(Emod,newdata=c(Eval[i],Cval[j]))
    }
}
#
# Fig 12.12(b)
#
pdf("ethMARSfig1.pdf",h=5,w=5) 
par(mfrow=c(1,1))
image(Eval,Cval,fitf,xlab="E",ylab="C",col=gray((0:32)/32))
dev.off()
# Perspective version: not in book
par(mfrow=c(1,1))
persp(Eval,Cval,fitf,theta=30,phi=30,ticktype="detailed")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Linear and quadratic discriminant analysis of BPD and birthweight data
# in Section 12.8.2
#
#----------------------------------------------------------
# Define expit functions to obtain fitted logistic linear regression on [0,1]
#
expit <- function(x,b0,b1){
      expit <- exp(b0+b1*x)/( 1+exp(b0+b1*x) )}
#----------------------------------------------------------
# Define expit functions to obtain fitted logistic quadratic regression on [0,1]
#
expit2 <- function(x,b0,b1,b2){
      expit <- exp(b0+b1*x+b2*x*x)/( 1+exp(b0+b1*x+b2*x*x) )}
bw <- read.table("birthweight.txt",header=T)
birthweight <- bw$birthweight
BPD <- bw$BPD
birthweight2 <- birthweight^2
#
# Linear and quadratic logistic regression
#
lrmod1 <- glm(BPD~birthweight,family=binomial)
lrmod2 <- glm(BPD~birthweight+birthweight2,family=binomial)
#
# Fig 12.13
#
pdf("bwlogitpredfig1.pdf",height=4,width=6)
plot(birthweight,BPD,pch="|",xlab="Birthweight (grams)")
x <- seq(min(birthweight),max(birthweight))
lines(x,expit(x,b0=lrmod1$coeff[1],b1=lrmod1$coeff[2]),lty=2)
lines(x,expit2(x,b0=lrmod2$coeff[1],b1=lrmod2$coeff[2],b2=lrmod2$coef[3]),lty=3)
abline(h=0.5,lty=1,col="grey")
abline(v=-lrmod1$coeff[1]/lrmod1$coeff[2],lty=2)
xquad1 <- (-lrmod2$coeff[2] + sqrt(lrmod2$coeff[2]**2-
          4*lrmod2$coeff[1]*lrmod2$coeff[3]))/(2*lrmod2$coeff[3])
xquad2 <- (-lrmod2$coeff[2] - sqrt(lrmod2$coeff[2]**2-
          4*lrmod2$coeff[1]*lrmod2$coeff[3]))/(2*lrmod2$coeff[3])
abline(v=xquad2,lty=3)
legend(x=1100,y=.95,legend=c("Linear Logistic","Quadratic Logistic"),lty=2:3,bty="n")
dev.off()
#
# Examination of normality within each disease classification
#
# Fig 12.14(a)
#
pdf("BPDqqBPD0.pdf",h=4,w=4)
qqnorm(birthweight[BPD==0],main="")
dev.off()
#
# Fig 12.14(b)
#
pdf("BPDqqBPD1.pdf",h=4,w=4)
qqnorm(birthweight[BPD==1],main="")
dev.off()
#
# Empirical estimates of prior probabilities of being in each group
#
pi1 <- sum(BPD)/length(BPD)
pi0 <- 1-pi1
#
# Empirical estimates of means and common variance (for LDA)
#
mu0 <- mean(birthweight[BPD==0])
mu1 <- mean(birthweight[BPD==1])
Sigma <- (sum((birthweight[BPD==0]-mu0)^2) + sum((birthweight[BPD==1]-mu1)^2))/(length(BPD)-2)
#
# 
LDAa0 <- -2*log(pi0) + mu0^2/Sigma
LDAa1 <- -2*log(pi1) + mu1^2/Sigma
LDAb0 <- -2*mu0/Sigma
LDAb1 <- -2*mu1/Sigma
xseq <- seq(min(birthweight),max(birthweight),1)
line0 <- LDAa0 + LDAb0*xseq
line1 <- LDAa1 + LDAb1*xseq
alpha0 <- log(pi1/pi0) + 0.5*mu0*mu0/Sigma - 0.5*mu1*mu1/Sigma
alpha1 <- mu1/Sigma - mu0/Sigma
cat("LDA boundary = ",-alpha0/alpha1,"\n") 
#
# Fig 12.15
#
pdf("BPDlindis.pdf",h=5,w=5)
plot(line0~xseq,type="n",ylab="-2 log Pr(Y=k | x)",xlab="Birthweight, x",
     ylim=c(min(line0,line1),max(line0,line1)))
lines(x=xseq,y=line0,lty=2)
lines(x=xseq,y=line1,lty=3)
abline(v=-alpha0/alpha1)
legend("topright",legend=c("k=0","k=1"),lty=2:3,bty="n")
dev.off()
#
# Variance estimates for QDA
#
Sigma0 <- sum((birthweight[BPD==0]-mu0)^2)/(length(birthweight[BPD==0])-1)
Sigma1 <- sum((birthweight[BPD==1]-mu1)^2)/(length(birthweight[BPD==1])-1)
# Quadratic discrimant
a <- (.5/Sigma0 - .5/Sigma1)
b <- (mu1/Sigma1 - mu0/Sigma0)
c <- log(pi1) - log(pi0) + .5*log(Sigma0) - .5*log(Sigma1) + .5*mu0^2/Sigma0 - .5*mu1^2/Sigma1
root1 <- .5*(-b+sqrt(b^2-4*a*c))/a
root2 <- .5*(-b-sqrt(b^2-4*a*c))/a
qdafun <- qda(BPD~birthweight) # The quick way of doing it!
#
# Fig 12.16
#
pdf("BPDquaddiscrim.pdf",h=4,w=6)
plot(birthweight,BPD,pch="|",xlab="Birthweight (grams)")
xseq <- seq(min(birthweight),max(birthweight),1)
scalenorm0 <- dnorm(xseq,mean=mu0,sd=sqrt(Sigma0))
max0 <- max(scalenorm0)
scalenorm1 <- dnorm(xseq,mean=mu1,sd=sqrt(Sigma1))
max1 <- max(scalenorm1)
points(xseq,pi0*scalenorm0/max0,type="l")
points(xseq,pi1*scalenorm1/max1,type="l")
abline(v=root2,lty=2)
dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# KDE and classification applied to BPD and birthweight in Section 12.8.3
#
library(sm)
myevalpts <- seq(min(birthweight),max(birthweight),1)
#
# Default smoothing parameter estimation is assuming normal distributions
#
sm0 <- sm.density(birthweight[BPD==0],display="none",eval.points=myevalpts)
sm1 <- sm.density(birthweight[BPD==1],display="none",eval.points=myevalpts)
sm0$h # Smoothing parameter for BPD=0
sm1$h # Smoothing parameter for BPD=1
#
# Fig 12.17(a)
#
pdf("KDEBPDfig1.pdf",h=5,w=5)
par(mfrow=c(1,1))
plot(birthweight,BPD,type="n",ylim=c(0,max(sm0$est,sm1$est)),
xlim=c(min(sm0$eval,sm1$eval),max(sm0$eval,sm1$eval)),ylab="",xlab="Birthweight, x")
points(sm0$eval,sm0$est,type="l")
points(sm1$eval,sm1$est,type="l",lty=2)
legend("topright",legend=c("k=0","k=1"),bty="n",lty=1:2,lwd=2)
dev.off()
pi1 <- sum(BPD)/length(BPD)
pi0 <- 1-pi1
#
yvals0 <- log(sm0$est/sm1$est)-log(pi0/pi1)
ind0 <- min(which(yvals0>0))
#
# Fig 12.17(b)
#
pdf("KDEBPDfig2.pdf",h=5,w=5)
plot(myevalpts,log(sm0$est/sm1$est)-log(pi0/pi1),type="l",xlab="Birthweight,x",
   ylab="log[Pr(Y=1|x)/Pr(Y=0|x)]")
abline(h=0,lty=2)
abline(v=myevalpts[ind0],lty=2)
dev.off()
#
# CV approach to smoothing parameter selection
#
h0hat <- hcv(birthweight[BPD==0],display="add") # CV estimate for BPD=0
h1hat <- hcv(birthweight[BPD==1],display="add") # CV estimate for BPD=1
sm0h <- sm.density(birthweight[BPD==0],h=h0hat,display="none",
    eval.points=myevalpts)
sm1h <- sm.density(birthweight[BPD==1],h=h1hat,display="none",
    eval.points=myevalpts)
#
# Fig 12.17(c)
#
pdf("KDEBPDfig3.pdf",h=5,w=5)
par(mfrow=c(1,1))
plot(birthweight,BPD,type="n",xlim=c(min(sm0h$eval,sm1h$eval),max(sm0h$eval,sm1h$eval)),
ylab="",xlab="Birthweight, x",ylim=c(0,max(sm0h$est,sm1h$est)))
points(sm0h$eval,sm0h$est,type="l")
points(sm1h$eval,sm1h$est,type="l",lty=2)
legend("topright",legend=c("k=0","k=1"),bty="n",lty=1:2,lwd=2)
dev.off()
#
# Fig 12.17(d)
#
pdf("KDEBPDfig4.pdf",h=5,w=5)
yvals1 <- log(sm0h$est/sm1h$est)-log(pi0/pi1)
ind1 <- min(which(yvals1>0))
plot(myevalpts,log(sm0h$est/sm1h$est)-log(pi0/pi1),type="l",xlab="Birthweight,x",
   ylab="log[Pr(Y=1|x)/Pr(Y=0|x)]")
abline(v=myevalpts[ind1],lty=2)
abline(h=0,lty=2)
dev.off()
#
# Now plug-in method of smoothing parameter estimation (Sheather-Jones)
#
sm0hsj <- sm.density(birthweight[BPD==0],method="sj",display="none",
    eval.points=myevalpts)
sm1hsj <- sm.density(birthweight[BPD==1],method="sj",display="none",
    eval.points=myevalpts)
sm0hsj$h # SJ estimation for BPD=0
sm1hsj$h # SJ estimation for BPD=1
#
# Fig 12.17(e)
#
pdf("KDEBPDfig5.pdf",h=5,w=5)
par(mfrow=c(1,1))
plot(birthweight,BPD,type="n",xlim=c(min(sm0hsj$eval,sm1hsj$eval),max(sm0hsj$eval,sm1hsj$eval)),
ylab="",xlab="Birthweight, x",ylim=c(0,max(sm0hsj$est,sm1hsj$est)))
points(sm0hsj$eval,sm0hsj$est,type="l")
points(sm1hsj$eval,sm1hsj$est,type="l",lty=2)
legend("topright",legend=c("k=0","k=1"),bty="n",lty=1:2,lwd=2)
dev.off()
#
# Fig 12.17(f)
#
pdf("KDEBPDfig6.pdf",h=5,w=5)
yvals2 <- log(sm0hsj$est/sm1hsj$est)-log(pi0/pi1)
ind2 <- min(which(yvals2>0))
myevalpts[ind2]
plot(myevalpts,log(sm0hsj$est/sm1hsj$est)-log(pi0/pi1),type="l",xlab="Birthweight,x",
   ylab="log[Pr(Y=1|x)/Pr(Y=0|x)]")
abline(v=myevalpts[ind2],lty=2)
abline(h=0,lty=2)
dev.off()
#
# Classification trees in Section 12.8.4
#
#
pseq <- seq(0.0001,.99999,.001)
misclass <- NULL
for (i in 1:length(pseq)){
    misclass[i] <- 1-max(pseq[i],1-pseq[i])}
#
# Cost complexity plot
#
# Fig 12.18
#
pdf("CostComplexBinary.pdf",h=5,w=5)
plot(misclass~pseq,xlab=expression(p[j]),type="l",ylab="Lack of Fit")
points(2*pseq*(1-pseq)~pseq,type="l",lty=2)
dev <- -pseq*log(pseq) - (1-pseq)*log(1-pseq)
dev <- 0.5*dev/max(dev)
points(dev~pseq,type="l",lty=3)
legend("bottom",legend=c("Misclassification Error","Gini Index","Deviance"),lty=c(1,2,3),bty="n")
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Motivation for bagging in 12.8.5: 3 random bootstrap samples -- will differ from
# book figures as different samples.
#
# First bootstrap sample
#
library(rpart)
library(lasso2)
data(Prostate)
attach(Prostate)
scramble <- sample(c(1:97),size=97,replace=T)
traini <- scramble[1:97]
traindat <- Prostate[traini,]
treefit1 <- rpart(lpsa~.,data=traindat,method="anova",
control=list(minsplit=15,minbucket=10,cp=0.0))
prunedtreefit1 <- prune(treefit1,cp=treefit$cptable[which.min(treefit$cptable[,"xerror"]),"CP"])
plotcp(treefit1,minline=F)
#
# Fig 12.19(a)
#
pdf("TreeProstateBAG1.pdf",h=6,w=6)
par(pty="s")
plot(prunedtreefit1,margin=.07,uniform=T)
text(prunedtreefit1,use.n=TRUE,cex=1)
dev.off()
# 
# Second bootstrap sample
#
scramble <- sample(c(1:97),size=97,replace=T)
traini <- scramble[1:97]
traindat2 <- Prostate[traini,]
treefit2 <- rpart(lpsa~.,data=traindat2,method="anova",control=list(minsplit=15,minbucket=10,cp=0.0))
prunedtreefit2 <- prune(treefit2,cp=treefit$cptable[which.min(treefit$cptable[,"xerror"]),"CP"])
#
# Fig 12.19(b)
#
pdf("TreeProstateBAG2.pdf",h=6,w=6)
par(pty="s")
plot(prunedtreefit2,margin=.07,uniform=T)
text(prunedtreefit2,use.n=TRUE,cex=1)
dev.off()
# 
# Third bootstrap sample
#
scramble <- sample(c(1:97),size=97,replace=T)
traini <- scramble[1:97]
traindat3 <- Prostate[traini,]
treefit3 <- rpart(lpsa~.,data=traindat3,method="anova",control=list(minsplit=15,minbucket=10,cp=0.0))
prunedtreefit3 <- prune(treefit3,cp=treefit$cptable[which.min(treefit$cptable[,"xerror"]),"CP"])
#
# Fig 12.19(bc
#
pdf("TreeProstateBAG3.pdf",h=6,w=6)
par(pty="s")
plot(prunedtreefit3,margin=.07,uniform=T)
text(prunedtreefit3,use.n=TRUE,cex=1)
dev.off()
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outcome after head injury in Section 12.8: Code written by Theresa Smith
#
library(MASS)
library(BMA)
library(rpart)
library(randomForest)
library(ipred)
#
##############Run once to set stuff up##########################
#
# Head injury data
#
yagg <- c(9,5,5,7,58,11,32,12,19,6,21,14,45,7,61,15,7,12,19,25,20,7,42,17)
zagg <- c(47,77,11,24,29,24,13,16,15,44,18,38,11,16,11,21,1,6,2,15,0,2,7,7)
haem <- (c(0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1))
pup <- (c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1))
coma <- (c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1))
agec <- as.factor(c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2))
#               
head <- data.frame(yagg,zagg,haem,pup,coma,agec)
#####un stack data######
####One row for each observation so that we can do train-test###
outcome <- factor(c(rep(1,sum(yagg)),rep(0,sum(zagg))),labels=c(0,1))
haem2 <- (c(rep(haem,yagg),rep(haem,zagg)))
pup2 <- (c(rep(pup,yagg),rep(pup,zagg)))
coma2 <- (c(rep(coma,yagg),rep(coma,zagg)))
age2 <- factor(c(rep(agec,yagg),rep(agec,zagg)), labels=c(0,1,2))
head2 <- data.frame(outcome,haem2,pup2,coma2,age2)
colnames(head2) <- c("outcome",colnames(head)[3:6])
#
# Create matrix of 0's and 1's to keep track of interactions and model elements#######
#
design.it<-function(k){
  design<-matrix(nrow=2^k,ncol=k)
  for(i in 1:k){
    design[,i]=rep(rep(c(0,1),each=2^(k-i)),2^(i-1))
  }
  return(design)
}
#
###is a subset of b? ##############
#
is.subset<-function(a,b,k){
  yes<-1
  for(i in 1:k)
    {
      yes=yes*(a[i]<=b[i])
    }
  return(yes)
}
#
#####is this a hierarchical  model? ########
#
is.hier<-function(int.mat,model,p,k){
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
##matrix to help identify all possible interactions####
#
interactions <- design.it(4)[-1,]
order <- apply(interactions,1,sum)
keep = order<3
interactions <- interactions[keep,] 
vars <- names(head)[-c(1,2)]
colnames(interactions) = vars
varints <- apply(interactions, 1, function(x){paste(colnames(interactions)[as.logical(x)], 
           collapse=":")})
#
##matrix to identify which interactions to include####
#
possible.model <- design.it(10)
rownames(possible.model) <- 1:1024
check.hier <- function(c){
  return(is.hier(interactions,c,10,4))
}
#
##keep models that are hierarchical
#
check.it <- apply(possible.model,1,check.hier)
keep <- which(check.it==1)
hier.models <- possible.model[keep,]
rownames(hier.models) <- 1:113
#
#####find aic for all heirachical models#########
#
get.aic<-function(model,var,hdata,bicflag){
  my.formula <- as.formula(paste(c("outcome ~ 1", var[as.logical(model)]), collapse=" + "))
  my.model <- glm(my.formula, family = binomial, data=hdata)
  return((1-bicflag)*extractAIC(my.model)[2]+bicflag*extractAIC(my.model,k=log(dim(hdata)[1]))[2])
}
#
######classification error########
#
class.error <- function(truth,pred){
Bsum00 <- sum((pred==0) & (truth==0)) 
Bsum01 <- sum((pred==1) & (truth==0)) 
Bsum10 <- sum((pred==0) & (truth==1)) 
Bsum11 <- sum((pred==1) & (truth==1)) 
Btesterr <- 100*(Bsum01+Bsum10)/length(truth)
return(Btesterr)
}

######needed for classification with logistic regression####
expit<-function(a){return(exp(a)/(1+exp(a)))}
###############################
prettytable<-function(visited,models,vars=varints,howmany=5){
	howmany=min(howmany,length(visited))
	topmod<-vector("list",length=howmany)
	for(i in 1:howmany){
	topmod[[i]]<-vars[as.logical(models[names(sort(visited))[i],])]
	}
	return(topmod)
}
#
###############End Run Once####################
#
B <- 100
null.err<-main.err<-sat.err<-aic.err<-bic.err<-aic.unres.err<-bic.unres.err<-tree.err<-bag.err<-rf.err<-rep(0,B)
aicmods<-bicmods<-uaicmods<-ubicmods<-rep(0,B)
for(i in 1:B){
#
######break up data 70% train 30% test########
 train<-sample(1:931,652)
 head.train<-head2[train,]
 head.test<-head2[-train,]
#####null model##########
nullmod<-glm(outcome~1, family = binomial, data=head.train)
nulp<-predict(nullmod,head.test[,-1])
nulp<-expit(nulp)
nulregp<-as.numeric(nulp>.5)
null.err[i]=class.error(head.test[,1],nulregp)
#######main only############
mainmod<-glm(outcome~., family = binomial, data=head.train)
mp<-predict(mainmod,head.test[,-1])
mp<-expit(mp)
mregp<-as.numeric(mp>.5)
main.err[i]=class.error(head.test[,1],mregp)
#####saturated#######



##calculate aic for each hierarchical model####
get.my.aic<-function(c){return(get.aic(c,varints,head.train,bicflag=0))}

aic.vals<-apply(hier.models,1,get.my.aic)

##which is the best?###
abest.model<-varints[as.logical(hier.models[which.min(aic.vals),])]
amy.formula<-as.formula(paste(c("outcome ~1", abest.model), collapse=" + "))
amy.model<-glm(amy.formula, family = binomial, data=head.train)
aicmods[i]=which.min(aic.vals)
#####test error for best model#########
alrp1<-predict(amy.model,head.test[,-1])
alrp<-expit(alrp1)
alogregp<-as.numeric(alrp>.5)
aic.err[i]=class.error(head.test[,1],alogregp)

#######Unrestricted Models############
aic.ur<-apply(possible.model,1,get.my.aic)
##which is the best?###
aubest.model<-varints[as.logical(possible.model[which.min(aic.ur),])]
aumy.formula<-as.formula(paste(c("outcome ~1", aubest.model), collapse=" + "))
aumy.model<-glm(aumy.formula, family = binomial, data=head.train)
uaicmods[i]=which.min(aic.ur)
#####test error for best model#########
aulrp1<-predict(aumy.model,head.test[,-1])
aulrp<-expit(aulrp1)
aulogregp<-as.numeric(aulrp>.5)
aic.unres.err[i]=class.error(head.test[,1],aulogregp)



##calculate aic for each hierarchical model####
get.my.bic<-function(c){return(get.aic(c,varints,head.train,bicflag=1))}

bic.vals<-apply(hier.models,1,get.my.bic)

##which is the best?###
best.model<-varints[as.logical(hier.models[which.min(bic.vals),])]
my.formula<-as.formula(paste(c("outcome ~1", best.model), collapse=" + "))
my.model<-glm(my.formula, family = binomial, data=head.train)

#####test error for best model#########
bicp<-predict(my.model,head.test[,-1])
bicp<-expit(bicp)
bicregp<-as.numeric(bicp>.5)
bic.err[i]=class.error(head.test[,1],bicregp)
bicmods[i]=which.min(bic.vals)


#######Unrestricted Models############
bic.ur<-apply(possible.model,1,get.my.bic)
##which is the best?###
ubest.model<-varints[as.logical(possible.model[which.min(bic.ur),])]
umy.formula<-as.formula(paste(c("outcome ~1", ubest.model), collapse=" + "))
umy.model<-glm(umy.formula, family = binomial, data=head.train)
ubicmods[i]=which.min(bic.ur)
#####test error for best model#########
ubicp<-predict(umy.model,head.test[,-1])
ubicp<-expit(ubicp)
ubicregp<-as.numeric(ubicp>.5)
bic.unres.err[i]=class.error(head.test[,1],ubicregp)

#####regression tree########
treefit <- rpart(outcome~haem+pup+coma+agec,data=head.train,method="class",control=list(minsplit=5,minbucket=3,cp=.001))  
prunedtreefit <- prune(treefit,cp=treefit$cptable[which.min(treefit$cptable[,"xerror"]),"CP"])
treep <- predict(prunedtreefit,head.test[,-1],type="class")
tree.err[i]=class.error(head.test[,1],treep)
#
######### Random Forests##########
#
RFmod <- randomForest(as.factor(outcome)~haem+pup+coma+agec,data=head.train,importance=TRUE, 
        type="classification")
RFp <- predict(RFmod,head.test[,-1],type="class")
rf.err[i] = class.error(head.test[,1],RFp)
#
######Bagging#########
#
bagmod <- bagging(outcome~.,data=head.train,method="class",coob=TRUE)
bagp <- predict(bagmod,head.test[,-1],type="class")
bag.err[i] = class.error(head.test[,1],bagp)
if (floor(i/5)==ceiling(i/5)) cat(i," iterations completed\n")
}
#
#
######average error all in one place: results are slightly different because seeds are different
#
Error <- rbind(null.err,main.err, aic.err,bic.err,aic.unres.err,bic.unres.err,tree.err,
         bag.err,rf.err)
apply(Error,1,mean)
sqrt(apply(Error,1,var))
#
aic.picks<-table(aicmods)
bic.picks<-table(bicmods)
uaic.picks<-table(uaicmods)
ubic.picks<-table(ubicmods)
####more like a list than a table########
prettytable(aic.picks,hier.models)
prettytable(bic.picks,hier.models)
prettytable(uaic.picks,possible.model)
prettytable(ubic.picks,possible.model)
#
# A couple of figures based on the last split: these will differ from the book versions
# because of different splits
#
# Fig 12.20
#
pdf("RFhead1.pdf")
varImpPlot(RFmod,main="",scale=F)
dev.off()
#
# Fig 12.21
#
pdf("RFhead2.pdf")
plot(RFmod$err[,1],typ="l",xlab="Trees",ylab="oob error")
dev.off()
