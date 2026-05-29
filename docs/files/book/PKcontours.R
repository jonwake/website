#
# R CODE FOR REPRODUCING CONTOUR FIGURE 6.10 IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 6.
# CODE WRITTEN BY RUTH SALWAY
##############################################
#
Plot.PKlik.comp<-function(data,range,fixed,n=10, log=F, lev=NULL,
    nlev=10,xlab=NULL,ylab=NULL,llmax=0)
{
# data is a data.frame with columns conc, dose and time
# range is  the range for each of the two parameters on the x and y axes
# fixed lists the values at which the other params are kept fixed
# both are lists with elements names V, ka and ke as required
#
# log=T plots on the log(V) scale for normal lik - this is the lik 
# actually maximised. Check all parameters are given
#
data$logconc=log(data$conc)
    pnames=c("V","ka","ke","sigma")
    fnames=names(fixed)
    rnames=names(range)
    if (any((c(fnames,rnames) %in% pnames)==F)) 
                stop(paste("fixed and range must contain ",
                paste(pnames,collapse=", ",sep="")))
    if (any(fnames %in% rnames)==T) 
       stop("Each parameter can appear in either fixed or range, not both.")
    if (length(fnames)!=2) stop("fixed must be a list with two parameters.")
    if (length(rnames)!=2) stop("range must be a list with two parameters.")
#
# Set up parameters
#
    if ("V" %in% fnames) {
        V=rep(fixed$V,n*n)
        if (log) V=log(V)
    }
    else {
        if (log) V=seq(log(range$V[1]), log(range$V[2]), len = n)
        else V=seq(range$V[1], range$V[2], len = n)
    }
    if ("ka" %in% fnames) {
        ka=rep(fixed$ka,n*n)
        if (log) ka=log(ka)
    }
    else {
        if (log) ka=seq(log(range$ka[1]), log(range$ka[2]), len = n)
        else ka=seq(range$ka[1], range$ka[2], len = n)
    }
    if ("ke" %in% fnames) {
        ke=rep(fixed$ke,n*n)
        if (log) ke=log(ke)
    }
    else {
        if (log) ke=seq(log(range$ke[1]), log(range$ke[2]), len = n)
        else ke=seq(range$ke[1], range$ke[2], len = n)
    }
    if ("sigma" %in% fnames) sigma=rep(fixed$sigma,n*n)
    else sigma=seq(range$sigma[1], range$sigma[2], len = n)
#
    eval(parse(text=paste(rnames[1],"=rep(",rnames[1],",n)",sep="")))
    eval(parse(text=paste(rnames[2],"=rep(",rnames[2],",rep(n,n))",sep="")))
#
# parameters should now be n*n vectors
#
    ll=rep(0,n*n)
    for (i in 1:(n*n)) {
        if (!log) {
            fit=log( data$dose*ka[i]/(V[i]*(ka[i]-ke[i])) * 
                  (exp(-ke[i]*data$time)-exp(-ka[i]*data$time)) )
            ll[i] = loglik.lnorm(data$conc,fit,sigma[i])
        }
        else {
            fit = log(data$dose) + ka[i] - V[i] - log(exp(ka[i])-exp(ke[i])) +
                    log( exp(-exp(ke[i])*data$time)-exp(-exp(ka[i])*data$time))
            ll[i] = loglik.norm(data$logconc,fit,sigma[i])
#if(ll[i] >llmax)cat("llmax lli V ka ke: ",llmax,ll[i],V[i],ka[i],ke[i],"\n")
        }
    }
ll <- ll - llmax
cat("max ll = ",max(ll),"\n")
    if (is.null(lev))  lev=pretty(range(ll, finite = TRUE),nlev)

    if (!log) {
if (is.null(xlab)) xlab=rnames[1]
if (is.null(ylab)) ylab=rnames[2]
eval(parse(text=paste("contour(seq(range$",rnames[1],"[1],range$",
      rnames[1],"[2],len=n) ,seq(range$",rnames[2],"[1],range$",
      rnames[2],"[2],len=n), matrix(ll,n,n), lev = lev,xlab='",xlab,"',
      ylab='",ylab,"')",sep="")))
#        title("Compartment")
    }
    else {
if (is.null(xlab)) xlab=paste("log(",rnames[1],")",sep="")
if (is.null(ylab)) ylab=paste("log(",rnames[2],")",sep="")
eval(parse(text=paste("contour(seq(log(range$",rnames[1],"[1]),
     log(range$",rnames[1],"[2]),len=n) ,seq(log(range$",rnames[2],"[1]),
     log(range$",rnames[2],"[2]),len=n), matrix(ll,n,n), lev = lev,xlab='",
     xlab,"',ylab='",ylab,"')",sep="")))
#        title("Compartment")
    }
}
#+++++++++++++++++++++++++++++++++++++
loglik.norm <- function(y,mu,sigma){
    n <- length(y)
    -n*log(sigma) - 0.5*n*log(2*pi) - 0.5 * sum((y-mu)^2)/(sigma^2)
}
#+++++++++++++++++++++++++++++++++++++
loglik.lnorm <- function(y,mu,sigma){
         n <- length(y)
         -n*log(sigma) - sum(log(y)) - 0.5*n*log(2*pi) - 
         0.5 * sum((log(y)-mu)^2)/(sigma^2)
    }
#+++++++++++++++++++++++++++++++++++++
loglik.gam <- function(y,alpha,beta){
        n <- length(y)
        # Y ~ Ga(a,b) with mean=a/b var=a/b^2
        -n*lgamma(alpha) + sum(alpha*log(beta) +(alpha-1)*log(y)-beta*y)
    } 
##############################################
# Code to produce plots
##############################################
#
# Plot likelihood contours for compartment model
#
y <- Theoph$conc[24:33];x <- Theoph$Time[24:33];dose <- Theoph$Dose[24] #Ind 3
z <- log(y)
pkdata <- data.frame(z=z,t=x,dose=dose)
modpk <- nls(z~log(dose)+theta1-theta0-log(exp(theta1)-exp(theta2))+
         log( exp(-exp(theta2)*t)-exp(-exp(theta1)*t) ),data=pkdata,
         start=list(theta0=log(1.14),theta1=log(0.565),theta2=log(0.109)))
#
# range over which to plot: nsd gives number of sds in each direction
#
nsd=3 # Change to 30 for wider range in book (Fig 6.10)
Vr=exp(coef(modpk)[1]+c(-1,1)*nsd*sqrt(vcov(modpk)[1,1]))
kar=exp(coef(modpk)[2]+c(-1,1)*nsd*sqrt(vcov(modpk)[2,2]))
ker=exp(coef(modpk)[3]+c(-1,1)*nsd*sqrt(vcov(modpk)[3,3]))
conc <- y
time <- x
dose <- dose
pkdata2 <- data.frame(conc,time,dose)
#These are the MLEs - your code will presumably vary
V.est=exp(coef(modpk)[1])
ka.est=exp(coef(modpk)[2])
ke.est=exp(coef(modpk)[3])
sig.est=0.06322
# fit at the max
fitmax = log(pkdata2$dose) + log(ka.est) - log(V.est) - log(ka.est-ke.est) +
         log( exp(-ke.est*pkdata2$time)-exp(-ka.est*pkdata2$time))
llmax = loglik.norm(log(pkdata2$conc),fitmax,sig.est)
#
# fix ke and sigma at estimated values, and plot V and ka on axes
#
levels30 <- c(-10000,-8000,-6000,-4000,-2000,-1000,-500,-250,-100,-50,-20)
pdf("logliklogVlogka3.pdf",height=4,width=4)
par(mfrow=c(1,1),pty="s")
#
# uncomment out the following if 30 sds wanted 
#
#Plot.PKlik.comp(data=pkdata2,range=list(V=Vr,ka=kar),fixed=list(ke=ke.est,sigma=sig.est),lev=levels30,n=50,log=T,xlab="",ylab="",llmax=llmax)
Plot.PKlik.comp(data=pkdata2,range=list(V=Vr,ka=kar),fixed=list(ke=ke.est,sigma=sig.est),n=50,log=T,xlab="",ylab="",llmax=llmax)
mtext(expression(paste("log ",V,,sep="")),side=1,line=3)
mtext(expression(paste("log ",k[a],sep="")),side=2,line=3)
points(coef(modpk)[1],coef(modpk)[2],pch=19)
dev.off()
pdf("logliklogVlogke3.pdf",height=4,width=4)
par(mfrow=c(1,1),pty="s")
#Plot.PKlik.comp(data=pkdata2,range=list(V=Vr,ke=ker),fixed=list(ka=ka.est,sigma=sig.est),lev=levels30,n=50,log=T,xlab="",ylab="",llmax=llmax)
Plot.PKlik.comp(data=pkdata2,range=list(V=Vr,ke=ker),fixed=list(ka=ka.est,sigma=sig.est),n=50,log=T,xlab="",ylab="",llmax=llmax)
mtext(expression(paste("log ",V,sep="")),side=1,line=3)
mtext(expression(paste("log ",k[a],sep="")),side=2,line=3)
points(coef(modpk)[1],coef(modpk)[3],pch=19)
dev.off()
pdf("logliklogkalogke3.pdf",height=4,width=4)
par(mfrow=c(1,1),pty="s")
#Plot.PKlik.comp(data=pkdata2,range=list(ka=kar,ke=ker),fixed=list(V=V.est,sigma=sig.est),lev=levels30,n=50,log=T,xlab="",ylab="",llmax=llmax)
Plot.PKlik.comp(data=pkdata2,range=list(ka=kar,ke=ker),fixed=list(V=V.est,sigma=sig.est),n=50,log=T,xlab="",ylab="",llmax=llmax)
mtext(expression(paste("log ",k[a],sep="")),side=1,line=3)
mtext(expression(paste("log ",k[e],sep="")),side=2,line=3)
points(coef(modpk)[2],coef(modpk)[3],pch=19)
dev.off()
#abline(h=coef(modpk)[3],lwd=2,lty=3)
stop()
#
# Plot likelihood contours for GLM model: these are not included in the
# book chapter
#
#range over which to plot: nsd gives number of sds in each direction
b0r=glm$parameters[1]+c(-1,1)*nsd*sqrt(glm$var[1,1])
b1r=glm$parameters[2]+c(-1,1)*nsd*sqrt(glm$var[2,2])
b1r[b1r>0]=0
b2r=glm$parameters[3]+c(-1,1)*nsd*sqrt(glm$var[3,3])
b2r[b2r>0]=0
#
# MLEs from fitted model
#
b2.est=glm$parameters[3]
phi.est=glm$phi
#
# fix beta2 and plot beta0 and beta1 on axes
#
Plot.PKlik.glm(data=pkdata,range=list(beta0=b0r,beta1=b1r),fixed=list(beta2=b2.est,phi=phi.est),lev=c(-10000,-8000,-6000,-4000,-2000,-1000,-500,-250,-100,-50),n=50,xlab="",ylab="")
mtext(expression(beta[0]),side=1,line=3)
mtext(expression(beta[1]),side=2,line=3)
points(glm$parameters[1],glm$parameters[2],pch=19)
#
# plot on beta0 scale
#
Plot.PKlik.glm<-function(data,range,fixed, n=10, lev=NULL,nlev=10,xlab=NULL,ylab=NULL)
{
#
# data is a data.frame with columns conc, dose and time
# range is  the range for each of the two parameters on the x and y axes
# fixed lists the values at which the other params are kept fixed
# both are lists with elements names beta0, beta1 and beta2 as required
#
data$invtime=1/data$time
#
# Check all parameters are given
#
    pnames=c("beta0","beta1","beta2","phi")
    fnames=names(fixed)
    rnames=names(range)
    if (any((c(fnames,rnames) %in% pnames)==F)) stop(paste("fixed and range must contain ",paste(pnames,collapse=", ",sep="")))
    if (any(fnames %in% rnames)==T) stop("Each parameter can appear in either fixed or range, not both.")
    if (length(fnames)!=2) stop("fixed must be a list with two parameters.")
    if (length(rnames)!=2) stop("range must be a list with two parameters.")
#
# Set up parameters
#
    if ("beta0" %in% fnames) beta0=rep(fixed$beta0,n*n)
        else beta0=seq(range$beta0[1], range$beta0[2], len = n)
    if ("beta1" %in% fnames) beta1=rep(fixed$beta1,n*n)
        else beta1=seq(range$beta1[1], range$beta1[2], len = n)
    if ("beta2" %in% fnames) beta2=rep(fixed$beta2,n*n)
        else beta2=seq(range$beta2[1], range$beta2[2], len = n)
    if ("phi" %in% fnames) phi=rep(fixed$phi,n*n)
        else phi=seq(range$phi[1], range$phi[2], len = n)
    eval(parse(text=paste(rnames[1],"=rep(",rnames[1],",n)",sep="")))
    eval(parse(text=paste(rnames[2],"=rep(",rnames[2],",rep(n,n))",sep="")))
#
# parameters should now be n*n vectors
#
    ll=rep(0,n*n)
    for (i in 1:(n*n)) {
        fit=data$dose*exp(beta0[i]+beta1[i]*data$time+beta2[i]*data$invtime)
        ll[i] <- loglik.gam(data$conc,1/phi[i],1/(fit*phi[i]))
    }
#    
    if (is.null(lev))  lev=pretty(range(ll, finite = TRUE),nlev) 
    if (is.null(xlab)) xlab=rnames[1]
    if (is.null(ylab)) ylab=rnames[2]
    eval(parse(text=paste("contour(seq(range$",rnames[1],"[1],range$",
         rnames[1],"[2],len=n) ,seq(range$",rnames[2],"[1],range$",
         rnames[2],"[2],len=n), matrix(ll,n,n), lev = lev,xlab='",
         xlab,"',ylab='",ylab,"')",sep="")))
#
    title("GLM")
}


