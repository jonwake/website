#
# R CODE FOR REPRODUCING THE PK FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 6.
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#=========================================================================
# Function to form product of lognormal priors for t_half, t_max, conc_max
# and the coefficient of variation.
#=========================================================================
#
prior.onecomp.open <- function(sumprior=0){
	print(paste("Median and 95% point for half-life = "))
	thalfm <- scan("",nmax=1)
	thalfu <- scan("",nmax=1)
	print(paste("Median and 95% point for tmax = "))
	tmaxm <- scan("",nmax=1)
	tmaxu <- scan("",nmax=1)
	print(paste("Median and 95% point for concmax = "))
	concmaxm <- scan("",nmax=1)
	concmaxu <- scan("",nmax=1)
	print(paste("Median and 95% point for CV = "))
	cvm <- scan("",nmax=1)
	cvu <- scan("",nmax=1)
	priormoms <- rep(0,8)
	priormoms[1] <- log(thalfm)
	priormoms[2] <- (log(thalfu)-priormoms[1])/1.645
	priormoms[3] <- log(tmaxm)
	priormoms[4] <- (log(tmaxu)-priormoms[3])/1.645
	priormoms[5] <- log(concmaxm)
	priormoms[6] <- (log(concmaxu)-priormoms[5])/1.645
	priormoms[7] <- log(cvm)
	priormoms[8] <- (log(cvu)-priormoms[7])/1.645
	if (sumprior == 1){
	   par(mfrow=c(2,2))
	   hist(rlnorm(100,priormoms[1],priormoms[2]),xlab="HALF-LIFE")
	   hist(rlnorm(100,priormoms[3],priormoms[4]),xlab="T MAX")
	   hist(rlnorm(100,priormoms[5],priormoms[6]),xlab="CONC MAX")
	   hist(rlnorm(100,priormoms[7],priormoms[8]),xlab="CV")
	}
	return(priormoms)
}
#=========================================================================
# Return the endpoints of a CI
#=========================================================================
#
confint2 <- function(esthat,sehat,coverage){
        zval <- qnorm(1-(1-coverage)/2)
        lo <- esthat - zval*sehat
        hi <- esthat + zval*sehat
        list(lo=lo,hi=hi)
}
#=========================================================================
# Normal log-likelihood
#=========================================================================
#
normloglik <- function(time,z,mu,sigma){
		normloglik <- -length(y)*log(sigma) - 0.5*length(y)*log(2*pi)
		for (i in 1:length(y)){
		  normloglik <- normloglik - 0.5*(z[i]-mu[i])^2/(sigma*sigma)
		}
		return(normloglik)
	}
#=========================================================================
# Lognormal log-likelihood
#=========================================================================
#
lognormloglik <- function(time,y,mu,sigma){
	lognormloglik <- -length(y)*log(sigma) - sum(log(y)) - 0.5*length(y)*log(2*pi)
	for (i in 1:length(y)){
	  lognormloglik <- lognormloglik - 0.5*(log(y[i])-mu[i])^2/(sigma*sigma)
	}
	return(lognormloglik)
}
#=========================================================================
# Gamma Log-likelihood Ga( alpha^{-1}, [alpha mu]^{-1} )
#=========================================================================
#
gamloglik <- function(y,mu,alpha){
	a <- 1/alpha
	gamloglik <- -length(y)*lgamma(a)
	for (i in 1:length(y)){
		 b <- 1/(alpha*mu[i])
		 gamloglik <- gamloglik + a*log(b) + (a-1)*log(y[i]) - b*y[i]
	}
	return(gamloglik)
}
#=========================================================================
# Function to return interval estimates for function so interest for the GLM.
# Needs the paramater estimates, beta, and the (asymptotic) var-cov matrix 
# for the beta estimates.
# Also alpha (the overdispersion parameter), and the asymptotic variance 
# of the estimate.
#=========================================================================
#
GLMderiveCI <- function(beta0,beta1,beta2,Vbeta,dose,phi,Vphi,CIprob){
#
#	delta method for thalf
#
  gfun <- log(log(2)) - log(-beta1)
  g1dhat <- 0
  g2dhat <- 1/beta1
  g3dhat <- 0
  varfun <- Vbeta[1,1]*g1dhat^2+2*Vbeta[1,2]*g1dhat*g2dhat+
          2*Vbeta[1,3]*g1dhat*g3dhat+Vbeta[2,2]*g2dhat^2+
          2*Vbeta[2,3]*g2dhat*g3dhat+Vbeta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  thalflo <- exp(out$lo)
  thalfhi <- exp(out$hi)
#
# 	  delta method for tmax
#
  gfun <- 0.5*log(-beta2) - 0.5*log(-beta1)
  g1dhat <- 0
  g2dhat <- -1/(2*beta1)
  g3dhat <- 1/(2*beta2)
  varfun <- Vbeta[1,1]*g1dhat^2+2*Vbeta[1,2]*g1dhat*g2dhat+
          2*Vbeta[1,3]*g1dhat*g3dhat+Vbeta[2,2]*g2dhat^2+
          2*Vbeta[2,3]*g2dhat*g3dhat+Vbeta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  tmaxlo <- exp(out$lo)
  tmaxhi <- exp(out$hi)
#
#	delta method for concmax
#
  gfun <- beta0 - 2*sqrt(beta1*beta2)
  g1dhat <- 1
  g2dhat <- sqrt(beta2/beta1)
  g3dhat <- sqrt(beta1/beta2)
  varfun <- Vbeta[1,1]*g1dhat^2+2*Vbeta[1,2]*g1dhat*g2dhat+
          2*Vbeta[1,3]*g1dhat*g3dhat+Vbeta[2,2]*g2dhat^2+
          2*Vbeta[2,3]*g2dhat*g3dhat+Vbeta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  concmaxlo <- exp(out$lo)
  concmaxhi <- exp(out$hi)
#
#	delta method for clearance
#
  temp <- 2*sqrt(beta1*beta2)
  gfun <- 0.5*log(beta1/beta2) -log(2) - beta0 -log(besselK(temp,1)) + log(dose)
  g1dhat <- -1
  g2dhat <- 1/beta1 + sqrt(beta2/beta1)*besselK(temp,0)/besselK(temp,1)
  g3dhat <- sqrt(beta1/beta2)*besselK(temp,0)/besselK(temp,1)
  varfun <- Vbeta[1,1]*g1dhat^2+2*Vbeta[1,2]*g1dhat*g2dhat+
          2*Vbeta[1,3]*g1dhat*g3dhat+Vbeta[2,2]*g2dhat^2+
          2*Vbeta[2,3]*g2dhat*g3dhat+Vbeta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  clearlo <- exp(out$lo)
  clearhi <- exp(out$hi)
#
#	delta method for CV = 1/sqrt(phi) 
#
  gfun <- -0.5*log(phi)
  gdhat <- -1/(2*phi)
  varfun <- Vphi*gdhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  cvlo <- 100*exp(out$lo)
  cvhi <- 100*exp(out$hi)
#
  list(thalflo=thalflo,thalfhi=thalfhi,tmaxlo=tmaxlo,tmaxhi=tmaxhi,concmaxlo=concmaxlo,concmaxhi=concmaxhi,clearlo=clearlo,clearhi=clearhi,cvlo=cvlo,cvhi=cvhi)
}
#=========================================================================
# Function to return interval estimates for function so interest for the GLM.
# Needs the parameter estimates, beta, and the (asymptotic) var-cov matrix 
# for the beta estimates.
#=========================================================================
#
COMPderiveCI <- function(theta1,theta2,theta3,Vtheta,dose,sigma,Vsigma,CIprob){
  Vd <- exp(theta1)
  ka <- exp(theta2)
  ke <- exp(theta3)
  xmax <- log(ka/ke)/(ka-ke)
#
#	delta method for max concentration
#
  concmax <- (dose*ka/(Vd*(ka-ke)))*(exp(-ke*xmax)-exp(-ka*xmax))
  concm2 <- (dose/Vd)*exp( (ke/(ka-ke))*log(ke/ka) )
  gfun <- log(concm2)
  temp1 <- exp(theta2-theta3)
  g1dhat <- -1
  g2dhat <- -1/(temp1-1) - (theta2-theta3)*temp1/(temp1-1)^2
  g3dhat <- 1/(temp1-1) - (theta2-theta3)*temp1/(temp1-1)^2
  varfun <- Vtheta[1,1]*g1dhat^2+2*Vtheta[1,2]*g1dhat*g2dhat+
            2*Vtheta[1,3]*g1dhat*g3dhat+Vtheta[2,2]*g2dhat^2+
            2*Vtheta[2,3]*g2dhat*g3dhat+Vtheta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  concmaxlo <- exp(out$lo);
  concmaxhi <- exp(out$hi)
#
#	delta method for clearance 
#
  gfun <- theta1+theta3
  g1dhat <- 1
  g2dhat <- 0
  g3dhat <- 1
  varfun <- Vtheta[1,1]*g1dhat^2+2*Vtheta[1,2]*g1dhat*g2dhat+
            2*Vtheta[1,3]*g1dhat*g3dhat+Vtheta[2,2]*g2dhat^2+
            2*Vtheta[2,3]*g2dhat*g3dhat+Vtheta[3,3]*g3dhat^2
  out <- confint2(gfun,sqrt(varfun),CIprob)
  clearlo <- exp(out$lo)
  clearhi <- exp(out$hi)
#
#	delta method for time to max 
#
  gfun <- log( xmax )
  g1dhat <- 0
  g2dhat <- -ka/(ka-ke) + 1/(theta2 - theta3)
  g3dhat <- ke/(ka-ke) - 1/(theta2 - theta3)
  varfun <- Vtheta[1,1]*g1dhat^2+2*Vtheta[1,2]*g1dhat*g2dhat+
            2*Vtheta[1,3]*g1dhat*g3dhat+Vtheta[2,2]*g2dhat^2+
            2*Vtheta[2,3]*g2dhat*g3dhat+Vtheta[3,3]*g3dhat^2
#cat("T TO MAX SE: ",sqrt(varfun),"\n")
  out <- confint2(gfun,sqrt(varfun),CIprob)  
  tmaxlo <- exp(out$lo)
  tmaxhi <- exp(out$hi)
#
#	delta method for x_1/2 
#
  gfun <- log(log(2)) - theta3
  g1dhat <- 0
  g2dhat <- 0
  g3dhat <- -1
  varfun <- Vtheta[1,1]*g1dhat^2+2*Vtheta[1,2]*g1dhat*g2dhat+
            2*Vtheta[1,3]*g1dhat*g3dhat+Vtheta[2,2]*g2dhat^2+
            2*Vtheta[2,3]*g2dhat*g3dhat+Vtheta[3,3]*g3dhat^2
cat("HALF-LIFE SE: ",sqrt(varfun),"\n")
  out <- confint2(gfun,sqrt(varfun),CIprob)  
  thalflo <- exp(out$lo)
  thalfhi <- exp(out$hi)
#
#	delta method for cv 
#
  gfun <- log(sigma)
  gdhat <- 1/sigma
  varfun <- Vsigma*gdhat^2
cat("CV SE: ",sqrt(varfun),"\n")
  out <- confint2(gfun,sqrt(varfun),CIprob)    
  cvlo <- 100*exp(out$lo)
  cvhi <- 100*exp(out$hi)
#
 list(thalflo=thalflo,thalfhi=thalfhi,tmaxlo=tmaxlo,tmaxhi=tmaxhi,concmaxlo=concmaxlo,concmaxhi=concmaxhi,clearlo=clearlo,clearhi=clearhi,cvlo=cvlo,cvhi=cvhi)
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#				MAIN PROGRAM     			   +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#	Bayesian and likelihood analyses of PK data using GLMs with
#	gamma errors or with one-compartment open model.
#       In the book, the data from individual 3 is analyzed.
#
library(MASS)
library(nlme)
data(Theoph)
#y <- Theoph$conc[2:11];x <- Theoph$Time[2:11];dose <- Theoph$Dose[2] #Ind 1
#y <- Theoph$conc[13:22];x <- Theoph$Time[13:22];dose <- Theoph$Dose[13] #Ind 2
y <- Theoph$conc[24:33];x <- Theoph$Time[24:33];dose <- Theoph$Dose[24] #Ind 3
#y <- Theoph$conc[24:26];x <- Theoph$Time[24:26];dose <- Theoph$Dose[24] #Ind 3
#y <- Theoph$conc[123:132];x<-Theoph$Time[123:132];dose<-Theoph$Dose[123]#Ind 12
## index to define analysis we want to perform
#index <- 1 # Lkhd GLM, as in Section 6.5 -- gives the results in rows 
           # 1 and 2 of Table 6.2 and the MLE and CI for clearance
#index <- 2 # Bayes GLM
#index <- 3 # Lkhd comp model
index <- 4 # Bayes comp model
#index <- 5 # Plot with combined model fit -- needs 1-4 to have been run first
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Lkhd GLM
#
if (index==1){
  z <- 1/x
  mod <-  glm(y~x+z,family=Gamma(link="log"))
  MLEgam <- gamma.shape(mod)
  summary(mod,dispersion=gamma.dispersion(mod))
  b0MLE <- mod$coeff[1]
  b1MLE <- mod$coeff[2]
  b2MLE <- mod$coeff[3]
  thalfMLE <- -log(2)/b1MLE
  tmaxMLE <- sqrt(b2MLE/b1MLE)
  concmaxMLE <- exp(b0MLE-2.0*sqrt(b1MLE*b2MLE))
  temp <- 2*sqrt(b1MLE*b2MLE)
  clearMLE <- dose*exp(0.5*log(b1MLE/b2MLE) -log(2) - b0MLE -
              log(besselK(temp,1)))
# Now let's get the MLE of alpha (and then convert to a CV)
# Note that the gamma,shape function returns the MLE of the reciprocal
# of the alpha parameter (in usual GLM notation).
  Vbeta <- vcov(mod)
  phiMLE <- gamma.shape(mod)$alpha
  Vphi <- gamma.shape(mod)$SE^2
  CIs <- GLMderiveCI(b0MLE,b1MLE,b2MLE,Vbeta,dose,phiMLE,Vphi,0.9)
#
# Results
#
  cat("GLM LKHD\n")
  cat("MLE and 90% CI: thalf: ",thalfMLE,CIs$thalflo,CIs$thalfhi,"\n")
  cat("MLE and 90% CI: tmax: ",tmaxMLE,CIs$tmaxlo,CIs$tmaxhi,"\n")
  cat("MLE and 90% CI: concmax: ",concmaxMLE,CIs$concmaxlo,CIs$concmaxhi,"\n")
  cat("MLE and 90% CI: Cl: ",clearMLE,CIs$clearlo,CIs$clearhi,"\n")
  cat("MLE and 90% CI: 100 x CV: ",100/sqrt(phiMLE),CIs$cvlo,CIs$cvhi,"\n")
# Sandwich estimation
  library(sandwich) 
  Vbetas <- sandwich(mod)
  CIsand <- GLMderiveCI(b0MLE,b1MLE,b2MLE,Vbetas,dose,phiMLE,Vphi,0.9)
# Give results
  cat("GLM SANDWICH ESTIMATION\n")
  cat("max and CIs: thalf: ",thalfMLE,CIsand$thalflo,CIsand$thalfhi,"\n")
  cat("max and CIs: tmax: ",tmaxMLE,CIsand$tmaxlo,CIsand$tmaxhi,"\n")
  cat("max and CIs: concmax: ",concmaxMLE,CIsand$concmaxlo,CIsand$concmaxhi,"\n")
  cat("max and CIs: Cl: ",clearMLE,CIsand$clearlo,CIsand$clearhi,"\n")
  cat("max and CIs: 100 x CV: ",100/sqrt(phiMLE),CIsand$cvlo,CIsand$cvhi,"\n")
#
# Now some diagnostics: Pearson residuals, as in Section 6.9
#
  vary <- fitted(mod)^2/sqrt(phiMLE)
  res <- (y-fitted(mod))/sqrt(vary)
#
# Fig 6.6(a)
#
  pdf("GLMassess1.pdf",height=4,width=4)
  plot( x,res,xlab="Time (hours)",ylab="Pearson residuals",
  ylim=c(-max(abs(res)),max(abs(res))) )
  abline(0,0,lty=2)
  dev.off()
#
# Fig 6.6(b)
#
   pdf("GLMassess2.pdf",height=4,width=4)
   plot(fitted(mod),abs(res),xlab="Fitted values",ylab="| Pearson residuals |")
   dev.off()
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bayes GLM: implementation is via a rejection algorithm that samples from
# the prior on the parameters of interest and then converts to the parameters
# that directly specify the likelihood.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
if (index==2){
#
#	First set up prior. 8 (12), 1.5 (3), 9 (12), 0.05 (0.1)
#       These values were obtained from the prior.onecomp.open function
#       with the above values being input (when asked for).
#
   moms <- c(2.0794415,0.2464833,0.4054651,0.4213661,2.1972246,0.1748827,
                -2.9957323,0.4213661)
#   sumprior <- 0
#  moms <- prior.onecomp.open(sumprior)
  priormoms <- moms
#
#	Maximize likelihood.
#
  x2 <- 1/x
  mod <-  glm(y~x+x2,family=Gamma(link="log"))
  fittedmax <- fitted(mod)
  phiMLE <- gamma.shape(mod)$alpha
  maxloglik <- gamloglik(y=y,mu=fittedmax,alpha=1/phiMLE)
#
  print(paste("No of MC samples = "))
  nsamp <- scan("",nmax=1)
  GLMsamples <- matrix(0,nrow=nsamp,ncol=4)
  naccept <- counter <- pyest <- 0
#
  while(naccept < nsamp){
	counter <- counter + 1
	thalf <- rlnorm(1,moms[1],moms[2])
	tmax <- rlnorm(1,moms[3],moms[4])
	concmax <- rlnorm(1,moms[5],moms[6])
	cv <- rlnorm(1,moms[7],moms[8])
	u <- runif(1)
	b1 <- -log(2)/thalf
	b2 <- b1*tmax^2
	b0 <- log(concmax) + 2*sqrt(b1*b2)
	alpha <- cv*cv
	fitted <- exp(b0+b1*x+b2/x)
	loglik <- gamloglik(y=y,mu=fitted,alpha=alpha)
 	if (is.na(loglik)==T){ 
	   counter <- counter-1
	   cat("WEIRDO POINT\n")
	}
	if (is.na(loglik)==F){ 
	   pyest <- pyest + exp(loglik)
	   test <- loglik - maxloglik
 #cat(maxloglik,loglik,"\n")
	   if (test > 0) print("MAX IS NOT MAX!")
 	   if (log(u) < test){
		naccept <- naccept + 1
		GLMsamples[naccept,1] <- thalf
		GLMsamples[naccept,2] <- tmax
		GLMsamples[naccept,3] <- concmax
		GLMsamples[naccept,4] <- cv
		if (ceiling(naccept/10)-floor(naccept/10)==0) 
                      cat("Accepted = ",naccept,"\n")
	}
       }
}
  cat("Counter = ",counter,"\n")
  paccept <- nsamp/counter
  cat("Acceptance prob = ",paccept,"\n")
  cat("p(y|GLM) from rejection = ",paccept*exp(maxloglik),"\n")
  cat("p(y|GLM) as average from prior = ",pyest/counter,"\n")
# 
  par(mfrow=c(1,1))
# 
#  Fig 6.4(a)
#
  pdf("GLMpost1.pdf",height=4,width=4)
  hist(GLMsamples[,1],xlab="Half-life",prob=TRUE,main="")
  curve(dlnorm(x,moms[1],moms[2]),lty=1,add=T)
  dev.off()
# 
#  Fig 6.4(b)
#
  pdf("GLMpost2.pdf",height=4,width=4)
  hist(GLMsamples[,2],xlab="Time to maximum",prob=TRUE,main="")
  curve(dlnorm(x,moms[3],moms[4]),lty=1,add=T)
  dev.off()
# 
#  Fig 6.4(c)
#
  pdf("GLMpost3.pdf",height=4,width=4)
  hist(GLMsamples[,3],xlab="Maximum concentration",prob=TRUE,main="")	
  curve(dlnorm(x,moms[5],moms[6]),lty=1,add=T)
  dev.off()
# 
#  Fig 6.4(d)
#
  pdf("GLMpost4.pdf",height=4,width=4)
  hist(100*GLMsamples[,4],xlab="Coefficient of variation",prob=TRUE,main="")
  curve(dlnorm(x,log(100)+moms[7],moms[8]),lty=1,add=T)
  dev.off()
#
# Get samples for other quantities, in particular the clearance.
#	
  b1 <- -log(2)/GLMsamples[,1]
  b2 <- b1*GLMsamples[,2]^2
  b0 <- log(GLMsamples[,3]) + 2.0*sqrt(b1*b2)
  b0med <- quantile(b0,0.5)
  b1med <- quantile(b1,0.5)
  b2med <- quantile(b2,0.5)
  temp <- 2.0*sqrt(b1*b2)
  cleargam <- dose*sqrt(b1/b2)/(2.0*exp(b0)*besselK(temp,1))
#
# Fig 6.5
#
  pdf("GLMpost6.pdf",height=4,width=4)
  hist(cleargam,xlab="Clearance",prob=TRUE,main="")
  dev.off()
#
# Bivariate posterior plots, if required
#
  par(mfrow=c(2,3))
  plot(GLMsamples[,1],GLMsamples[,2],xlab="HALF-LIFE",ylab="T MAX")
  plot(GLMsamples[,1],GLMsamples[,3],xlab="HALF-LIFE",ylab="CONC MAX")	
  plot(GLMsamples[,2],GLMsamples[,3],xlab="T MAX",ylab="CONC MAX")	
  plot(GLMsamples[,1],cleargam,xlab="HALF-LIFE",ylab="Cl")
  plot(GLMsamples[,2],cleargam,xlab="T MAX",ylab="Cl")
  plot(GLMsamples[,3],cleargam,xlab="C MAX",ylab="Cl")
#
# Table of results: To give results in Table 6.2
#
  cat("GLM Bayes\n")
  cat("Median and 90% CI: thalf: ",quantile(GLMsamples[,1],c(0.5,0.05,0.95)),
      "\n")
  cat("Median and 90% CI: tmax: ",quantile(GLMsamples[,2],c(0.5,0.05,0.95)),
      "\n")
  cat("Median and 90% CI: concmax: ",quantile(GLMsamples[,3],c(0.5,0.05,0.95)),
      "\n")
  cat("Median and 90% CI: Cl: ",quantile(cleargam,c(0.5,0.05,0.95)),"\n")
  cat("Median and 90% CI: 100 x CV: ",quantile(100*GLMsamples[,4],
      c(0.5,0.05,0.95)),"\n")	
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Lkhd Comp: basic model is log y = log f + eps with eps ~ N(0,sigma2)	
#
if (index == 3){
   z <- log(y)
   pkdata <- data.frame(z=z,t=x,dose=dose)
   modpk <- nls(z~log(dose)+theta1-theta0-log(exp(theta1)-exp(theta2))+
   log( exp(-exp(theta2)*t)-exp(-exp(theta1)*t) ),
   data=pkdata,start=list(theta0=log(1.14),theta1=log(0.565),theta2=log(0.109)))
   library(sandwich)
   cat("Model-based ses: ",sqrt(diag(vcov(modpk))),"\n")
   cat("Sandwich ses: ",sqrt(diag(sandwich(modpk))),"\n")
#
   theta1 <- coef(modpk)[1]
   theta2 <- coef(modpk)[2]
   theta3 <- coef(modpk)[3]
   Vd <- exp(theta1)
   ka <- exp(theta2)
   ke <- exp(theta3)
# Summarize results
  thalfMLE <- log(2)/ke
  tmaxMLE <- log(ka/ke)/(ka-ke)
  concmaxMLE <- (dose*ka/(Vd*(ka-ke)))*(exp(-ke*tmaxMLE)-exp(-ka*tmaxMLE))
  clearMLE <- Vd*ke
  Vtheta <- vcov(modpk)
  sigma <- sqrt(sum((z-fitted(modpk))^2)/(length(z)-3))
  Vsigma <- sigma*sigma/(2*length(z))
  CIs <- COMPderiveCI(theta1,theta2,theta3,Vtheta,dose,sigma,Vsigma,0.9)
#
# Give results: appear in rows 3 and 4 of Table 6.2
#
  cat("COMPARTMENTAL MODEL LKHD\n")
  cat("max and CIs: thalf: ",thalfMLE,CIs$thalflo,CIs$thalfhi,"\n")
  cat("max and CIs: tmax: ",tmaxMLE,CIs$tmaxlo,CIs$tmaxhi,"\n")
  cat("max and CIs: concmax: ",concmaxMLE,CIs$concmaxlo,CIs$concmaxhi,"\n")
  cat("max and CIs: Cl: ",clearMLE,CIs$clearlo,CIs$clearhi,"\n")
  cat("max and CIs: 100 x CV: ",100*sigma,CIs$cvlo,CIs$cvhi,"\n")
# Sandwich results
  Vthetas <- sandwich(modpk)
  CIsand <- COMPderiveCI(theta1,theta2,theta3,Vthetas,dose,sigma,Vsigma,0.9)
# Give results
  cat("COMPARTMENTAL MODEL SANDWICH ESTIMATION\n")
  cat("max and CIs: thalf: ",thalfMLE,CIsand$thalflo,CIsand$thalfhi,"\n")
  cat("max and CIs: tmax: ",tmaxMLE,CIsand$tmaxlo,CIsand$tmaxhi,"\n")
  cat("max and CIs: concmax: ",concmaxMLE,CIsand$concmaxlo,CIsand$concmaxhi,"
      \n")
  cat("max and CIs: Cl: ",clearMLE,CIsand$clearlo,CIsand$clearhi,"\n")
  cat("max and CIs: 100 x CV: ",100*sigma,CIsand$cvlo,CIsand$cvhi,"\n")
# Now some diagnostics
# Pearson residuals
# Note: these are on a log scale, while the GLM residuals are on the original
# scale.
  vary <- sigma^2 # log scale!
  fittedpk <- fitted(modpk)
  respk <- (z-fittedpk)/sqrt(vary)
#
# Fig 6.6(c)
#
  pdf("COMPassess1.pdf",height=4,width=4)
  plot(x,respk,xlab="Time (hours)",ylab="Pearson residuals",
       ylim=c(-max(abs(respk)),max(abs(respk))) )
  abline(0,0,lty=2)
  dev.off()
#
# Fig 6.6(d)
#
  pdf("COMPassess2.pdf",height=4,width=4)
  plot(fittedpk,abs(respk),xlab="Fitted values",ylab="| Pearson residuals |")
  dev.off()
#
# Fig 6.15
#
  pdf("COMPqqnorm.pdf",height=4,width=4)
  par(mfrow=c(1,1),pty="s")
  qqnorm(respk,xlab="Expected",ylab="Observed",xlim=c(-1.5,1.5),
         ylim=c(-1.5,1.5),main="")
  abline(0,1)
  dev.off()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Alternative models:
# Now let's do gls and maximum likelihood with a constant 
# cv variance model on the original scale
# These results appear in Table 6.3
#
  library(nlme)
  glsmod <- gnls(y~exp(log(dose)+theta1-theta0-log(exp(theta1)-exp(theta2))+
        log( exp(-exp(theta2)*t)-exp(-exp(theta1)*t) )),weights=varPower(fixed=1),
        data=pkdata,start=list(theta0=log(1.14),theta1=log(0.565),theta2=log(0.109)))
# Now MLE
  fn <- function(theta,y,x,dose){
     mu <- exp(log(dose)+theta[2]-theta[1]-log(exp(theta[2])-exp(theta[3]))+
           log( exp(-exp(theta[3])*x)-exp(-exp(theta[2])*x) ))
     sigma2 <- exp(theta[4])
     varf <- sigma2*mu
     ll <- 0
     for (i in 1:length(y)){
     	 ll <- ll - 0.5*log(varf[i]) - 0.5*(y[i]-mu[i])^2/varf[i]
    }
     -ll
   }
  fn1 <- function(theta,y,x,dose){
     mu <- exp(log(dose)+theta[2]-theta[1]-log(exp(theta[2])-exp(theta[3]))+
           log( exp(-exp(theta[3])*x)-exp(-exp(theta[2])*x) ))
     sigma2 <- exp(theta[4])
     varf <- sigma2*sqrt(mu)
     ll <- 0
     for (i in 1:length(y)){
     	 ll <- ll - 0.5*log(varf[i]) - 0.5*(y[i]-mu[i])^2/varf[i]
    }
     -ll
   }
#   out <- nlm(fn,c(coef(glsmod),log(0.08)),y=y,x=x,dose=dose,hessian=TRUE,print.level=2)
   out <- optim(c(coef(glsmod),log(0.08)),fn, y=y,x=x,dose=dose,hessian = TRUE)
#
# var-cov matrix
#
   vout <- solve(out$hessian)
   ses <- sqrt(diag(vout))
   cat("\nPARAMETER ESTIMATES\n")
   cat("MLE LOG: ",coef(modpk),"\n")
   cat("GLS ORIGINAL: ",coef(glsmod),"\n")
   cat("MLE ORIGINAL: ",out$par[1:3],"\n")
   cat("\nSTANDARD ERRORS\n")
   cat("MLE LOG: ",sqrt(diag(vcov(modpk))),"\n")
   cat("GLS ORIGINAL: ",sqrt(diag(vcov(glsmod))),"\n")
   cat("MLE ORIGINAL: ",sqrt(diag(solve(out$hessian)))[1:3],"\n")
#
# Fitted curves
#
   plot(x,y,ylim=c(0,10))
   xseq <- seq(0,max(x),.1)
   Vd1 <- exp(coef(glsmod)[1]); ka1 <- exp(coef(glsmod)[2]); ke1 <- exp(coef(glsmod)[3])
   Vd2 <- exp(out$par[1]); ka2 <- exp(out$par[2]); ke2 <- exp(out$par[3])
   fseq <- (dose*ka)/(Vd*(ka-ke))*(exp(-ke*xseq)-exp(-ka*xseq))
   fseq1 <- (dose*ka1)/(Vd1*(ka1-ke1))*(exp(-ke1*xseq)-exp(-ka1*xseq))
   fseq2 <- (dose*ka2)/(Vd2*(ka2-ke2))*(exp(-ke2*xseq)-exp(-ka2*xseq))
   lines(xseq,fseq,lty=1)
   lines(xseq,fseq1,lty=2)
   lines(xseq,fseq2,lty=3)
   legend("topright",legend=c("MLE LOG SCALE","GLS ORIGINAL","MLE ORIGINAL"),
    lty=c(1:3),bty="n")
   }
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#	Compartmental model with Bayes rejection algorithm
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
if (index == 4){
   z <- log(y)
   pkdata <- data.frame(z=z,t=x,dose=dose)
   modpk <- gnls(z~log(dose)+theta1-theta0-log(exp(theta1)-exp(theta2))+
   log( exp(-exp(theta2)*t)-exp(-exp(theta1)*t) ),
   data=pkdata,start=list(theta0=log(1.14),theta1=log(0.565),theta2=log(0.109)))
#
   VdMLE <- exp(coef(modpk)[1])
   kaMLE <- exp(coef(modpk)[2])
   keMLE <- exp(coef(modpk)[3])
   fitted <- fitted(modpk)
   sigmaMLE <- sqrt(sum((z-fitted)^2)/length(z))
#maxnorm <- normloglik(time=x,z=z,mu=fitted,sigma=sigmaMLE)
   maxlognorm <- lognormloglik(time=x,y=y,mu=fitted,sigma=sigmaMLE)
#if (length(y)>3) maxlognorm <- lognormloglik(time=x,y=y,mu=fitted,sigma=sigmaMLE)
#if (length(y)<=3) maxlognorm <- 5
#
# Set up prior.
#	First set up prior. 8 (12), 1.5 (3), 9 (12), 0.05 (0.1)
   moms <- c(2.0794415,0.2464833,0.4054651,0.4213661,2.1972246,0.1748827,
             -2.9957323,0.4213661)
   priormoms <- moms
#
# Start simulation -- rejection algorithm sampling from the prior
#
   print(paste("No of MC samples for normal = "))
   nsamp <- scan("",nmax=1)
   COMPsamples <- matrix(0,nrow=nsamp,ncol=7)
   naccept <- counter <- pyest <- 0
#
   f <- function(theta,tmax,ke) { tmax - (log(exp(theta)+ke)-log(ke))/exp(theta) }
   while(naccept < nsamp){
	counter <- counter + 1
# Sample from prior for derived parameters
	thalf <- rlnorm(1,moms[1],moms[2])
	tmax <- rlnorm(1,moms[3],moms[4])
	concmax <- rlnorm(1,moms[5],moms[6])
	cv <- rlnorm(1,moms[7],moms[8])
	u <- runif(1)
	ke <- log(2)/thalf
# To solve for ka we need to find the root of an equation
#	f <- function(x,ke,tmax){ tmax - log(x/ke)/(x-ke) }
#	ka <- uniroot(f,interval=c(ke,10),ke=ke,tmax=tmax)
if (kaMLE-ke > 0) mid <- log(kaMLE-ke) else mid <- 0
   	flo <- f(mid-5,ke=ke,tmax=tmax)
        fhi <- f(mid+5,ke=ke,tmax=tmax)
# This next check is to see if the function we solve has the same sign at 
# our initial values - 
# this needs some work - if not many points being lost then OK...
  if (flo*fhi > 0) {
     cat("TROUBLE WITH ROOT-FINDING ",thalf,tmax,concmax,"\n")
     test <- -10^10
  }
  else {	
     ka <- ke + exp(uniroot(f,interval=c(mid-5,mid+5),ke=ke,tmax=tmax )$root)
     V <- (dose*ka/(concmax*(ka-ke)))*(exp(-ke*tmax)-exp(-ka*tmax))
     sigma <- cv
     fitted <- log((dose*ka/(V*(ka-ke)))*(exp(-ke*x)-exp(-ka*x)))
#	test <- normloglik(time=x,z=z,mu=fitted,sigma=sigma) - maxnorm
     loglik <- lognormloglik(time=x,y=y,mu=fitted,sigma=sigma)
     test <- loglik - maxlognorm
     pyest <- pyest + exp(loglik)
     if (test > 0) {
     	print("MAX IS NOT MAX!")
        cat(lognormloglik(time=x,y=y,mu=fitted,sigma=sigma),maxlognorm,V,ka,ke,sigma,"\n")}
      }
      if (log(u) < test){
	naccept <- naccept + 1
	COMPsamples[naccept,1] <- thalf
	COMPsamples[naccept,2] <- tmax
	COMPsamples[naccept,3] <- concmax
	COMPsamples[naccept,4] <- cv
	COMPsamples[naccept,5] <- V
	COMPsamples[naccept,6] <- ka
	COMPsamples[naccept,7] <- ke
#	cat("Point accepted: ",V,ka,ke,V*ke*1000,cv,test,naccept,"\n")
        if (ceiling(naccept/10)-floor(naccept/10)==0) cat("Accepted = ",naccept,"\n")
	    }
}
cat("Counter = ",counter,"\n")
paccept <- nsamp/counter
cat("Acceptance prob = ",paccept,"\n")
# Estimates of normalizing constant
cat("p(y|comp) from rejection = ",paccept*exp(maxlognorm),"\n")
cat("p(y|comp) as average from prior = ",pyest/counter,"\n")
#	
  par(mfrow=c(1,1))
  pdf("COMPpost1.pdf",height=4,width=4)
  hist(COMPsamples[,1],xlab="Half-life",prob=TRUE,main="",ylim=c(0,.22))
  curve(dlnorm(x,moms[1],moms[2]),lty=1,add=T)
  dev.off()
  pdf("COMPpost2.pdf",height=4,width=4)
  hist(COMPsamples[,2],xlab="Time to maximum",prob=TRUE,main="")
  curve(dlnorm(x,moms[3],moms[4]),lty=1,add=T)
  dev.off()
  pdf("COMPpost3.pdf",height=4,width=4)
  hist(COMPsamples[,3],xlab="Maximum concentration",prob=TRUE,main="")	
  curve(dlnorm(x,moms[5],moms[6]),lty=1,add=T)
  dev.off()
  pdf("COMPpost4.pdf",height=4,width=4)
  hist(100*COMPsamples[,4],xlab="Coefficient of variation",prob=TRUE,main="")
  curve(dlnorm(x,log(100)+moms[7],moms[8]),lty=1,add=T)
  dev.off()
#  
# Fig 6.13
# 
  pdf("COMPpost5.pdf",height=4,width=4)
  plot(COMPsamples[,6],COMPsamples[,7],xlab=expression(k[a]),ylab=expression(k[e]))
  image(kde2d(COMPsamples[,6],COMPsamples[,7],n=50),col=grey(seq(1,0,length=50)),
        xlab=expression(k[a]),ylab=expression(k[e]),main="")#,xlim=c(1.5,3));
  box()
  clearsamps <- V*ke
  dev.off()
  Vdmed <- quantile(COMPsamples[,5],0.5)
  kamed <- quantile(COMPsamples[,6],0.5)
  kemed <- quantile(COMPsamples[,7],0.5)
#
# Results that appear in the last line of Table 6.2
#
  cat("GLM Bayes\n")
  cat("Median and 90% CI: thalf: ",quantile(COMPsamples[,1],c(0.5,0.05,0.95)),"\n")
  cat("MLE and 90% CI: tmax: ",quantile(COMPsamples[,2],c(0.5,0.05,0.95)),"\n")
  cat("MLE and 90% CI: concmax: ",quantile(COMPsamples[,3],c(0.5,0.05,0.95)),"\n")
  cat("MLE and 90% CI: Cl: ",quantile(clearsamps,c(0.5,0.05,0.95)),"\n")
  cat("MLE and 90% CI: 100 x CV: ",quantile(100*COMPsamples[,4],c(0.5,0.05,0.95)),"\n")	
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
if (index == 5){
#
# Fig 6.1
#
   pdf("pkfits-theoph.pdf",height=6,width=6)
   par(mfrow=c(1,1))
   xseq <- seq(0.01,max(x),.01)
   glmMLE <- exp(b0MLE+b1MLE*xseq+b2MLE/xseq)
   glmBAYES <- exp(b0med+b1med*xseq+b2med/xseq)
   compMLE <- (dose*ka)/(Vd*(ka-ke))*(exp(-ke*xseq)-exp(-ka*xseq))
   compBAYES <- (dose*kamed)/(Vdmed*(kamed-kemed))*
                (exp(-kemed*xseq)-exp(-kamed*xseq))
   plot(x,y,xlab="Time (hours)",ylab="Concentration (mg/liter)",ylim=c(0,max(glmMLE,glmBAYES,compMLE,compBAYES,y)))
   lines(xseq,glmMLE,lty=2)
   lines(xseq,glmBAYES,lty=3)
   lines(xseq,compMLE,lty=4)
   lines(xseq,compBAYES,lty=5)
   legend("topright",legend=c("GLM MLE","GLM BAYES","Nonlinear MLE","Nonlinear BAYES"),lty=c(2:5),bty="n")
   dev.off()
}