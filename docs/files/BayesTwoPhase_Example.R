#####################################################################################
# Example analysing a data set using Bayesian two-phase approach
# Simple random sampling at phase I
# Data should be laid out as in Table 1 of Ross and Wakefield (2012)
# Function "Bayes2phase.func" takes input arguments:
#   - Nvec: vector of length 2m_z of phase I sample sizes
#   - nvec: 2 x m_x x m_z matrix containing phase II data
#   - maxiter: maximum number of iterations MCMC should perform (not including iterations for tuning)
# Function returns a list containing beta parameter estimates and 95% credible intervals
# Requires "tps.q"
#####################################################################################

# Use population data from simulation study so we can compare results
# Randomly sampled from this population to obtain a phase II data set

library(mvtnorm)

source("BayesTwoPhase.R")

Ntruth<-matrix(c(3930,430,94,10,109,47,212,168), nrow=2, ncol=4, byrow=TRUE)
Nvec<-apply(Ntruth, 2, sum)

Y<-c(rep(0, Nvec[1]), rep(1, Nvec[2]), rep(0, Nvec[3]), rep(1, Nvec[4]))
X<-c(rep(0, Ntruth[1,1]), rep(1,Ntruth[2,1]), rep(0, Ntruth[1,2]), rep(1, Ntruth[2,2]), rep(0, Ntruth[1,3]), rep(1,Ntruth[2,3]), rep(0, Ntruth[1,4]), rep(1, Ntruth[2,4]))
Z<-c(rep(0, Nvec[1]+Nvec[2]), rep(1, Nvec[3]+Nvec[4]))

# true model
options(contrasts=c("contr.treatment", "contr.poly"))
fit.true<-glm(Y~X+Z, family=binomial)
fit.true.CI<-confint.default(fit.true)

nvec<-matrix(c(238,92,72,5,12,8,178,95), nrow=2, ncol=4)

results<-Bayes2phase.func(Nvec, nvec, 250000)
results

##############################################################
## To compare results to NPML method
ntot<-apply(nvec, 2, sum)

y<-c(rep(0, ntot[1]), rep(1, ntot[2]), rep(0, ntot[3]), rep(1, ntot[4]))
x<-c(rep(0, nvec[1,1]), rep(1,nvec[2,1]), rep(0, nvec[1,2]), rep(1, nvec[2,2]), rep(0, nvec[1,3]), rep(1,nvec[2,3]), rep(0, nvec[1,4]), rep(1, nvec[2,4]))
z<-c(rep(0, ntot[1]+ntot[2]), rep(1, ntot[3]+ntot[4]))

# fit NPML method to get starting values and proposal variance-covariance matrix
options(contrasts=c("contr.treatment", "contr.poly"))
gr<-z+1
fitML<-tps(y~factor(x)+factor(z), nn0=Nvec[c(1,3)], nn1=Nvec[c(2,4)], group=gr, method="ML")
fitML$coef
##############################################################