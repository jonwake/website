#####################################################################################
# Example analysing a data set using Bayesian two-phase approach
# Case-control sampling at Phase I
# Data should be laid out as in Table 1 of Ross and Wakefield (2012)
# Function "Bayes2phaseCC.func" takes input arguments:
#   - Nvec: vector of length 2m_z of phase I sample sizes
#   - nvec: 2 x m_x x m_z matrix containing phase II data
#   - maxiter: maximum number of iterations MCMC should perform (not including iterations for tuning)
# Function returns a list containing beta parameter estimates and 95% credible intervals
# Requires "tps.q"
#####################################################################################

# Use population data from case-control example from paper using sex as confounder (collapse across age groups)
# Phase II data set includes smoking status (x)

library(mvtnorm)

source("BayesTwoPhase_CaseControl.R")

Nvec<-c(6310, 281, 5576, 528)
nvec<-matrix(c(1053,34,1319,30, 91,19,645,69, 77,56,137,50, 88,125,350,280), nrow=4, ncol=4)

results<-Bayes2phaseCC.func(Nvec, nvec, 250000)
results

##############################################################
## To compare results to NPML method
ntot<-apply(nvec, 2, sum)

y<-c(rep(0, ntot[1]), rep(1, ntot[2]), rep(0, ntot[3]), rep(1, ntot[4]))
x<-c(rep(0, nvec[1,1]), rep(1,nvec[2,1]), rep(2,nvec[3,1]), rep(3,nvec[4,1]), rep(0, nvec[1,2]), rep(1, nvec[2,2]), rep(2, nvec[3,2]), rep(3, nvec[4,2]), rep(0, nvec[1,3]), rep(1,nvec[2,3]), rep(2,nvec[3,3]), rep(3,nvec[4,3]), rep(0, nvec[1,4]), rep(1, nvec[2,4]), rep(2, nvec[3,4]), rep(3, nvec[4,4]))
z<-c(rep(0, ntot[1]+ntot[2]), rep(1, ntot[3]+ntot[4]))

# fit NPML method to get starting values and proposal variance-covariance matrix
options(contrasts=c("contr.treatment", "contr.poly"))
gr<-z+1
fitML<-tps(y~factor(x)+factor(z), nn0=Nvec[c(1,3)], nn1=Nvec[c(2,4)], group=gr, method="ML")
fitML$coef
##############################################################
