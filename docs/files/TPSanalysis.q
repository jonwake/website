##
source("tps.q")
source("infantsAgg.dat")
infants <- infantsAgg
#
# fitWL, fitPL, fitML gives the output
#
## Function to generate a random deviate from a multivariate hypergeometric distribution
##
rmhyper <- function(Mk, m)
{
  k <- length(Mk)
  M <- sum(Mk)
  mk <- rep(0, k)
	
  ##	
  if(m > M) mk <- Mk
  else
  {
    mk[1] <- rhyper(1, Mk[1], M - Mk[1], m)
    if(k > 2)
    {
      for(j in 2:(k-1))
      {
        mk[j] <- rhyper(1, Mk[j], M - sum(Mk[1:j]), m - sum(mk[1:(j-1)]))
      }
    }
    mk[k] <- m - sum(mk[1:(k-1)])
  }
  return(mk)
}


##
expit <- function(x) exp(x) / (1 + exp(x))


## "True" disease model based on the complete data
##
fit <- glm(cbind(Y, N-Y) ~ as.factor(region) + sex + race*lbw, data=infants, family=binomial)
eta <- predict(fit)


## Phase I stratification
##
## CHOOSE ONE OF THESE:
##
scheme <- 2
##
if(scheme == 0) infants$strata <- rep(1, nrow(infants))
if(scheme == 1) infants$strata <- infants$region
if(scheme == 2) infants$strata <- (10*infants$region) + infants$sex
if(scheme == 3) infants$strata <- (10*infants$region) + infants$race
if(scheme == 4) infants$strata <- (10*infants$region) + infants$lbw
if(scheme == 5) infants$strata <- (100*infants$region) + (10*infants$sex) + infants$race
if(scheme == 6) infants$strata <- (100*infants$region) + (10*infants$sex) + infants$lbw
if(scheme == 7) infants$strata <- (100*infants$region) + (10*infants$race) + infants$lbw


## Ensure levels of phase I stratification variable are consecutive (required by the tps() code)
##
lvls <- sort(unique(infants$strata))
for(i in 1:length(lvls)) infants$strata[infants$strata == lvls[i]] <- i



## Phase II sample sizes
##
n  <- 1000
m1 <- n / (2 * length(lvls))
m0 <- n / (2 * length(lvls))
if(round(m0) != m0)
{
  m0 <- round(m0)
  m1 <- (n/length(lvls)) - m0
}


##
infants$alive <- infants$N - infants$Y
infants$death <- infants$Y


## Phase One data
nn0     <- aggregate(infants$alive, list(infants$strata), FUN=sum)
nn0[,1] <- as.numeric(as.vector(nn0[,1]))
nn1     <- aggregate(infants$death, list(infants$strata), FUN=sum)
nn1[,1] <- as.numeric(as.vector(nn1[,1]))
  
##
okStrata <- intersect(nn0$Group[nn0$x > 0], nn1$Group[nn1$x > 0])
nn0      <- nn0[is.element(nn0$Group, okStrata),]
nn1      <- nn1[is.element(nn1$Group, okStrata),]
  
##
phaseTwo <- NULL
for(i in 1:length(okStrata))
{
  ##
  temp <- infants[(infants$strata == okStrata[i]),]

  ##
  N1 <- sum(temp$death)
  if(N1 <= m1)
    temp$alive <- rmhyper(temp$alive, m0 + (m1-N1))

  ##
  if(N1 > m1){
    temp$alive <- rmhyper(temp$alive, m0)
    temp$death <- rmhyper(temp$death, m1)
  }

  ##
  phaseTwo <- rbind(phaseTwo, temp)
}

if(length(okStrata) != length(lvls))
{
  okStrataSort <- sort(okStrata)
  for(i in 1:length(okStrataSort)) phaseTwo$strata[phaseTwo$strata == okStrataSort[i]] <- i
}

  
## Fit model
##
fitWL <- try(tps(cbind(death, alive) ~ factor(region) + sex + race*lbw, phaseTwo, nn0$x, nn1$x, phaseTwo$strata, method="WL"))
fitPL <- try(tps(cbind(death, alive) ~ factor(region) + sex + race*lbw, phaseTwo, nn0$x, nn1$x, phaseTwo$strata, method="PL"))
fitML <- try(tps(cbind(death, alive) ~ factor(region) + sex + race*lbw, phaseTwo, nn0$x, nn1$x, phaseTwo$strata, method="ML"))

