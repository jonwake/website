##Set your working directory to one you would like to use.
##R will read files from and write files to this directory only.
#setwd("~")

# Check out 
# http://faculty.washington.edu/jonno/SISMIDmaterial
# for files you need

##Load appropriate libraries.
##You may need to download these
##using the install.packages() function

library(geostatsp)
library(geoRglm)
library(INLA)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(dplyr)
library(ggplot2)
library(SpatialEpi)

##Functions we will use throughout.
logit <- function(x){log(x/(1-x))}
expit <- function(x) {exp(x)/(1+exp(x))}

##Set seed.
set.seed(123456)

##Load GPS locations at which to simulate data.
##Load Kenya SpatialPolygons at different admin levels.
##KEN_adm0 and KEN_adm1 available from GADM
kenya.coord <- readRDS(gzcon(url('http://faculty.washington.edu/jonno/SISMIDmaterial/kenya.coord.sim.rds')))
adm0 <- readRDS(gzcon(url('http://faculty.washington.edu/jonno/SISMIDmaterial/KEN_adm0.rds')))
adm1 <- readRDS(gzcon(url('http://faculty.washington.edu/jonno/SISMIDmaterial/KEN_adm1.rds')))

##Get a boundary for Kenya.
points.kenya <- spsample(adm0,1000,"regular")
boundary.kenya <- inla.nonconvex.hull(points.kenya@coords,convex=0.1,
                                      resolution = c(120,120))


##Plot admin 1 polygons with their indices.
plot(adm1, axes=TRUE)
invisible(polygonsLabel(adm1, labels= adm1@data$ID_1, 
                        method = "centroid", cex = 0.75))

##Build and plot a mesh.
mesh.true <- inla.mesh.2d(boundary = boundary.kenya, 
                          offset = c(1,5),max.edge = c(0.1,100))
plot(mesh.true,main="")
plot(points.kenya,col="red",add=T)


##Make an spde object
spde <- inla.spde2.matern(mesh.true)

##Simulate some data from the field,
##assuming a national prevalence mean level of 7%
##(mu = logit(.07) ).
mu <- logit(0.07)
Q.true <- inla.spde.precision(spde,c(-1/2,1/2))
field.true <- as.numeric(inla.qsample(1,Q.true, seed=5L)) 
mean <- rep(mu,length(field.true))

##Plot true prevalence field.
##inla.mesh.projector() inla.mesh.project() reparamaterize information.
proj <- inla.mesh.projector(mesh.true,dims=c(300,300))
field.plot <- inla.mesh.project(proj,field=field.true,dims=c(400,400))

x <- proj$x
y <- proj$y
true.prev.all <- (expand.grid(x,y,KEEP.OUT.ATTRS = F))
names(true.prev.all) = c("x","y")
true.prev.plotdata <- true.prev.all %>%mutate(field = as.numeric(field.plot), 
                                              prev = exp(mu + field)/(1+exp(mu + field)) )
coordinates(true.prev.plotdata)=~x+y ##Make data.frame a SpatialPointsDataFrame

##Make sure plot data and SpatialPolygon have same CRS proj4string.
true.prev.plotdata@proj4string = adm0@proj4string 

##Subset SpatialPointsDataFrame by the SpatialPolygon of interest.
true.prev.plot.final <- true.prev.plotdata[adm0,]
true.prev.plot.final <- as.data.frame(true.prev.plot.final)

##Use ggplot2 to make nice plot of true underlying prevalence surface.
ggplot(true.prev.plot.final, aes(x,y)) + geom_tile(aes(fill = prev)) +
  scale_fill_gradient(limits = c(0,0.20))


##Now simulate point level data.
##Two-stage cluster design: 1) Clusters 2) Househoulds.

A <- inla.spde.make.A(mesh.true,(kenya.coord@coords))

## The number of households sampled at each EA
house_prob <- c(31 ,47 ,51 ,63 ,65 ,62 ,29 ,16)
house_prob <- house_prob/sum(house_prob)
house_prob <- data.frame(values = c(4:11), prob = house_prob)
replicates <- round(house_prob$prob*400)

m <- c( rep(4, replicates[1]), rep(5, replicates[2]), 
        rep(6, replicates[3]), rep(7, replicates[4]), 
        rep(8, replicates[5]), rep(9, replicates[6]), 
        rep(10, replicates[7]), rep(11, replicates[8]) )

## The number of individuals tested in heach household
test_prob <- c(1376, 1223 , 418,  165 ,  76  , 22 ,
               7 ,   1 ,   1  ,  1,  2   , 1 )
test_prob <- test_prob/sum(test_prob)
test_prob <- data.frame(value=c(1:12), prob = test_prob)

##N: number tested in each household,
##y: number tested positive in each household,
##v: household random effect
N <- sample(x = test_prob$value, sum(m), replace=T, prob=test_prob$prob)
v <- rnorm(sum(m), sd = 0.1)
logit.prev <- rep( rep(logit(0.07), 
                       length(kenya.coord))+as.numeric(A%*%field.true), times = m) + v
y <- rbinom(sum(m), N, expit(logit.prev))


##Create dataframe for use in plotting and INLA.
obs <- as.data.frame(kenya.coord)[rep(seq_len(length(kenya.coord)),times=m),] %>% mutate(observed = y/N,
cluster1 = rep(c(1:length(kenya.coord)),times=m),N=N, present = as.factor(y))
colnames(obs)[1:2] <- c("LONGNUM", "LATNUM")
##Let's plot our simulated data on the true field.
last_plot() + geom_point(data=obs,aes(x=(LONGNUM),y=(LATNUM),colour=observed)) +
  scale_colour_gradient(limits=c(0,1), low = "white", high = "black")

# ggsave("plots/prevalence_field.pdf",device="pdf")

##**Changed plot.clip.true to true.prev.plot.final


## Fit the model in INLA to estimate the prevalence.
##NB: lots of things recomputed so that this section still works when
##  we don't simulate the data!

loc <- kenya.coord@coords[rep(seq_len(length(kenya.coord)),times=m),]
mesh <-  inla.mesh.2d(loc.domain =loc, offset=c(2,5), max.edge=c(0.5,10))
spde <- inla.spde2.matern(mesh)
A.est <- inla.spde.make.A(mesh,loc)

##Specify the formula in INLA
formula = y ~ -1 + intercept + f(field, model=spde) + f(eps,model="iid") 

##Get index for mesh. Will be used in creating
##inla.stack object for estimation.
spde.index <- inla.spde.make.index("field",spde$n.spde)

##Create the inla.stack object for model estimation.
stack.est <- inla.stack(data = list(y = as.numeric(obs$present)-1, N=N),
                        A = list(A.est,1), effects=list(c(spde.index, list(intercept=1)), 
                                                        list(eps = seq_len(sum(m)), cluster1=obs$cluster1 )), tag="est")



## We want the average prevelance per polygon,
## which we are going to compute by pure Monte Carlo.
## We're going to do this for every single region in adm1.
# The "true prevalence average" vs MC histogram is for each region is
# in the "plots/" directory

## Speed up the inference by having a warm start (with now prediction)
result.keep <- inla(formula, family="binomial", Ntrials=N, 
                    data=inla.stack.data(stack.est), 
                    control.predictor = list(compute=TRUE,
                                             A=inla.stack.A(stack.est)),
                    verbose=FALSE,control.compute = list(config=TRUE))

##Again reparameterize for interpretability/plotability
result.keep.field <- inla.spde2.result(result.keep,'field',
                                       spde,do.transform = TRUE)


##Compute posterior means of params defining field.
inla.emarginal(function(x) x,result.keep.field$marginals.kappa[[1]])
inla.emarginal(function(x) x,result.keep.field$marginals.variance.nominal[[1]])
inla.emarginal(function(x) x,result.keep.field$marginals.range.nominal[[1]])

## Plot the predicted median with related stuff...
proj = inla.mesh.projector(mesh,dims=c(300,300))
S = result.keep.field$summary.values$`0.5quant`
S25 = result.keep.field$summary.values$`0.025quant`
S75 = result.keep.field$summary.values$`0.975quant`

field.plot.pred = inla.mesh.project(proj,field=S,dims=c(200,200))
field.plot.pred25 = inla.mesh.project(proj,field=S25,dims=c(200,200))
field.plot.pred75 = inla.mesh.project(proj,field=S75,dims=c(200,200))

x = proj$x
y = proj$y
all = (expand.grid(x,y,KEEP.OUT.ATTRS = F))
names(all) = c("x","y")
field.pred.plotdata = all %>%mutate(field = as.numeric(field.plot.pred), 
                                    prev = exp(field)/(1+exp(field)) )
field.pred.plotdata25 = all %>%mutate(field25 = as.numeric(field.plot.pred25), 
                                      prev25 = exp(field25)/(1+exp(field25)) )
field.pred.plotdata75 = all %>%mutate(field75 = as.numeric(field.plot.pred75), 
                                      prev75 = exp(field75)/(1+exp(field75)) )
##Make these SpatialPointsDataFrame s
coordinates(field.pred.plotdata)=~x+y 
coordinates(field.pred.plotdata25)=~x+y
coordinates(field.pred.plotdata75)=~x+y

##adm0 is a SpatialPolygonsDataFrame
##We want to assign the same CRS() proj4string
##to our plot data.

field.pred.plotdata@proj4string = adm0@proj4string 

##We need to subset the field by the area we want to look at.

field.pred.plotdata.keep = field.pred.plotdata[adm0,] 
field.pred.plotdata25@proj4string = adm0@proj4string 
field.pred.plotdata25.keep = field.pred.plotdata25[adm0,]
field.pred.plotdata75@proj4string = adm0@proj4string 
field.pred.plotdata75.keep = field.pred.plotdata75[adm0,]

ggplot(as.data.frame(field.pred.plotdata.keep),aes(x,y)) + geom_tile(aes(fill=field)) #+scale_fill_gradient(limits=10*c(-0.1,0.1))
tmp = as.data.frame(kenya.coord@coords)
colnames(tmp) <- c("LONGNUM", "LATNUM")
last_plot() + geom_point(data=tmp,aes(x=LONGNUM,y=LATNUM))
# ggsave("plots/spatial_field_estimated.png",device="png")

ggplot(as.data.frame(field.pred.plotdata25.keep),aes(x,y)) + geom_tile(aes(fill=field25)) #+scale_fill_gradient(limits=10*c(-0.1,0.1))
tmp25 = as.data.frame(kenya.coord)
colnames(tmp25) <- c("LONGNUM", "LATNUM")
last_plot() + geom_point(data=tmp25,aes(x=LONGNUM,y=LATNUM))
#ggsave("plots/spatial_field_estimated25",device="png")

ggplot(as.data.frame(field.pred.plotdata75.keep),aes(x,y)) + geom_tile(aes(fill=field75)) #+scale_fill_gradient(limits=10*c(-0.1,0.1))
tmp75 = as.data.frame(kenya.coord)
colnames(tmp75) <- c("LONGNUM", "LATNUM")
last_plot() + geom_point(data=tmp75,aes(x=LONGNUM,y=LATNUM))
#ggsave("plots/spatial_field_estimated75",device="png")


### Compute adm1 estimates
medprev <- matrix(NA, nrow = length(adm1), ncol=2)
colnames(medprev) <- c("MC", "True.MC")
obs$weight <- obs$prob <- NA

#for(region.number in c(1:length(adm1))){
for(region.number in c(1:3)){
  ##Subset by region
  region.of.interest = adm1[region.number,]
  
  ##sample GPS points within region.of.interest
  n.pred =100
  points.mc = spsample(region.of.interest,type = "random",n = n.pred)
  
  ##Make our predict stack
  A.pred = inla.spde.make.A(mesh,points.mc@coords)
  stack.pred = inla.stack( data = list(y = NA,N=1), A = list(A.pred,1), 
                           effects=list(c(spde.index,list(intercept=1)),
                                        list(eps=rep(NA,n.pred),cluster1=rep(NA,n.pred)) ),tag="pred")
  
  ##Remember we already created an inla.stack object for estimation
  stack.join = inla.stack(stack.pred, stack.est)
  
  ##Fit model & do prediction with INLA call
  ##restart = FALSE means we will not redo estimation
  
  result = inla(formula, family="binomial", Ntrials=N, 
                data=inla.stack.data(stack.join), 
                control.predictor = list(compute=TRUE,
                                         A=inla.stack.A(stack.join)),
                verbose=FALSE,control.compute = list(config=TRUE),
                control.mode = list(theta=result.keep$mode$theta,
                                    restart=FALSE))
  
  ## extract the indices for the 1000 prediction points
  index=inla.stack.index(stack.join,"pred")
  
  mean.pred = result$summary.fitted.values[index$data,"mean"]
  
  ## Sample from the posterior 100 times and do a MC estimate of the 
  ## area average prevalence
  ## Note: control.compute = list(config=T) must be used
  ## in INLA call to sample from posteriors.
  samps = inla.posterior.sample(100,result)
  means = lapply(X=samps, FUN= function(x) { mean(expit(x$latent[index$data])) })
  
  ## For convenience, redefine kenya.coord@proj4string
  ## to make plots at adm1 level
  kenya.coord@proj4string = adm1@proj4string
  
  
#  fname = paste("plots/region",region.number,".png",sep="")
#  png(filename=fname)
  
  A.true.mc = inla.spde.make.A(mesh.true, points.mc@coords)
  tmp = length(kenya.coord[adm1[region.number,],])
  
  hist(unlist(means),breaks=12,prob=T,xlim=c(0,0.2),xlab="Prevalence",main = paste(tmp," of 400 observations in region ",region.number,sep=""))
  abline(v = mean(expit(as.numeric(mu+A.true.mc%*%field.true))),col="red",lwd=2,lty=2)
  abline(v=mean(unlist(means)),col="blue",lty=1)
  legend(x="topright",legend=c("True average","Posterior MC estimate"),
         col = c('red', 'blue'), lty = c(2,1))
#  dev.off()
  medprev[region.number, 'MC'] <- median(unlist(means))
  medprev[region.number, 'True.MC'] <- median(expit(mu+as.numeric(A.true.mc%*%field.true)))
  
}  

#pdf("plots/PosteriorMC.pdf")
mapvariable(medprev[,"MC"], adm1, main = "", nlevels=5,lower = 0.03, 
            upper = 0.1)
#dev.off()

#pdf("plots/PosteriorMCtrue.pdf")
mapvariable(medprev[,"True.MC"], adm1, main = "", nlevels=5,lower = 0.03,
            upper = 0.1)
#dev.off()

#############
##Design+ICAR models
##############

## Get weights for simulated data.

for(i in 1:400){
  p <- rep((400/46034)*m[i]/100, m[i])
  obs[obs$cluster1 == i, 'prob'] <- p
  obs[obs$cluster1 == i, 'weight']<- 1/p
}

print(paste("Total of Weights is ",sum(obs$weight), sep =""))

## Lay polygon over points
library(geosphere)
library(survey)
coordinates(obs) <- ~LONGNUM + LATNUM
adm1.poly <- SpatialPolygons(adm1@polygons,
                             proj4string=CRS("+proj=longlat +ellps=WGS84"))
proj4string(obs) <- proj4string(adm1.poly)
adm1.points2poly.2003 <- over(obs,adm1.poly)

obs.adm1 <- obs
obs.adm1$Polygon <- adm1.points2poly.2003


##Assign all points to polygon
miss.index <- is.na(obs.adm1@data$Polygon)
if(sum(miss.index != 0)){
miss.poly <- dist2Line(unique(coordinates(obs.adm1)[miss.index,]), adm1.poly)
obs.adm1@data[miss.index, 'Polygon'] <- miss.poly[ ,'ID']
}

##Get household ID.
for(i in 1:400){
  index <- which(obs.adm1@data$cluster1 == i)
  obs.adm1@data[index, 'HH_id'] <- 1:m[i]
}

##Create survey design

obs.hh <- matrix(NA, nrow= sum(obs.adm1$N), ncol = 6)
obs.hh <- as.data.frame(obs.hh)
colnames(obs.hh) <- colnames(obs.adm1@data)[-c(3,4)]

house.no <- 0
person.no <- 0
counter <- 0
for(i in 1:400){
  n.house <- m[i] 	## Number households in cluster i
  
  for(j in 1:n.house){
    house.no <- house.no + 1 ## unique household number
    household.data <- obs.adm1@data[house.no,]
    
    for(k in 1:household.data$N){
      counter <- counter + 1 ## Counter for prevalence
      
      
      person.no <- person.no + 1 ## Unique person number
      
      obs.hh$cluster1[person.no] <- household.data$cluster1
      obs.hh$Polygon[person.no] <- household.data$Polygon
      obs.hh$HH_id[person.no] <- household.data$HH_id
      obs.hh$prob[person.no] <- household.data$prob
      obs.hh$weight[person.no] <- household.data$weight
      
      if(counter <= (as.numeric(household.data$present)-1)){
        obs.hh$observed[person.no] <- 1
      }else{ obs.hh$observed[person.no] <- 0}
      
    }
    counter <- 0				
  }
}

## Get survery design object.
## Calculated design-based small area estimate.
options(survey.lonely.psu = "adjust")
adm1.des.2003 <- svydesign(ids= ~cluster1 + HH_id, weights = ~weight, data = obs.hh)

res.adm1.2003w <- svyby(~observed, ~Polygon, adm1.des.2003, svymean)
which(!(1:47 %in% res.adm1.2003w$Polygon))
res.adm1.2003w <- rbind(res.adm1.2003w[1:8,],
                        c(9, NA, NA), res.adm1.2003w[9:46,])

for(i in 1:47){
  if(!is.na(res.adm1.2003w$observed[i])){
    if(res.adm1.2003w$observed[i] %in% c(0) ){
      N <- length(obs.hh[obs.hh$Polygon== i,]$weight)
      Y.raw <- sum(obs.hh[obs.hh$Polygon==i, ]$observed)
      p.fix <- (Y.raw + 0.5)/(N + 1)
      res.adm1.2003w[i, 'observed'] <- p.fix
      res.adm1.2003w[i, 'se'] <- sqrt(p.fix*(1-p.fix)/N)
    }
  }
}

## Fixed variance, spatial smoothing
n.poly.adm1 <- 47
prev.adm1.poly.2003 <- matrix(NA, nrow=n.poly.adm1, ncol = 1)
prev.adm1.poly.2003 <- as.data.frame(prev.adm1.poly.2003)
colnames(prev.adm1.poly.2003) <- "p.i"

prev.adm1.poly.2003$p.i <- res.adm1.2003w$observed
prev.adm1.poly.2003$logit.pi <- logit(prev.adm1.poly.2003$p.i)
prev.adm1.poly.2003$w.i <- res.adm1.2003w$se^2
prev.adm1.poly.2003$v.i <- prev.adm1.poly.2003$w.i/(prev.adm1.poly.2003$p.i^2*(1-prev.adm1.poly.2003$p.i)^2)
prev.adm1.poly.2003$logit.prec <- 1/prev.adm1.poly.2003$v.i
prev.adm1.poly.2003$struct <- 1:n.poly.adm1
prev.adm1.poly.2003$unstruct <- 1:n.poly.adm1

##Make a .graph file so that INLA knows the neighborhood
##structure of your regions.

#neighb <- poly2nb(SpatialPolygons(adm1@polygons, proj4string = adm1@proj4string))
#nb2INLA(file = 'Kenyaadm1.graph', neighb)

formula = logit.pi ~ 1 +f(unstruct,model='iid',param=c(0.5,0.008)) + 
  f(struct,model='besag',adjust.for.con.comp=TRUE,constr=TRUE,
    graph='Kenyaadm1.graph')
mod1.adm1.2003 <- inla(formula,family="gaussian",
                       data=prev.adm1.poly.2003, control.family=
                         list(hyper=list(prec=list(initial=log(1),fixed=TRUE))), 
                       scale=logit.prec,control.predictor = list(compute=TRUE))

fixed.med.2003 <- rep(mod1.adm1.2003$summary.fixed[1,4],n.poly.adm1)
random.iid.2003 <- mod1.adm1.2003$summary.random$unstruct[,5]
random.besag.2003 <- mod1.adm1.2003$summary.random$struct[,5]

odds.2003 <- pred.2003 <- NA
for(i in 1:47){
  tmp <- inla.rmarginal(1000, mod1.adm1.2003$marginals.linear.predictor[[i]])
  pred.2003[i] <- median(expit(tmp))
  odds.2003[i] <- median(exp(tmp))
}


prev.adm1.poly.res.2003 <- cbind(prev.adm1.poly.2003,fixed.med.2003,
                                 random.iid.2003,random.besag.2003, pred.2003, odds.2003)

#pdf("plots/PosteriorICAR.pdf")
mapvariable(pred.2003, adm1.poly, ncut = 1000, nlevels = 10, lower = 0, upper = 0.15,
            main = "")
#dev.off()
#pdf("plots/comparison.pdf")
par(mfrow=c(2,2))
plot(medprev[,"MC"]~medprev[,"True.MC"],xlim=c(0.03,0.1),
     ylim=c(0.03,0.1),xlab="Truth",ylab="SPDE")
abline(0,1)
plot(pred.2003~medprev[,"True.MC"],xlim=c(0.03,0.1),
     ylim=c(0.03,0.1),xlab="Truth",ylab="Smoothed Design")
abline(0,1)
plot(unlist(res.adm1.2003w["observed"])~medprev[,"True.MC"],
     xlim=c(0.03,0.1),
     ylim=c(0.03,0.1),xlab="Truth",ylab="Raw Design")
abline(0,1)
plot(pred.2003~medprev[,"MC"],xlim=c(0.03,0.1),
     ylim=c(0.03,0.1),ylab="Smoothed Design",xlab="SPDE")
abline(0,1)
#dev.off()
