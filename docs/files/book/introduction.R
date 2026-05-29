#
# R CODE FOR REPRODUCING THE FIGURES AND ANALYSES IN JW'S 
#       "BAYESIAN AND FREQUENTIST REGRESSION ANALYSIS" CHAPTER 1
# CODE WRITTEN BY JON WAKEFIELD, UNLESS OTHERWISE STATED
#
#
# Prostate cancer example
#
library(lasso2)
data(Prostate)
attach(Prostate)
y <- Prostate$lpsa
lcavol <- Prostate$lcavol
lweight <- Prostate$lweight
age <- Prostate$age
lbph <- Prostate$lbph
svi <- Prostate$svi
lcp <- Prostate$lcp
gleason <- Prostate$gleason
pgg45 <- Prostate$pgg45
#
# Fig 1.1
#
pdf("prostfig1.pdf",height=7,width=4.5)
par(mar=c(4,4,2,1)+.1)  # bottom/left/top/right
par(mfrow=c(4,2))
plot(lcavol,lpsa,xlab="log(can vol)",ylab="y")
lines(lowess(lcavol,lpsa))
plot(lweight,lpsa,xlab="log(weight)",ylab="y")
lines(lowess(lweight,lpsa))
plot(age,lpsa,xlab="age",ylab="y")
lines(lowess(age,lpsa))
plot(lbph,lpsa,xlab="log(BPH)",ylab="y")
lines(lowess(lbph,lpsa))
boxplot(lpsa~svi,xlab="SVI",ylab="y")
plot(lcp,lpsa,xlab="log(cap pen)",ylab="y")
lines(lowess(lcp,lpsa))
plot(gleason,lpsa,xlab="gleason",ylab="y")
lines(lowess(gleason,lpsa))
plot(pgg45,lpsa,xlab="PGS45",ylab="y")
lines(lowess(pgg45,lpsa))
dev.off()
#
# Fig 1.2
#
pdf("prostfig2.pdf",height=7,width=4.5)
par(mfrow=c(4,3))
par(mar=c(4,4,2,1)+.1)  # bottom/left/top/right
boxplot(lcavol~svi,ylab="log(can vol)",xlab="SVI")
plot(lcp,lcavol,ylab="log(can vol)",xlab="log(cap pen)")
lines(lowess(lcp,lcavol))
plot(gleason,lcavol,ylab="log(can vol)",xlab="gleason")
lines(lowess(gleason,lcavol))
plot(pgg45,lcavol,ylab="log(can vol)",xlab="PGS45")
lines(lowess(pgg45,lcavol))
plot(lweight,age,ylab="age",xlab="log(weight)")
lines(lowess(lweight,age))
plot(lbph,lweight,xlab="log(BPH)",ylab="log(weight)")
lines(lowess(lbph,lweight))
plot(lbph,age,ylab="age",xlab="log(BPH)")
lines(lowess(lbph,age))
boxplot(lcp~svi,ylab="lcp",ylab="log(cap pen)",
xlab="SVI")
boxplot(pgg45~svi,ylab="PGS45",xlab="SVI")
plot(gleason,lcp,ylab="log(cap pen)",xlab="gleason")
lines(lowess(gleason,lcp))
plot(lcp,pgg45,ylab="PGS45",xlab="log(cap pen)")
lines(lowess(lcp,pgg45))
plot(gleason,pgg45,ylab="PGS45",xlab="gleason")
lines(lowess(gleason,pgg45))
dev.off()
#
# Minnesota lung cancer and radon example
#
library(maps)
lung<-read.table("MNlung.txt", header=TRUE, sep="\t")
radon<-read.table("MNradon.txt", header=TRUE)
#
# This function allows the mapping of data at the county level in Minnesota
#
MNmap <- function(data, ncol=5, figmain="", digits=5, type="e",
                   lower=NULL, upper=NULL,udig=2) {
  if (is.null(lower)) lower <- min(data)
  if (is.null(upper)) upper <- max(data) 
  if (type=="q"){p <- seq(0,1,length=ncol+1)
                 br <- round(quantile(data,probs=p),2)}
  if (type=="e"){br <- round(seq(lower,upper,length=ncol+1),udig)}
  shading <- gray((ncol-1):0/(ncol-1))
  data.grp <- findInterval(data,vec=br,rightmost.closed=T,all.inside=T)
  data.shad <- shading[data.grp]
  map("county", "minnesota", fill=TRUE, col=data.shad)
  leg.txt<-paste("[",br[ncol],",",br[ncol+1],"]",sep="")
  for(i in (ncol-1):1){
    leg.txt<-append(leg.txt,paste("[",format(br[i],digits=2),",",
             format(br[i+1],digits=2),")",sep=""),)
  }
  leg.txt<-rev(leg.txt)
  legend(-91.9,46.7, legend=leg.txt,fill=shading,bty="n",ncol=1,text.width=1)
  title(main=figmain,cex=1.5)
  invisible()
}
#
# Form observed counts and expected numbers
#
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
SMR <- Obs/Exp
#
# Fig 1.3 
#
pdf("SMRlung_radon.pdf", height=5.5, width=4.5)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1)+.1)  # bottom/left/top/right
MNmap(SMR, ncol=8, type="e", figmain="", lower=min(SMR), upper=max(SMR))
dev.off() 
#
# Map the radon
#
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i] <- mean(radon[radon$county==i,2])
}
rad.avg[26] <- 0
rad.avg[63] <- 0
#
# Fig 1.4 
#
pdf("radonMN.pdf", height=5.5, width=4.5)
par(mar=c(1,1,1,1)+.1)  # bottom/left/top/right
MNmap(rad.avg, ncol=8, type="e", figmain="", lower=min(rad.avg), upper=max(rad.avg),udig=1)
dev.off()
#
# Fig 1.5 
#
pdf("MNassoc.pdf",width=4.5,height=4)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(Obs/Exp~rad.avg,xlab="Average radon (pCi/liter)",ylab="SMR")
lines(lowess(Obs/Exp~rad.avg))
dev.off()
x <- rad.avg
x[26] <- NA
x[63] <- NA
poismod <- glm(Obs~offset(log(Exp))+rad.avg,family="poisson")
quasimod <- glm(Obs~offset(log(Exp))+rad.avg,family=quasipoisson(link="log"))
#
# Pharmacokinetic Theophylline example 
#
library(nlme)
data(Theoph)
conc1 <- Theoph$conc[24:33]
time1 <- Theoph$Time[24:33]
#
# Fig 1.6
#
pdf("theoph-fig1.pdf",width=4.5,height=4)
plot(conc1~time1,xlab="Time (hours)",ylab="Concentration (mg/liter)",
     ylim=c(0,max(conc1)))
dev.off()
#
# Fig 1.7
#
pdf("pkcomp1.pdf",width=4.5,height=4.5)
library(diagram)
par(mar=c(.8,.8,.8,.8))
openplotmat()
elpos  <- coordinates(c(2,2))
textellipse(mid=elpos[1,],0.16,0.08,lab=c("Compartment 0"))
textempty(mid=elpos[3,],lab=c("DOSE"),box.col="white")
textellipse(mid=elpos[2,],0.16,0.08,lab=c("Compartment 1"))
# textempty(mid=elpos[4,],lab=c("CONCENTRATION"))
straightarrow(from=elpos[1,]+c(0.16,0),to=elpos[2,]+c(-0.185,0),arr.pos=1)
straightarrow(from=elpos[3,]+c(0,0.03),to=elpos[1,]+c(0,-0.11),arr.pos=1)
straightarrow(from=elpos[2,]+c(0,-0.085),to=elpos[4,]+c(0,0.055),arr.pos=1)
textplain(mid=c(0.48,.78),lab=expression(k[a]))
textplain(mid=c(0.75,.25),lab=expression(k[e])) 
dev.off()
#
# Dental example
#
library(nlme)
data(Orthodont)
attach(Orthodont)
age <- Orthodont$age
distace <- Orthodont$distance
#
# Fig 1.8
#
pdf("dental.pdf",height=5,width=4.5)
plot(age,distance,type="n",xlab="Age (years)",ylab="Distance (mm)")
for (i in 1:16){
    index <- (i-1)*4+seq(1:4)
    lines(age[index],distance[index])
}
for (i in 1:11){
    index <- 64+(i-1)*4+seq(1:4)
    lines(age[index],distance[index],lty=2)
}
legend("topleft",bty="n",legend=c("Boys","Girls"),lty=c(1,2))
dev.off()
#
# Spinal bone density example
#
bone <- read.table("spinalbonedata.txt", header=T)
bone <- subset(bone, sex=="fem")
bone$ethnic <- as.factor(bone$ethnic)
bone$groupVec <- 1
bone$seqno <- 1:nrow(bone)
library(lattice)
pdf(file="spinalbone.pdf",height=5,width=4.5)
lattice.options(default.theme = "col.whitebg")
xyplot(spnbmd~age|ethnic, data = bone, groups=idnum, layout=c(2,2), aspect=1,
       xlab="Age (years)",ylab="Spinal Bone Marrow Density",
       panel=function(x,y,subscripts,groups) {
         lpoints (x,y, pch=20, col=1, cex=.3)
         ids=unique(groups[subscripts])
         for (id in ids) {
            subset=groups[subscripts]==id
            llines(x[subset], y[subset], col=1, lty=1)
         }
       }
)
dev.off()



