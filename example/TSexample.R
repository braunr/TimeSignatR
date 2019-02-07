#####################################################################
# TimeSignature example / Fig 1 reproduction
# (c) Rosemary Braun <rbraun@northwestern.edu>
######################################################################


#=====================================================================
# Loading data and packages -- 
#   assumes TimeStampR/example is the current working directory
#=====================================================================

set.seed(194) 

library(limma)
library(glmnet)

source("../R/TimeStampFns.R")

# This file contains expression data for all three datasets,
# time, subject, and condition metadata, as well as previously
# obtained ZeitZeiger and PLSR predictions for comparison
load("DATA/TSexampleData.Rdata")

#=====================================================================
# convenience functions 
#=====================================================================

Rint <- function(...){
	Reduce(intersect,list(...))
}

asNum <- function(z){
	as.numeric(as.character(z))
}

TScoef <- function(timestamp,s=timestamp$cv.fit$lambda.min){
	out <- as.matrix(do.call(cbind,coef(timestamp$cv.fit, s=s)))
	out <- out[rowSums(out!=0)>0,][-1,]
	colnames(out) <- c("sunX","sunY")
	return(out)
}

nwhich <- function(x){
	if(is.null(names(x))){
		which(x)
	}else{
		names(which(x))
	}
}

asTime <- function(x){
	hrs <- (x%%24)%/%1; 
	min <- round(60*(x%%1)); 
	hrs[min==60] <- (hrs+1)%%24
	min[min==60] <- 0
	sprintf("%2i:%02i",hrs,min)
}

asTime0 <- function(x){
	hrs <- (x%%24)%/%1; 
	min <- round(60*(x%%1)); 
	hrs[min==60] <- (hrs+1)%%24
	min[min==60] <- 0
	sprintf("%02i:%02i",hrs,min)
}

#=====================================================================
# new plotting routines
#=====================================================================

timeErrPlot <- function(trueTimes,predTimes,...){
	errplot(trueTimes,predTimes,...)
	timeDiffs <- timeErr(trueTimes,predTimes)
	return(timeDiffs)
}

timeErrSignedPlot <- function(trueTimes,predTimes,...){
	errplot(trueTimes,predTimes,...)
	timeDiffs <- timeErrSigned(predTimes,trueTimes)
	return(timeDiffs)
}

tolplot <- function(trueHr,predHr,add=FALSE,col=1,...){
	# plot % correct by tolerance, without lines 
	predTimes <- predHr%%24
	trueTimes <- trueHr%%24
	hrerr <- abs(predTimes-trueTimes)
	hrerr <- hrerr[!is.na(hrerr)]
	hrerr <- pmin(hrerr,24-hrerr)
	hrsoff <- seq(0,12,length=49)
	fracacc <- sapply(hrsoff,function(hrtol){
		100*sum(abs(hrerr)>hrtol)/length(hrerr)
	})
	if(!add){
		col=1
		plot(hrsoff,100-fracacc,xlim=c(0,12), type="n",
			main="Absolute error CDF", xlab="",ylab="")
		mtext("correct to within (hrs)",side=1,line=2.2,cex=0.8)
		mtext(paste("% correct (N = ",length(hrerr),")",sep=""),side=2,line=2.4,cex=0.8)
		abline(a=0,b=100/12,col="grey")
		asTime <- function(x){
			hrs <- (x%%24)%/%1; 
			min <- round(60*(x%%1)); 
			if(min==60){hrs <- (hrs+1)%%24; min<-0}
			sprintf("%2i:%02i",hrs,min)
		}
	}
	lines(hrsoff, 100-fracacc, col=col,lwd=1.5,...)
	norm.fracacc <- (100-fracacc)/100
	norm.hrsoff <- hrsoff/12
	auc <- sum(norm.fracacc[-1]*diff(norm.hrsoff))
	return(list(auc=auc,mederr=median(abs(hrerr))) )
}

predplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)",...){
	# do both of the above -- set par(mfrow=c(2,1)) or (1,2) first!
	opar <- par(xpd=F,mar=c(4,4,3,1))	
	on.exit(par(opar))
	out <- timeErrPlot(trueHr,predHr,col,pch,main) 
	out <- tolplot(trueHr,predHr)
	invisible(out)
}

#=====================================================================
# Training with Moller subset data
#=====================================================================

# Note - the subset chosen for training is already indicated in
# the data, but we will refresh it here for completeness:
	set.seed(194)
	trainFrac <- 0.5
	moller.char <- all.meta[all.meta$study=="TrTe",]
	trainSet <- unique(moller.char$ID)
	trainSet <- sample(trainSet,size=round(trainFrac*length(trainSet)),replace=FALSE)
	moller.char$train <- as.numeric(moller.char$ID%in%trainSet)
	all.meta[all.meta$study=="TrTe","train"] <- moller.char$train
# When fixing the foldid's for opimizing alpha, the following 
# randomly-generate foldids were used; we repeat this here too
# for illustration and reproducibility:
	train.foldid <- sample(trainSet)
	train.foldid <- sample(rep(seq(10),length=sum(moller.char$train)))

trainDat <- all.expr[,all.meta$train==1]
trainSubjs <- all.meta[all.meta$train==1,"ID"]
trainTimes <- all.meta[all.meta$train==1,"LocalTime"]
# within-subject normalization using all timepoints
trainWSN <- recalibrateExprs(trainDat, trainSubjs)

TSorig <- trainTimeStamp(
	expr=trainWSN, # use within-subject normalized data
	subjIDs=trainSubjs,
	times=trainTimes,
	trainFrac=1, # no need to subset training samples; already done!
	recalib=FALSE, # no need to within-subj normalize; already done!
	a=0.5, s=exp(-1.42), # penalty params as used in the paper
	foldid=train.foldid, # foldIDs as in the paper
	plot=FALSE 
)


#=====================================================================
# Testing with 2-Point Calibration
#     Rather than using the recalibrateExprs function, which uses
#     all data for a given subject, we mimic the case where each 
#     subject has only two timepoints by attempting to find the 
#     "antipodal" point within each subject for each sample.
#=====================================================================


# First, get all the times at which each subject has data
all.times.by.subj <- split(all.meta[,c("samp","LocalTime")],all.meta$sID)

# Function to identify the antipodal point in that list
findAntipode <- function(time.by.subj){
	antipodes <- lapply(rownames(time.by.subj),function(thisSamp){
		thisTime <- time.by.subj[thisSamp,"LocalTime"]
		availTimes <- time.by.subj$LocalTime
		names(availTimes) <- rownames(time.by.subj)
		availTimes[thisSamp] <- NA # don't choose self
		antipTime <- thisTime+12
		antipInd <- which.min(timeErr(availTimes,antipTime))
		antipSamp <- rownames(time.by.subj)[antipInd]
		out <- list(
			samp=thisSamp,
			asamp= antipSamp,
			LocalTime=thisTime,
			aLocalTime= availTimes[antipSamp]
		)
		out$LocalTimediff <- timeErr(out$LocalTime,out$aLocalTime)
		return(data.frame(out,stringsAsFactors=F))
	})
	return(do.call(rbind,antipodes))
}

# Put the antipodes together and name them appropriately
allAntipodes <- do.call(rbind,lapply(all.times.by.subj,findAntipode))
rownames(allAntipodes) <- allAntipodes$samp
# Put it back in the original order
allAntipodes <- allAntipodes[rownames(all.meta),]
# Get a vector with the name of the sample antipode
sampAntipodes <- allAntipodes$asamp
names(sampAntipodes) <- allAntipodes$samp 

# And now calibrate!
all.calib <- sapply(colnames(all.expr),function(samp){
	centerExpr <- rowMeans(all.expr[,c(samp,sampAntipodes[samp])],na.rm=T)
	return(all.expr[,c(samp)]-centerExpr)
})
# When NA, set to 0 -- ie, the gene is assumed to be flat, ie same as 
# its circadian mean at that point (NB: does not imply unexpressed)
all.calib.NA0 <- all.calib
all.calib.NA0[is.na(all.calib.NA0)] <- 0

#=====================================================================
# Making TimeSignature predictions:
#	In paper: a=0.5, s=exp(-1.42), 2-sample antipodal calibration
#=====================================================================

all.meta$TSpred <- predTimeStamp(TSorig,newx=all.calib.NA0,s=exp(-1.42))

#=====================================================================
# Making the Fig 1 plot
#=====================================================================

statLegend <- function(){
	legend("bottomright", bty="n", text.font=c(3,2,2,2), lwd=1.5, cex=0.9, lty=c(0,1,2,3), col=c("white","black","purple","cyan4"), text.col=c("black","black","purple","cyan4"),
		legend=c("      AUC,  med",sprintf(
			"%s: %.2f, %s",
			c("TS","PL","ZZ"),
			c(TSauc$auc,PLSauc$auc,ZZauc$auc),
			asTime(c(TSauc$mederr,PLSauc$mederr,ZZauc$mederr))
		)))
}

par(mfcol=c(2,4))

# split up the results by study
all.meta.split <- split(all.meta, all.meta$study)

# Train/Test data (Moller) - TEST SAMPLES ONLY 
TSauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
	predplot(LocalTime,TSpred,plot=T,main="Test Set (Moller &al)",col=cond))
PLSauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
	tolplot(LocalTime, PLSpredLT,add=T,col="purple",lty=2))
ZZauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
	tolplot(LocalTime, ZZpredLT,add=T,col="cyan4",lty=3))
statLegend()

# V1 validation data (Archer) 
TSauc <- with(all.meta.split$V1,
	predplot(LocalTime,TSpred,plot=T,main="Validation V1 (Archer &al)",col=cond))
PLSauc <- with(all.meta.split$V1,
	tolplot(LocalTime, PLSpredLT,add=T,col="purple",lty=2))
ZZauc <- with(all.meta.split$V1,
	tolplot(LocalTime, ZZpredLT,add=T,col="cyan4",lty=3))
statLegend()

# V2 validation data (Arnardottir &al) 
TSauc <- with(all.meta.split$V2,
	predplot(LocalTime,TSpred,plot=T,main="Validation V2 (Arnardottir &al)",col=cond))
PLSauc <- with(all.meta.split$V2,
	tolplot(LocalTime, PLSpredLT,add=T,col="purple",lty=2))
ZZauc <- with(all.meta.split$V2,
	tolplot(LocalTime, ZZpredLT,add=T,col="cyan4",lty=3))
statLegend()

# V3 validation data (new RNA-seq) 
TSauc <- with(all.meta.split$V3,
	predplot(LocalTime,TSpred,plot=T,main="Validation V3 (new RNA-seq)",col=1))
PLSauc <- with(all.meta.split$V3,
	tolplot(LocalTime, PLSpredLT,add=T,col="purple",lty=2))
ZZauc <- with(all.meta.split$V3,
	tolplot(LocalTime, ZZpredLT,add=T,col="cyan4",lty=3))
statLegend()

