#####################################################################
# Script for TimeStamp algorithm
# (c) Rosemary Braun <rbraun@northwestern.edu>
######################################################################


#=====================================================================
# time conversion functions
#=====================================================================

time2dectime <- function(time){
	if(any(c("POSIXct","POSIXt")%in%class(time))){
		# if we have POSIX input, make it a string to meet our (bad)
		# assumptions later...
		time <- strftime(time,"%H:%M")
	}
	if (is.numeric(time)) {
		# horrible unchecked assumption that if input is numeric 
		# it's already decimal 24h time and can be returned
		dectime <- time
	} else {
		# horrible unchecked assumption that otherwise we have
		# times in %H:%M string format -- note that 2:30PM will fail!
		hr_min <- strsplit(as.character(time),":")
		hr <- as.numeric(sapply(hr_min,`[`,1))
		min <- as.numeric(sapply(hr_min,`[`,2))
		min[is.na(min)] <- 0
		dectime <- hr+min/60
	}
	return(dectime)
}

time2angle <- function(time){
	dectime <- time2dectime(time)
	timeAngle <- (dectime%%24)*2*pi/24
	return(timeAngle)
}

time2XY <- function(time){
	a <- time2angle(time)
	XYtime <- zapsmall(cbind(
		timeX=sin(a),
		timeY=cos(a)
	))
	return(XYtime)
}

XY2dectime <- function(XYtime){
	(atan2(XYtime[,1],XYtime[,2])%%(2*pi))*(24/(2*pi))
}

#=====================================================================
# Splitting & centering the expression data for each subject
#=====================================================================

recalibrateExprs <- function(exprMat,subjectIDs){
	# exprMat is a genes * subjects matrix
	# subjectIDs a vector of subject IDs
	exprs <- as.data.frame(t(exprMat))
	exprs <- split(exprs,subjectIDs)
	FCs <- lapply(exprs,scale,center=T,scale=F)
	FCs <- t(do.call(rbind,FCs))[,colnames(exprMat)]
	return(FCs)
}	


#=====================================================================
# FNS FOR MODELING X&Y CLOCK COORDS
#=====================================================================

trainTimeStamp <- function(expr,subjIDs,times,trainFrac=0.5,a=0.5,s=NULL,plot=FALSE,recalib=FALSE,...){
	stopifnot(require(glmnet))
	out <- list()
	# select trainFrac subjects to be the training set
	subjects <- unique(subjIDs)
	trainSubj <- sample(subjects,size=round(trainFrac*length(subjects)),replace=FALSE)
	train <- subjIDs%in%trainSubj
	out$train <- train
	if (recalib) {
		# normalize exprs to FC from mean on a per-subject basis
		x <- t(recalibrateExprs(expr, subjIDs))
	} else {
		x <- t(expr)
	}
	y <- time2XY(times)
	# do the fit
	out$cv.fit <- cv.glmnet(x[train,],y[train,],keep=T,alpha=a,family="mgaussian",...)
	if(is.null(s)){
		s <- out$cv.fit$lambda.min
	}
	out$coef <- as.matrix(do.call(cbind,coef(out$cv.fit, s=s)))
	out$coef <- out$coef[rowSums(out$coef!=0)>0,][-1,]
	out$pred <- XY2dectime(predict(out$cv.fit,x,s=s)[,,1])
	if(plot){
		dectime <- time2dectime(times)
		errplot(dectime,out$pred,col=2-train)
	}
	return(out)
}

predTimeStamp <- function(timestamp,newx=NULL,s=timestamp$cv.fit$lambda.min){
	stopifnot(require(glmnet))
	if(is.null(newx)){
		return(timestamp$pred)
	}else{
		newx <- t(newx)
		XY2dectime(predict(timestamp$cv.fit,newx,s=s)[,,1])
	}
}

#=====================================================================
# FNS FOR CALCULATING PREDICTION ERROR MODULO 24
#=====================================================================

timeErr <- function(trueTimes,predTimes,...){
	predTimes <- predTimes%%24
	trueTimes <- trueTimes%%24
	timeDiffs <- abs(predTimes-trueTimes)
	timeDiffs <- pmin(timeDiffs,24-timeDiffs)
	return(timeDiffs)
}

timeErrSigned <- function(trueTimes,predTimes,...){
	predTimes <- predTimes%%24
	trueTimes <- trueTimes%%24
	timeDiffs <- predTimes-trueTimes
	timeDiffs[which(timeDiffs>12)]<-timeDiffs[which(timeDiffs>12)]-24
	timeDiffs[which(timeDiffs<(-12))]<-timeDiffs[which(timeDiffs<(-12))]+24	
	return(timeDiffs)
}

#=====================================================================
# FNS FOR PLOTTING PREDICTED CLOCK COORDS
#=====================================================================

errplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)"){
	# plot predicted vs. true time of day , with grey bands at 2 & 4h MOE
	plot(trueHr,predHr,xlab="",ylab="",main=main,yaxs='i', xaxp=c(0,24,6),yaxp=c(0,24,6),col=col,pch=pch,xlim=c(0,24),ylim=c(0,24))
	trueHr <- (trueHr+240)%%24
	predHr <- (predHr+240)%%24
	mtext("True",side=1,line=2.2,cex=0.8)
	mtext("Predicted",side=2,line=2.4,cex=0.8)
	abline(a=0,b=1,col="grey");
	abline(a=24,b=1,col="grey");
	abline(a=-24,b=1,col="grey");
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4),border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2),border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4)+24,border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2)+24,border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4)-24,border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2)-24,border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	points(trueHr,predHr,col=col,pch=pch)
}


tolplot <- function(trueHr,predHr,add=FALSE,col=1){
	# plot % correct by tolerance, dropping lines at median, 90quantile, & 1sd 
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
		plot(hrsoff,100-fracacc,xlim=c(0,12),
			main="Absolute error CDF", xlab="",ylab="")
		mtext("correct to within (hrs)",side=1,line=2.2,cex=0.8)
		mtext(paste("% correct (N = ",length(hrerr),")",sep=""),side=2,line=2.4,cex=0.8)
		abline(h=c(50,80,100),v=c(median(abs(hrerr)),quantile(abs(hrerr),0.8)),col="grey")
		asTime <- function(x){
			hrs <- (x%%24)%/%1; 
			min <- round(60*(x%%1)); 
			if(min==60){hrs <- (hrs+1)%%24; min<-0}
			sprintf("%2i:%02i",hrs,min)
		}
		text(quantile(abs(hrerr),0.8),30,asTime(quantile(abs(hrerr),0.8)),cex=0.9)
		text(median(abs(hrerr)),10, asTime(median(abs(hrerr))),cex=0.9)
	}
	points(hrsoff, 100-fracacc, col=col)
	lines(hrsoff, 100-fracacc, col=col)
	norm.fracacc <- (100-fracacc)/100
	norm.hrsoff <- hrsoff/12
	auc <- sum(norm.fracacc[-1]*diff(norm.hrsoff))
	if(!add){
		text(10,20, sprintf("nAUC=%.2f",auc),cex=1,font=2)
	}
	invisible(auc)
}


predplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)",...){
	# do both of the above -- set par(mfrow=c(2,1)) or (1,2) first!
	opar <- par(xpd=F,mar=c(4,4,3,1))	
	on.exit(par(opar))
	out <- timeErr(trueHr,predHr,plot=T,col,pch,main) 
	tolplot(trueHr,predHr)
	invisible(out)
}

