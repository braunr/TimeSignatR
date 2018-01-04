#####################################################################
# Initial attempt at TimeStamp algorithm
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
# NB: erroneous error model!
#=====================================================================

trainTimeStamp <- function(expr,subjIDs,times,trainFrac=0.5,plot=FALSE,recalib=FALSE,...){
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
	out$cv.fit <- cv.glmnet(x[train,],y[train,],family="mgaussian")
	out$coef <- as.matrix(do.call(cbind,coef(out$cv.fit, s=out$cv.fit$lambda.min)))
	out$coef <- out$coef[rowSums(out$coef!=0)>0,][-1,]
	out$pred <- XY2dectime(predict(out$cv.fit,x,s=out$cv.fit$lambda.min)[,,1])
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
# FNS FOR PLOTTING PREDICTED CLOCK COORDS
#=====================================================================
#--- "CLOCK" PLOT ---
clockplot <- function(trueHr,predHr,IDs=1,col=rainbow(25)[round(trueHr%%24)+1],noweb=F){
	#opar <- par(mfrow=c(1,2),xpd=T)	
	opar <- par(xpd=T,mar=c(2,1.5,3,1.5))	
	on.exit(par(opar))
	web <- function(rads){
		a <- seq(0,2*pi,length=180)
		sapply(rads,function(r){lines(r*sin(a),r*cos(a),col="grey80",lwd=0.5)})
		hrs <- (0:23)*2*pi/24
		r <- 1.05*max(rads)
		sapply(hrs,function(a){segments(0,0,r*sin(a),r*cos(a),col="grey80",lwd=0.5)})
	}			
	rads <- as.numeric(factor(IDs))
	radmax <- max(rads)
	tx <- rads*sin(trueHr*2*pi/24)
	ty <- rads*cos(trueHr*2*pi/24)
	px <- rads*sin(predHr*2*pi/24)
	py <- rads*cos(predHr*2*pi/24)
	plot(tx,ty,type="n",axes=F,bg=col,xlab="",ylab="",main="True time (24h)",pch=21,xlim=c(-radmax,radmax),ylim=c(-radmax,radmax))
	text(
		sin((0:23)*2*pi/24)*(radmax*1.1),
		cos((0:23)*2*pi/24)*(radmax*1.1),
		labels=as.character(0:23)
	)
	if(!noweb){web(rads)}
	points(tx,ty,bg=col,pch=21)
	plot(px,py,type="n",axes=F,bg=col,xlab="",ylab="",main="Predicted time (24h)",pch=21,xlim=c(-radmax,radmax),ylim=c(-radmax,radmax))
	text(
		sin((0:23)*2*pi/24)*(radmax*1.1),
		cos((0:23)*2*pi/24)*(radmax*1.1),
		labels=as.character(0:23)
	)
	if(!noweb){web(rads)}
	points(px,py,bg=col,pch=21)
}

#--- PRED PLOT ---

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

timeErr <- function(trueTimes,predTimes,plot=F,...){
	if(plot){
		errplot(trueTimes,predTimes,...)
	}
	predTimes <- predTimes%%24
	trueTimes <- trueTimes%%24
	timeDiffs <- abs(predTimes-trueTimes)
	timeDiffs <- pmin(timeDiffs,24-timeDiffs)
	return(timeDiffs)
}

timeErrSigned <- function(trueTimes,predTimes,plot=F,...){
	if(plot){
		errplot(trueTimes,predTimes,...)
	}
	predTimes <- predTimes%%24
	trueTimes <- trueTimes%%24
	timeDiffs <- predTimes-trueTimes
	#timeDiffs <- pmin(timeDiffs,24-timeDiffs)
	timeDiffs[which(timeDiffs>12)]<-timeDiffs[which(timeDiffs>12)]-24
	timeDiffs[which(timeDiffs<(-12))]<-timeDiffs[which(timeDiffs<(-12))]+24	
	return(timeDiffs)
}



tolplot <- function(trueHr,predHr,add=FALSE,col=1){
	# plot % correct by tolerance, dropping lines at median, 90quantile, & 1sd 
	predTimes <- predHr%%24
	trueTimes <- trueHr%%24
	hrerr <- abs(predTimes-trueTimes)
	hrerr <- hrerr[!is.na(hrerr)]
	hrerr <- pmin(hrerr,24-hrerr)
	#hrerr <- cbind(predHr-trueHr, predHr-trueHr+24, predHr-trueHr-24)
	#hrerr <- apply(hrerr,1,function(z){z[which.min(abs(z))]})
	hrsoff <- seq(0,12,length=49)
	fracacc <- sapply(hrsoff,function(hrtol){
		100*sum(abs(hrerr)>hrtol)/length(hrerr)
	})
	if(!add){
		col=1
		plot(hrsoff,100-fracacc,xlim=c(0,12),
			#main="% correct by tolerance", xlab="",ylab="")
			main="Absolute error CDF", xlab="",ylab="")
		mtext("correct to within (hrs)",side=1,line=2.2,cex=0.8)
		mtext(paste("% correct (N = ",length(hrerr),")",sep=""),side=2,line=2.4,cex=0.8)
		#abline(h=c(50,90,100),v=c(median(abs(hrerr)),sd(hrerr),quantile(abs(hrerr),0.9)),col="grey")
		abline(h=c(50,80,100),v=c(median(abs(hrerr)),quantile(abs(hrerr),0.8)),col="grey")
		#text(sd(hrerr),0,round(sd(hrerr),2),cex=0.7)
		#text(quantile(abs(hrerr),0.9),0,round(quantile(abs(hrerr),0.9),2),cex=0.7)
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
	#invisible(hrerr)
	invisible(auc)
}

predplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)",...){
	# do both of the above -- set par(mfrow=c(2,1)) or (1,2) first!
	opar <- par(xpd=F,mar=c(4,4,3,1))	
	on.exit(par(opar))
	out <- timeErr(trueHr,predHr,plot=T,col,pch,main) #errplot(trueHr,predHr,col,pch,main)
	tolplot(trueHr,predHr)
	invisible(out)
}

#--- STARBURST "GLOW" PLOT ---
glowPlot <- function(predX,predY,trueHr,rad=2,main=NULL){
	xlim=rad*c(-1,1);ylim=rad*c(-1,1)
	plot(predX,predY,xlim=xlim,ylim=ylim,main=main,
		type="n",xlab="",ylab="",axes=F)
	invisible(sapply(0:23,function(hr){
		wedgecol=adjustcolor(rainbow(24)[hr+1],alpha=0.4)
		polygon(
			c(0,10*sin(hr*2*pi/24),10*sin((hr+1)*2*pi/24),0),
			c(0,10*cos(hr*2*pi/24),10*cos((hr+1)*2*pi/24),0),
			col=wedgecol,border=0
		)
	}))
	lines(sin(seq(0,2*pi,length=180)),cos(seq(0,2*pi,length=180)))
	points(predX,predY,pch=21,bg=rainbow(25)[round(trueHr%%24)+1])
}


