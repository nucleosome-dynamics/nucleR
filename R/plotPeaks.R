setMethod("plotPeaks", signature(peaks="numeric"),
  function(peaks, data, threshold=0, scores=NULL, start=1, end=length(data), 
		 xlab="position", type="l", col.points="red", thr.lty=1, thr.lwd=1, thr.col="darkred",
		 scor.col=col.points, scor.font=2, scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2, ...) {

	data = data[start:end]
	names(data) = start:end

	subset = which(peaks >= start & peaks <= end)
	peaks = peaks[subset]
	scores = scores[subset]

	plot(start:end, data, type=type, xlab=xlab, ...)
	if(is.null(scores)) points(peaks, data[as.character(peaks)], col=col.points)
	if(threshold != 0) abline(h=quantile(data, threshold, na.rm=TRUE), col=thr.col, lty=thr.lty, lwd=thr.lwd)

	if(!is.null(scores))
		 text(peaks, data[as.character(peaks)], labels=format(scores, digits=scor.digits), 
		 cex=scor.cex, adj=scor.adj, col=scor.col, font=scor.font)
	}
)

setMethod("plotPeaks", signature(peaks="data.frame"),
  function(peaks, data, ...) {
		plotPeaks(peaks=peaks$peak, data=data, scores=peaks$score, ...)
  }
)

setMethod("plotPeaks", signature(peaks="RangedData"),
	function(peaks, data, ...){
		if(length(unique(space(peaks))) > 1) stop("Only uni-spatial RangedData is supported")
		plotPeaks(peaks=ranges(peaks)[[1]], data=data, scores=peaks$score, ...)
	}
)

setMethod("plotPeaks", signature(peaks="IRanges"),
  function(peaks, data, threshold=0, scores=NULL, start=1, end=length(data), dyn.pos=TRUE,
		 xlab="position", type="l", col.points="red", thr.lty=1, thr.lwd=1, thr.col="darkred",
		 rect.thick=2, rect.lwd=1, rect.border="black", 
     scor.col=col.points, scor.font=2, scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2, ...) {

		#Subset data
		subset = which(start(peaks) >= start & end(peaks) <= end)
  	data = data[start:end]
	  names(data) = start:end
	  peaks = peaks[subset]
  	scores = scores[subset]

		#One percentile for vertical measures
		pc = (0.01*max(data))

		#Dynamic positioning (on top of the peaks) or all the same
		midpoints = floor((start(peaks) + end(peaks)) / 2)
		if(dyn.pos) ybottom = data[as.character(midpoints)] + pc
		else 				ybottom = quantile(data, threshold, na.rm=TRUE)

		#Overlap correction
		if(!dyn.pos)
		{
	  	bins = disjointBins(IRanges(start(peaks), end(peaks) + 10))
			ybottom = ((bins-1) * pc * 2) + ybottom
		}

		#Plot lines and ranges
		plot(start:end, data, type=type, xlab=xlab, ...)
		rect(start(peaks), ybottom+pc, end(peaks), ybottom+pc*rect.thick, lwd=rect.lwd, col=col.points, border=rect.border)
  	if(threshold != 0) abline(h=quantile(data, threshold, na.rm=TRUE), col=thr.col, lty=thr.lty, lwd=thr.lwd)

		#Add text
		if(!is.null(scores))
		{
			text(midpoints, ybottom+pc*rect.thick+pc*3, labels=format(scores, digits=scor.digits),
					 cex=scor.cex, adj=scor.adj, col=scor.col, font=scor.font)
		}
	}
)

