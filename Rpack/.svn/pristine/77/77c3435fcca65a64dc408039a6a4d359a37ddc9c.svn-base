setMethod("plotPeaks", signature(peaks="numeric"),
  function(peaks, data, threshold=0, scores=NULL, start=1, end=length(data), 
		 xlab="position", type="l", col.points="red", thr.lty=1, thr.lwd=1, thr.col="darkred",
		 scor.col=col.points, scor.font=2, scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2, ...) {

	#Calculate the ranges in threshold and get the coverage
	#If threshold is given as a string with percentage, convert it
	if(!is.numeric(threshold)) if(grep("%", threshold) == 1)
  {
		threshold = quantile(data, as.numeric(sub("%","", threshold))/100, na.rm=TRUE)
	}

	data = data[start:end]
	names(data) = start:end

	subset = which(peaks >= start & peaks <= end)
	peaks = peaks[subset]
	if(!is.null(scores)) scores = scores[subset]

	plot(start:end, data, type=type, xlab=xlab, ...)
	if(is.null(scores)) points(peaks, data[as.character(peaks)], col=col.points)
	if(threshold != 0) abline(h=threshold, col=thr.col, lty=thr.lty, lwd=thr.lwd)

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
		scoreMatrix = as.data.frame(values(peaks)[[1]])
		if(ncol(scoreMatrix) == 0) scoreMatrix=NULL
		plotPeaks(peaks=ranges(peaks)[[1]], data=data, scores=scoreMatrix, ...)
	}
)

setMethod("plotPeaks", signature(peaks="IRanges"),
  function(peaks, data, threshold=0, scores=NULL, start=1, end=length(data), dyn.pos=TRUE,
		 xlab="position", type="l", col.points="red", thr.lty=1, thr.lwd=1, thr.col="darkred",
		 rect.thick=2, rect.lwd=1, rect.border="black", 
     scor.col=col.points, scor.font=2, scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2, 
		 indiv.scores=TRUE, ...) {

		if(!is.null(scores) & is.numeric(scores)) scores = data.frame(score=scores)

		#Calculate the ranges in threshold and get the coverage
		#If threshold is given as a string with percentage, convert it
	  if(!is.numeric(threshold)) if(grep("%", threshold) == 1)
  	{
	   threshold = quantile(data, as.numeric(sub("%","", threshold))/100, na.rm=TRUE)
	  }

		#Subset data
		subset = which(end(peaks) >= start & start(peaks) <= end)
  	data = data[start:end]
	  names(data) = start:end
	  peaks = peaks[subset]
		if(!is.null(scores)) scores = scores[subset,]

		#One percentile for vertical measures
		pc = (0.01*max(data))

		#Dynamic positioning (on top of the peaks) or all the same
		win_m = round(min(width(peaks)) / 2) / 2 #This takes a widow |half|
		midpoints = .mid(peaks)
		if(dyn.pos) ybottom = sapply(midpoints, function(x) max(data[(x-win_m):(x+win_m)]))
		else 				ybottom = quantile(data, threshold, na.rm=TRUE)

		#Overlap correction
		if(!dyn.pos)
		{
	  	bins = disjointBins(IRanges(start(peaks), end(peaks) + 10))
			ybottom = ((bins-1) * pc * 2) + ybottom
		}

		#Plot coverage
		plot(start:end, data, type=type, xlab=xlab, ...)

		#We have information about merged calls
		if(!is.null(scores) & "nmerge" %in% names(scores))
		{
			#Plot WP nucleosomes
			p1 = peaks[scores$nmerge == 1,]
			ybot = ybottom[scores$nmerge == 1]
			if(length(p1) > 0) 	rect(start(p1), ybot+pc, end(p1), ybot+pc*rect.thick, lwd=rect.lwd,
													 col=col.points, border=rect.border)

			#Plot fuzzy nucleosomes
			p2 = peaks[scores$nmerge > 1,]
			ybot = ybottom[scores$nmerge > 1]
			if(length(p2) > 0) rect(start(p2), ybot+pc, end(p2), ybot+pc*rect.thick, lwd=rect.lwd,
													 col=col.points, border=rect.border, density=30)

		#All are normal calls
		}else{
			rect(start(peaks), ybottom+pc, end(peaks), ybottom+pc*rect.thick, lwd=rect.lwd,
                           col=col.points, border=rect.border)
		}

  	if(threshold != 0) abline(h=threshold, col=thr.col, lty=thr.lty, lwd=thr.lwd)

		#Add text
		if(!is.null(scores))
		{
			.f <- function(x) format(x, digits=scor.digits)

			#Composite or simple score
			if("score_w" %in% names(scores) & "score_h" %in% names(scores) & indiv.scores)
			{
				scores = paste(.f(scores$score), " (", .f(scores$score_h), "h | ", .f(scores$score_w), "w)", sep="")
			}else{
				scores = .f(scores$score)
			}	

			text(midpoints, ybottom+pc*rect.thick+pc*3, labels=scores,
					 cex=scor.cex, adj=scor.adj, col=scor.col, font=scor.font)
		}
	}
)

