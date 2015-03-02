################################################################################################
# Internal function for detection plotting
################################################################################################
.doPlot <- function(res, calls, cover, window, plot.limy)
{
  chr = as.character(res$chrom)
	pos = res$pos
	id  = as.character(res$id)

	#Plot peaks and calls
  plotPeaks(data=as.numeric(cover[[chr]]), peaks=calls[chr], xlim=pos + c(-window,window), main=id, ylim=c(0, plot.limy))
  abline(v=pos, col="red", lwd=2)

	#Add +1 nucleosome
  if(!is.na(res$p1.pos))
  {
    text(res$p1.pos, plot.limy/4*3, paste("+1 (", res$p1.class, ")", sep=""), cex=1, font=2)
    abline(v=res$p1.pos, col="darkgreen", lty=3)
  }

	#Add -1 nucleosome
  if(!is.na(res$m1.pos))
  {
    text(res$m1.pos, plot.limy/4*3, paste("-1 (", res$m1.class, ")", sep=""), cex=1, font=2)
    abline(v=res$m1.pos, col="darkorchid", lty=3)
  }

  #Add arrow for strand
	if(!is.na(res$p1.pos) & !is.na(res$m1.pos))
	{
	  arrows(res$m1.pos, plot.limy/4*3.3, res$p1.pos, plot.limy/4*3.3, lwd=1, length=0.10, col="black")
	}else{
		str = ifelse(res$strand=="+", 1, -1)
		arrows(pos+(-25*str), plot.limy/4*3.3, pos-(-25*str), plot.limy/4*3.3, lwd=2, length=0.10, col="darkgrey")
	}

	#Add distance
  if(!is.na(res$p1.pos) & !is.na(res$m1.pos))
    text((res$m1.pos+res$p1.pos)/2, plot.limy/4*3.35, paste("dist =", res$dist), cex=0.8, font=3, adj=c(0.5,0))

	#Add message
  if(res$msg != "") text(pos, plot.limy/8*7, adj=c(0.5,0), substr(x=res$msg, start=4, stop=nchar(res$msg)), font=4)
}

################################################################################################
#Multiple execution using data.frame
################################################################################################
nucleosomePatternsDF <- function(calls, cover=NULL, df, col.id="name", col.chrom="chrom", col.pos="pos",
												col.strand="strand", show.plot=FALSE, mc.cores=1, ...)
{
	#Check if nmerge is in the column
	if(!"nmerge" %in% colnames(calls))
	{
		warning("Column 'nmerge' not found in 'calls'. Running mergeCalls with default parameters")
		calls = mergeCalls(calls)
	}

	#Cannot plot if mc.cores>1, parallel output not supported
	if(mc.cores > 1 & show.plot==TRUE) stop("'show.plot=TRUE' cannot work with 'mc.cores > 1'")	

	#We need the coverage for plotting
	if(is.null(cover) & show.plot) stop("'cover' should be provided if 'show.plot=TRUE'")

	#mc.cores = .check.mc(mc.cores) #TODO: Uncomment after testing

	#Parallel execution		
	if(mc.cores > 1)
	{
		res = mclapply(1:nrow(df), function(x) nucleosomePatterns(calls=calls, cover=cover, id=df[x, col.id],
																							chrom=df[x, col.chrom], pos=df[x, col.pos], strand=df[x, col.strand],...),
																							 mc.cores=mc.cores)

	#Serial execution
	}else{
		res = lapply(1:nrow(df), function(x) nucleosomePatterns(calls=calls, cover=cover, id=df[x, col.id],
      chrom=df[x, col.chrom], pos=df[x, col.pos], strand=df[x, col.strand], show.plot=show.plot, ...))
	}


	#Join results in a single data.frame
	res = do.call("rbind", res)

	return(res)
}


################################################################################################
# Pattern detection in a single locus
################################################################################################
nucleosomePatterns <- function(calls, cover=NULL, id, chrom, pos, strand="+", window=300, show.plot=FALSE,
											plot.limy=20, p1.max.merge=3, p1.max.downstream=20, open.thresh=215,
											wp.thresh.score_w=0.6, wp.thresh.score_h=0.4, max.uncovered=150)
{
	#Error if asking for plot without coverage
	if(is.null(cover) & show.plot) stop("'cover' should be provided if 'show.plot=TRUE'")

  #Check if nmerge is in the data.frame
  if(!"nmerge" %in% colnames(calls)) 
  {
    warning("Column 'nmerge' not found in 'calls'. Running mergeCalls with default parameters")
    calls = mergeCalls(calls)
  }

  #The default information, we will fill it in the next steps
  res = data.frame(id=id, chrom=chrom, strand=strand, pos=pos, p1.pos=NA, p1.score=NA, p1.score_w=NA, 
									 p1.score_h=NA, p1.class=NA, m1.pos=NA, m1.score=NA, m1.score_w=NA, m1.score_h=NA,
									 m1.class=NA, p1.nmerge=NA, m1.nmerge=NA, dist=NA, dist.class=NA, descr=NA, msg="")
	chrom = as.character(chrom)
	strand = as.character(strand)

  #Select nucleosomes nearby the given position
  start_in_bounds = unlist(.mid(calls) > (pos - window), use.names=FALSE)
  end_in_bounds = unlist(.mid(calls) < (pos + window), use.names=FALSE)
	sel = start_in_bounds & end_in_bounds & space(calls) == chrom

	#Uncovered regions detection
	if(!is.null(cover))
	{
		if((pos-window) < 0 | (pos+window) > length(cover[[chrom]]))
		{
        res$msg = paste(res$msg, "WINDOW OUT OF CHROMOSOME BOUNDS", sep=" | ")
        if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
        return(res)

		}else{
			if(length(which(cover[[chrom]][(pos-window):(pos+window)] == 0)) > max.uncovered)
			{
			  res$msg = paste(res$msg, "UNCOVERED REGION", sep=" | ")
	  	  if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
				return(res)
			}
		}
	}

	#Check if there's any nucleosome around
  if(!any(which(sel)))
  {
    res$msg = paste(res$msg, "NO NEARBY NUC.", sep=" | ")
    if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
    return(res)
  }

	nearby = calls[sel,]

  #Detect +1 nucleosome, strand sensitive
	if(strand=="+")
	{
		p1 = nearby[.mid(nearby) > (pos - p1.max.downstream),]
		if(nrow(p1) > 1) p1 = p1[.mid(p1) == min(.mid(p1)),]

	}else{ #- strand
    p1 = nearby[.mid(nearby) < (pos + p1.max.downstream),]
    if(nrow(p1) > 1) p1 = p1[.mid(p1) == max(.mid(p1)),]
	}

	#No +1 nucleosome upstream ref. point
	if(nrow(p1) < 1)
  {
		res$msg = paste(res$msg, "+1 NUC MISSING", sep=" | ")
		res$descr = "+1_missing"
		if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
		return(res)
	}

  res$p1.pos = .mid(p1)
  res$p1.score = score(p1)
	res$p1.score_h = p1$score_h
	res$p1.score_w = p1$score_w
	res$p1.nmerge = p1$nmerge

	#Check if +1 is a low complexity zone
	if(p1$nmerge > 3)
	{
		res$p1.class = "F+"
		res$descr = "+1_too_fuzzy"
		res$msg = paste(res$msg, "+1 NUC. TOO FUZZY", sep=" | ")
    if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
    return(res)		
	}

	#Detect -1 nucleosome, strand sensitive
  if(strand == "+")
  {
    m1 = nearby[start(p1) - .mid(nearby) > 0,]
		if(nrow(m1) > 0)
		{
	    m1 = m1[start(m1) == max(start(m1)),]
#  	  dist = start(p1) - end(m1)
			dist = .mid(p1) - .mid(m1)
		}
  }else{
    m1 = nearby[end(p1) - .mid(nearby) < 0,]
		if(nrow(m1) > 0)
		{
    	m1 = m1[start(m1) == min(start(m1)),]
#    	dist = start(m1) - end(p1)
			dist = .mid(p1) - .mid(m1)
		}
  }

	#No -1 nucleosome is detected
  if(nrow(m1) < 1)
  {
    res$msg = paste(res$msg, "-1 NUC MISSING", sep=" | ")
		res$descr="-1_missing"
    res$m1.class="M"
 #  if(show.plot) .doPlot(res, calls, cover, window, plot.limy)
 #  return(res)

  }else{

	  res$m1.pos = .mid(m1)
  	res$m1.score = score(m1)
		res$m1.score_h = m1$score_h
		res$m1.score_w = m1$score_w
		res$m1.nmerge = m1$nmerge

	  res$dist = abs(dist)
	}

	#Evaluate class and add description
  if(is.na(res$p1.class))
		res$p1.class = ifelse(p1$score_w > wp.thresh.score_w & p1$score_h > wp.thresh.score_h, "W", "F")
  if(is.na(res$m1.class))
		res$m1.class = ifelse(m1$score_w > wp.thresh.score_w & m1$score_h > wp.thresh.score_h, "W", "F")

	.gimmeDist <- function(x){
		 if(is.na(x)) return("-")
		 if(x > open.thresh) return("open")
;		 if(x < 120) return("overlap")
		 return("close")
	}

  res$dist.class = sapply(res$dist, .gimmeDist)
	res$descr = paste(res$m1.class, res$dist.class, res$p1.class, sep="-")
  res$descr[is.na(res$dist.class)] = NA

  if(show.plot) .doPlot(res, calls, cover, window, plot.limy)

	#Remove " | " from the message begining (should be done after plotting)
	if(res$msg != "") res$msg = substr(res$msg, 4, nchar(res$msg))

  return(res)
}

plotPatterns <- function(coverage, classif, window=300, ...)
{
	res = classif #alias for c&p
	win = window

	res2 = res[!is.na(res$descr),] #rm na's

	#Get unique classes and assign colors
	classes = unique(res2$descr)
	colors = rainbow(length(classes))
	names(colors) = classes

	#Plot each class
	for(class in classes)
	{

	  #Get rows from current class and center nuc. +1 
	  tmp = res2[res2$descr == class,]
  	tmp$p1.pos.zer = tmp$p1.pos + (tmp$pos - tmp$p1.pos)

	  #Calculate coverage, average and individual
	  c1 = coverAroundDF(coverage, df=tmp, window=win, avg=TRUE, col.pos="p1.pos.zer", col.names="id")
	  c1_ind = coverAroundDF(coverage, df=tmp, window=win, avg=FALSE, col.pos="p1.pos.zer", col.names="id")

	  #Main plot
	  plot(-win:win, c1, type="l", main=class, sub=paste(nrow(tmp), "genes"), col=colors[class], ...)

	  #Sd calculation
	  sds = apply(c1_ind, 2, sd)
	  c1_p_sd = c1 + sds
	  c1_m_sd = c1 - sds

	  #Sd lines
	  lines(-win:win, c1_p_sd, lty=2, col="grey", lwd=2)
	  lines(-win:win, c1_m_sd, lty=2, col="grey", lwd=2)
	  abline(v=0, col="red")

	  legend("bottomright", c("coverage", "+/- std dev"), lty=c(1,2), col=c(colors[class], "grey"), bty="n")
	}
}
