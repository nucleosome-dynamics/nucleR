#################################################################################
#  Average coverage for the given experiments in the given positions 
#################################################################################
coverAroundDF<-function(exps, df, window, avg=TRUE, col.pos="txStart", col.chr="chrom", col.strand="strand", col.row.names="name", cores=1)
{
	if(!(col.pos %in% names(df))) simpleError("Position column not found in reference data frame")
	if(!(col.chr %in% names(df))) simpleError("Chromosome column not found in reference data frame")
	if(!(col.strand %in% names(df))) simpleError("Strand column not found in reference data frame")

	row.names=1:nrow(df)
	if(col.row.names %in% names(df)) row.names = as.character(df[,col.row.names])

	return(coverAround(exps,chrom=as.character(df[,col.chr]), start=as.numeric(df[,col.pos]),
				 end=as.numeric(df[,col.pos]), strand=as.character(df[,col.strand]),
				 names=as.character(df[,col.row.names]), window=window, avg=avg, cores=cores))
}

#Wrapper for SacCer1
coverAroundSC1<-function(exps, window, avg=TRUE, data="/mmb/homes/oflores/workbench/data/SacCer1_genes.RData", feature="txStart", cores=1, names=NULL)
{
	genes = load(data)
	genes = get(genes)

	if(!is.null(names)) genes = genes[genes$name %in% toupper(names),]
	
	return(coverAroundDF(exps, df=genes[genes$chrom != "chrM",], window=window, avg=avg, col.pos=feature, col.chr="chrom", col.strand="strand", col.row.names="name", cores=cores))
}

#CoverAroundIRanges
coverAroundIR<-function(exps, ir, row.names=1:length(ir), window=0, avg=TRUE, cores=1)
{
	return(coverAround(exps, chrom=ir$names, start=ir$start, end=ir$end, strand="+", 
					names=row.names, window=window, avg=avg, cores=cores))
}

#Wrapper for intelligent parallel aplication
coverAround<-function(exps, chrom, start, end, strand="+", names=1:length(start), window=0, avg=TRUE, cores=1)
{
	if(length(exps) > 1 & cores > 1)
	{
		print("Parallel execution mode")
		res = mclapply(exps, .coverAround_1cor, chrom, start, end, strand, names, window, avg, cores=1, mc.cores=cores)
		names(res) = exps
		return(res)
	}
	else
	{
		print("Serial execution mode")
		return(.coverAround_1cor(exps, chrom, start, end, strand, names, window, avg, cores=cores))
	}
}

#####This function does the work
.coverAround_1cor<-function(exps, chrom, start, end, strand="+", names=1:length(start), window=0, avg=TRUE, cores=1)
{
	if(length(strand) == 1) strand = rep(strand, length(start))
	if(length(start) != length(end) | length(start) != length(chrom) | length(start) != length(strand) | length(start) != length(names))
		simpleError("Error: chrom, start, end, and strand must have the same length")

	serial = length(exps) > 1
		
	swap = which(end < start)
	if(length(swap) > 0) { start_tmp = start[swap]; start[swap]=end[swap]; end[swap]=start_tmp }

  win_ini = start - window + 1 #+1 for the "0"
  win_end = end + window -1

	if(serial) { writeLines(paste("  Calculating coverage for", length(exps), "samples: ", sep=" "), sep=" ") }
	else { print(paste("  Calculating coverage for", exps[[1]])) }

  coverage=list()
  for(exp in exps)
  {
    if(serial) writeLines(as.character(length(coverage)+1), sep=" ")

		VEC = cbind(as.character(chrom), win_ini, win_end, strand)
		names(VEC) = NULL
		VEC_L = list(); for(i in 1:nrow(VEC)) VEC_L[[i]] = c(chr=VEC[[i,1]],
										 wstart=VEC[[i,2]], wend=VEC[[i,3]], strand=VEC[[i,4]])
		
		sample=get(exp) #Now it supports Rle lists
		if(length(grep("rle", class(sample), ignore.case=TRUE)) > 0)
		{
			library(IRanges)
			sample = mclapply(as.list(sample), as.numeric, mc.cores=cores)
		}

		max_len = window*2 + max(end-start)
		
		indRange<-function(chr, start, end, strand)
		{
			start = as.numeric(start)
			end = as.numeric(end)
			vec = try(sample[[chr]][start:(end+1)], silent=TRUE)
			if(class(vec) == "try-error") tmp = rep(NA, max_len)
			vec = as.numeric(vec)
			if(strand == "-") vec = rev(vec)
			return(c(vec, rep(0, max_len - length(vec))))
		}		

		tmp = mclapply(VEC_L, function(x) indRange(x["chr"], x["wstart"], x["wend"], x["strand"]), mc.cores=cores)
		coverage[[exp]] = matrix(unlist(tmp), nrow=length(start), ncol=max_len, byrow=TRUE, dimnames=list(names, -window:(max_len - window -1)))

	}

	if(serial) { writeLines("   done.") }
	else { print(paste("   ", exps[[1]], "done")) }

	if(length(coverage)==1)
	{
		if(!avg) return(coverage[[1]])
		return(colMeans(coverage[[1]], na.rm=TRUE))
	}
	else
	{
  	if(!avg) return(coverage)
	  return(lapply(coverage, colMeans, na.rm=TRUE))
	}
}

#################################################################################
#  Plots the coverage overlap version (coverArround results)
#################################################################################
plotCoverOverlap<-function(covers, smooth.pc=NULL, norm=FALSE, cols=NULL, ltys=NULL, lwds=NULL, leg.x="bottomright", ylim=NULL, xlim=NULL, ...)
{
	if(norm)
	{
		for(cov in names(covers))
		{
			tmp = covers[[cov]]
			tmp = tmp - min(tmp)
			tmp = tmp / max(tmp)
			covers[[cov]] = tmp
		}

	}

	if(is.null("ylim"))
	{
		min = min(unlist(lapply(covers, min)))
		max = max(unlist(lapply(covers, max)))
		ylim=c(min-abs(min*0.01), max+abs(max*0.01))
	}

	if(is.null("xlim"))
	{
		labels = as.numeric(names(covers[[1]]))
		xlim=c(min(labels), max(labels))
	}

	if(is.null(cols))
	{
		n = length(covers)
		if(n > 2 && n <= 10) { cols = RColorBrewer::brewer.pal(n, "Paired")}else{ cols = rainbow(n) }
	} 

	if(is.null(ltys)) ltys = rep(1, length(covers))
	if(is.null(lwds)) lwds = rep(2, length(covers))

  if(!is.null(smooth.pc))
  {
    covers = mclapply(covers, filterFFT0, pcKeepComp=smooth.pc)
  }

	plot(as.numeric(names(covers[[1]])), covers[[1]], ylim=ylim, xlim=xlim, type="l", lwd=lwds[[1]], col=cols[[1]], lty=ltys[[1]], ...)

	for(i in 2:length(covers))
	{
		lines(as.numeric(names(covers[[i]])), covers[[i]], type="l", lwd=lwds[[i]], col=cols[[i]], lty=ltys[[i]])
	}

	if(leg.x != "n") legend(leg.x, names(covers), col=cols, lty=ltys, lwd=lwds, cex=0.8, bty="n")
}

#################################################################################
#  Plots the coverage individual (non-avg) version (coverArround results, avg=FALSE)
#################################################################################
plotCoverIndiv<-function(cover, col=NULL, breaks=NULL,  ylab="Genes", xlab="Position from TSS", file="", ...)
{
	if(file != "")
	{
		png(paste(file, "png", sep="."), width=6000, height=6000); 
		par(cex=10)
	}

	if(is.null(breaks)) breaks = c((0:8)/8)
	if(is.null(col)) col = RColorBrewer::brewer.pal(8, "Blues")

	cov = apply(cover, 1, function(x) (x-min(x)) / max((x-min(x))))
	cov = t(cov)
	cov = apply(cov, 1, function(x) rev(x))
	
	image(x=rev(as.numeric(rownames(cov))), y=1:ncol(cov), z=cov, zlim=c(0,1), breaks=breaks, col=col, ylab=ylab, xlab=xlab, ...)

	if(file != "")
	{
		dev.off()
	}
}
