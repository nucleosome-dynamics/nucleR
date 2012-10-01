mergeCalls <- function(calls, min.overlap=50, discard.low=0.2, mc.cores=1, verbose=TRUE)
{
	res = lapply(calls, .mergeSpace, min.overlap, discard.low=discard.low, mc.cores=mc.cores, verbose=verbose)
	return(do.call(c, unname(res)))
}

.mergeSpace <-function(calls, min.overlap, discard.low, mc.cores, verbose)
{
	if(verbose) message("* Starting space: ", names(calls))

	mc.cores = .check.mc(mc.cores)

	#Set only one state
	calls$nmerge = 1

	#Efficient call to find overlaps
	if(verbose) message(" - Finding overlapped reads")
	ovlps = findOverlaps(calls, minoverlap=min.overlap, type="any",
											 select="all", ignoreSelf=TRUE, ignoreRedundant=TRUE)

	#Select those reads wich are overlapped (by construction with the n+1 read)
	hits = queryHits(ovlps[[1]])
	selhits = unique(sort(c(hits, hits+1))) #This is the rownumber of ALL the overlapped reads

	#No overlaped reads
	if(length(selhits) == 0) return(calls)

  #Make a list of the id of those reads wich are overlapped and with
  #how many following reads they are overlapped
  red = reduce(IRanges(start=hits, width=1))

	#Make list of "grouping" calls
	#So [[1]] = 23, 24   means that the first "merged" call is the overlap of rows (calls) 23 and 24
	if(verbose) message(" - Constructing merge list")
	xs = mclapply(1:length(red), function(i) seq.int(from=start(red[i]), length.out=width(red[i])+1),
								mc.cores=mc.cores)

	#This saves a lot of time later, just create vectors
	dfcalls = as.data.frame(calls)
	df_start   = dfcalls$start
	df_end     = dfcalls$end
	df_score   = dfcalls$score
	df_score_w = dfcalls$score_w
	df_score_h = dfcalls$score_h

	if(verbose) message(" - Merging calls")
	#Join function
	.join <- function(xi)
	{
			#This is the heigth selection, to avoid low nucleosomes be merged with big ones
			x = xi[which(df_score_h[xi] > discard.low)]
			if(length(x) == 0) x = xi

			start   = min(df_start[x])
			end     = max(df_end[x])
			score   = mean(df_score[x])
			score_w = mean(df_score_w[x])
			score_h = mean(df_score_h[x])
			nmerge = length(x)
			return(c(start, end, score, score_w, score_h, nmerge))
	}

  #Apply in parallel by precalculated groups
  if(mc.cores > 1)
  {
    res = mclapply(xs, .join, mc.cores=mc.cores)
  }else{
    res = lapply(xs, .join)
  }

	if(verbose) message (" - Formatting results")
	#Join the results as a dataframe
	resdf = data.frame(matrix(unlist(res), ncol=6, byrow=TRUE))
  names(resdf) = c("start", "end", "score", "score_w", "score_h", "nmerge")

	#Add other derivate data and order according what expects RangedData(...) call
	resdf$space = names(ovlps)
	resdf$width = resdf$end - resdf$start
	order = c("space", "start", "end", "width", "score", "score_w", "score_h", "nmerge")
	fuz = resdf[,order]
	
	#Select WP nucleosomes
	wp = as.data.frame(calls[-selhits,])

	#Join and order (the last maybe is not needed, but is nice)
	all = rbind(wp, fuz)
	all = all[order(all$start),]

	#Return all of them
	if(verbose) message(" - Done (", nrow(wp), " non-overlapped | ", nrow(fuz), " merged calls)")
	return(RangedData(all))
}
