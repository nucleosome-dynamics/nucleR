mergeCalls <- function(calls, dyad.dist=70, discard.low=0.2, mc.cores=1)
{
	res = lapply(calls, .mergeSpace, dyad.dist, discard.low=discard.low, mc.cores=mc.cores)
	return(do.call(c, unname(res)))
}

.mergeSpace <-function(calls, dyad.dist, discard.low, mc.cores)
{
	mc.cores = .check.mc(mc.cores)

	#Set only one state
	calls$nmerge = 1

	#Calculate distace between middles (only upper diagonal)
	dists = sapply(.mid(calls), function(x) abs(.mid(calls) - x))
	dists[lower.tri(dists, diag=TRUE)] = NA

	#Elements closer than the thresdhold
	tomerge = which(dists < dyad.dist, arr.ind=TRUE)

	#This returs the call to merge with the following "n"
	merger = IRanges(1:max(tomerge) %in% tomerge[,1])

	#Create a matrix with [,1] = start call   [,2] = number of following calls to merge
	merger_mat = matrix(c(start(merger), width(merger)), ncol=2, byrow=FALSE)

	#Joins two or more calls given and index of the before matrix	
	.joinCallsMat <- function(idx)
	{
		start = merger_mat[idx,1]
		end = merger_mat[idx,1] + merger_mat[idx,2]
		cr0 = calls[start:end,]
		cr = cr0[cr0$score_h > discard.low,] #Discard low peaks in fuzzy (force separation)
		if(nrow(cr) == 0) cr = cr0 #If all peaks discarded, use all instead
		tmp = data.frame(as.data.frame(reduce(ranges(cr)[[1]])), as.list(colMeans(as.data.frame(values(cr)[[1]]))))
		tmp$nmerge = nrow(cr)
		return(tmp)
	}
	
	#Apply in parallel
	if(mc.cores > 1)
	{
		tmp = mclapply(1:nrow(merger_mat), .joinCallsMat, mc.cores=mc.cores)
	}else{
		tmp = lapply(1:nrow(merger_mat), .joinCallsMat) 
	}
	
	#Join and add the space name
	tmp2 = do.call("rbind", tmp)
	tmp2$space = levels(space(calls))

	#Select WP nucleosomes
	wp = as.data.frame(calls[-as.vector(tomerge),])

	#Return all of them
	return(RangedData(rbind(wp, tmp2)))
}
