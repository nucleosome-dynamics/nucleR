setMethod("peakDetection", signature(data="list"),
	function(data, threshold=0.25, width=1, score=TRUE, mc.cores=1){

  #Check if multicore is supported or set to 1
  mc.cores = .check.mc(mc.cores) 
	
		if(mc.cores > 1){
			res = mclapply(data, peakDetection, threshold=threshold,
									 width=width, score=score, mc.cores=mc.cores)
		}
		else
			res = lapply(data, peakDetection, threshold=threshold, width=width, score=score)

		#Process the result, case with ranges
		if(width > 1){
		
			#res is a list of IRanges, conversion is direct
			if(score==FALSE){ 
				return(IRangesList(unlist(res)))

			#res is a list of RangedData
			}else{
				#Put the correct space in RangedData by changing range names
				for(name in names(res)) names(res[[name]]@ranges) = name

				#Combine RangedData objects into single one
				if(length(res) > 1)
				{
					temp = res[[1]]
					for(i in 2:length(res)) temp = c(temp, res[[i]])
				}
				else temp = res[[1]]

				return(temp)
			}
		}

		#Case without ranges, anyway return a list
		else
			return(res)
	}
)

setMethod("peakDetection", signature(data="numeric"),
	function(data, threshold=0.25, width=1, score=TRUE, mc.cores=1) {

	if(width < 1) stop("'width' attribute should be greater than 1")

  #Check if multicore is supported or set to 1
  mc.cores = .check.mc(mc.cores) 

	#Calculate the ranges in threshold and get the coverage
	ranges = IRanges(!is.na(data) & data > quantile(data, threshold, na.rm=TRUE))
	covers = lapply(ranges, function(x) data[x])

	#For each range, look for changes of trend and keep the starting position of trend change
	if(mc.cores > 1)
	{
		pea = mclapply(covers, function(x) if(length(x)==1){return(1)}else{return(start(IRanges(x[2:length(x)] < x[1:(length(x)-1)])))}, mc.cores=mc.cores)
	}	else {
		pea = lapply(covers, function(x) if(length(x)==1){return(1)}else{return(start(IRanges(x[2:length(x)] < x[1:(length(x)-1)])))})
	}

	#Some peaks can have only one trend, correct them
	unitrend = which(sapply(pea, function(x) length(x) == 0))
	pea[unitrend] = sapply(covers[unitrend], function(x) which(x==max(x)))[1]

	#Add start offset to peaks relative to the start of the range
	starts = start(ranges)
	res = unlist(sapply(1:length(starts), function(i) pea[[i]] + starts[[i]]))

	#Extension
	if(width > 1)
	{
		ext = floor(width/2)
		starts = res - ext
		ends 	 = res + ifelse(width %% 2 == 0, ext-1, ext) #Odd/pair correction
		res = IRanges(start=starts, end=ends)
		res = res[start(res) > 1 & end(res) < length(data)] #Remove out of bounds
	}

	if(score) return(peakScoring(peaks=res, data=data, threshold=threshold))
	else return(res)
}
)

