setMethod("peakScoring", signature(peaks="list"),
	function(peaks, data, threshold=0.25, mc.cores=1) {

  #Check if multicore is supported or set to 1
  mc.cores = .check.mc(mc.cores) 
	
	if(mc.cores > 1) {
			res = mclapply(peaks, peakScoring, data=data, threshold=threshold, mc.cores=mc.cores)
		}else{
			res = lapply(peaks, peakScoring, data=data, threshold=threshold)
		}
		return(res) # Return the list directly
	}
)

setMethod("peakScoring", signature(peaks="IRangesList"),
	function(peaks, data, threshold=0.25, dyad.length=38, mc.cores=1) {

	  #Check if multicore is supported or set to 1
	  mc.cores = .check.mc(mc.cores) 

		if(mc.cores > 1) {
			peaks_list = lapply(peaks, identity) #mclapply doesn't work with S4 classes
			res = mclapply(peaks_list, peakScoring, data=data, threshold=threshold, dyad.length=dyad.length, mc.cores=mc.cores)  
		}else{
			res = lapply(peaks, peakScoring, data=data, threshold=threshold, dyad.length=dyad.length)
		}

		#Result should be returned as a single RangedData object
    #Put the correct space in RangedData
	  for(name in names(res)) space(res[[name]]) = name

    #Combine RangedData objects into single one
    temp = res[[1]]
    for(i in 2:length(res)) temp = c(temp, res[[i]])

    return(temp)
	}
)

setMethod("peakScoring", signature(peaks="numeric"),
	function(peaks, data, threshold=0.25) {

		mean = mean(data[!is.na(data) & data > quantile(data, threshold, na.rm=TRUE)], na.rm=TRUE)
		sd   = sd(data[!is.na(data) & data > quantile(data, threshold, na.rm=TRUE)], na.rm=TRUE)
	
		res = pnorm(data[peaks], mean=mean, sd=sd, lower.tail=TRUE)
		return(data.frame(peak=peaks, score=res))
	}
)

setMethod("peakScoring", signature(peaks="IRanges"),
  function(peaks, data, threshold=0.25, dyad.length=38) {

		#Hack for TA, that could have negative values in range calculations
		if(min(data, na.rm=TRUE) < 0) data = data+abs(min(data, na.rm=TRUE))

    mean = mean(data[data > quantile(data, threshold, na.rm=TRUE)], na.rm=TRUE)
    sd   = sd(data[data > quantile(data, threshold, na.rm=TRUE)], na.rm=TRUE)

		#Calculate dyad range
		dyad.middl = start(peaks) + floor(width(peaks)/2)
		dyad.start = dyad.middl - floor((dyad.length/2))
		dyad.end   = dyad.start + (dyad.length-1)
		dyads = IRanges(start=dyad.start, end=dyad.end)

	  sums.range = lapply(peaks, function(x) mean(data[x], na.rm=TRUE))
		sums.dyad  = lapply(dyads, function(x) mean(data[x], na.rm=TRUE))
	
		#Score the heigh of the peak
   	scor.heigh = pnorm(data[dyad.middl], mean=mean, sd=sd, lower.tail=TRUE)
		
		#Score the width (dispersion) of the peak
		scor.width = unlist(sums.dyad) / unlist(sums.range)
		scor.width = scor.width / max(scor.width)
	
		#Final score
		scor.final = scor.heigh * scor.width

    return(RangedData(peaks, score=scor.final))
  }
)
