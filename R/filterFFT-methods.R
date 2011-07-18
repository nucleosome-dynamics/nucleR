setMethod("filterFFT", signature(data="SimpleRleList"),
	function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, mc.cores=1, ...)  {
	temp = lapply(data, as.vector)
	return(filterFFT(temp, pcKeepComp, showPowerSpec, useOptim, mc.cores, ...))
})

setMethod("filterFFT", signature(data="Rle"),
	function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)  {
	return(filterFFT(as.vector(data), pcKeepComp, showPowerSpec, useOptim, ...))
})

setMethod("filterFFT", signature(data="list"),
  function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, mc.cores=1, ...)  {
	
	if(length(data) > 1 & showPowerSpec)
		stop("showPowerSpec only can be applyied to lists of length = 1")

	#Check if multicore is supported or set to 1
	mc.cores = .check.mc(mc.cores)

	if(mc.cores > 1)
	{
		return(mclapply(data, filterFFT, pcKeepComp, showPowerSpec, useOptim, mc.cores=mc.cores, ...))
	}
	else
	{ 
		return(lapply(data, filterFFT, pcKeepComp, showPowerSpec, useOptim, ...))	
	}
})

#This is the filter itself
.filterFFT <- function(data, pcKeepComp)
{
	data[is.na(data)] = 0
	temp = fft(data)
	keep = round(length(temp) * pcKeepComp)
	temp[keep:(length(temp)-keep)] = 0
	return(Re(fft(temp, inverse=TRUE))/length(temp))
}

setMethod("filterFFT", signature(data="numeric"),
	function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)	{
		
		#Check the pcKeepComp
		if(is.numeric(pcKeepComp))
		{
		 	if(pcKeepComp > 1 | pcKeepComp < 0) stop("numeric pcKeepComp must be in range 0:1")
		}else{
			pcKeepComp = pcKeepCompDetect(data, ...)
		}

		if(showPowerSpec)
		{
			data[is.na(data)] = 0
			upper = min(250000, length(data))
			if(length(data) > 250000)
				warning("warning: only first 250,000bp are represented in the power spectrum")
			temp = fft(data[1:upper])
			keep = round(upper * pcKeepComp)
			plot(Re(temp[2:round(length(temp)/2)]), type="l", xlab="Components", ylab="Power",
					 sub="Selected components threshold marked as red line")
			abline(v=keep, col="red", lwd=1, lty=2)
		}

		#Bypass all optimizations, very much slower
		if(!useOptim) return(.filterFFT(data, pcKeepComp))

		#Partition the data for available values
		ranges = IRanges(!is.na(data))
		ranges = ranges[width(ranges) > 100] # Discard short regions
		defVal = NA

		#Split, instead of NAs, the large 0 rows
		if(width(ranges[1]) == length(data))
		{
			r = Rle(data==0)
			r@values[r@values==TRUE & r@lengths>500] = "SPLIT"	
			ranges = IRanges(as.vector(r) != "SPLIT")
			ranges = ranges[width(ranges) > 100]
			defVal = 0
		}
		
		#Define FFT by regions for avoid large amount of memory and drop in performance
		.fftRegion<-function(data2, pcKeepComp)
		{
			#FFT works best with lengths power of 2
			if(length(data2) > (2^19)+(2^17)-5000) #This is ~650,000
			{
				res1 = .filterFFT(data2[1:(2^19)], pcKeepComp)
				res2 = .fftRegion(data2[((2^19)-4999):length(data2)], pcKeepComp) #Recursive call
				res = c(res1[1:((2^19)-2500)], res2[2501:length(res2)])
			} else {
				#This is a trick, append 0 till nearest power of 2 length
				#It doesn't affect too much the periodicity and increase dramatically the speed
				pad = rep(0, 2^ceiling(log2(length(data2))) - length(data2))
				res = .filterFFT(c(data2, pad), pcKeepComp)
				res = res[1:length(data2)]
			}
			return(res)
		}

		fft_ranges = lapply(ranges, function(x) .fftRegion(data[x], pcKeepComp))
		res = rep(defVal, length(data)) #Create a vector of default values and fill it
		for(i in 1:length(ranges)) res[ranges[[i]]] = fft_ranges[[i]]
		
		#Set to default values the positions that have them in the input
		#(remove strange periodicities from large uncovered regions)
		if(is.na(defVal)) 
		{
			rtmp = IRanges(is.na(data))
			rtmp = rtmp[width(rtmp) > 15]
			for(i in 1:length(rtmp)) res[rtmp[[i]]] = NA
		}
		else if (defVal == 0)
		{
			rtmp = IRanges(data == 0)
			rtmp = rtmp[width(rtmp) > 15]
			for(i in 1:length(rtmp)) res[rtmp[[i]]] = 0
			res[res<0] = 0
		}

		return(res)
})

