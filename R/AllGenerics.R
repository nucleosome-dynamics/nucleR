setGeneric("export.bed", function(ranges, score=NULL, chrom, name, description, filepath=name, splitByChrom=TRUE)  standardGeneric("export.bed"))
setGeneric("filterFFT", function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...) standardGeneric("filterFFT"))
setGeneric("fragmentLenDetect", function(reads, samples=1000, window=5000, min.shift=1, max.shift=200, mc.cores=1) standardGeneric("fragmentLenDetect"))
setGeneric("peakDetection", function(data, threshold=0.25, width=1, score=TRUE, mc.cores=1) standardGeneric("peakDetection"))
setGeneric("peakScoring", function(peaks, data, threshold=0.25, ...) standardGeneric("peakScoring"))
setGeneric("plotPeaks", function(peaks, data, ...) standardGeneric("plotPeaks"))
setGeneric("processReads", function(data, type="single", fragmentLen, trim) standardGeneric("processReads"))
setGeneric("controlCorrection", function(exp, ctr, ...) standardGeneric("controlCorrection"))

#Check if there's support for multicore or use only one
.check.mc<-function(mc.cores)
{
  if(mc.cores > 1)
  {
    succ.mc = 'multicore' %in% loadedNamespaces()
    if(!succ.mc) succ = library("multicore", logical.return=TRUE)
    if(!succ.mc)
    {
      warning("'multicore' library not available, switching to mc.cores=1")
      mc.cores=1
    }
  }

	return(mc.cores)
}

