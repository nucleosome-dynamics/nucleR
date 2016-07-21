setGeneric(
    "export.bed",
    function(ranges, score=NULL, chrom, name, desc=name, filepath=name,
             splitByChrom=TRUE)
        standardGeneric("export.bed")
)
setGeneric(
    "filterFFT",
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)
        standardGeneric("filterFFT")
)
setGeneric(
    "fragmentLenDetect",
    function(reads, samples=1000, window=5000, min.shift=1, max.shift=100,
             mc.cores=1, as.shift=FALSE)
        standardGeneric("fragmentLenDetect")
)
setGeneric(
    "peakDetection",
    function(data, threshold=0.25, width=1, score=TRUE, min.cov=2, mc.cores=1)
        standardGeneric("peakDetection")
)
setGeneric(
    "peakScoring",
    function(peaks, data, threshold=0.25, ...)
        standardGeneric("peakScoring")
)
setGeneric(
    "plotPeaks",
    function(peaks, data, ...)
        standardGeneric("plotPeaks")
)
setGeneric(
    "processReads",
    function(data, type="single", fragmentLen, trim, ...)
        standardGeneric("processReads")
)
setGeneric(
    "controlCorrection",
    function(exp, ctr, ...)
        standardGeneric("controlCorrection")
)

.check.mc <- function(mc.cores)
{   # Check if there's support foor multicore or use only one
    lib <- 'parallel'
    if (mc.cores > 1 && !lib %in% loadedNamespaces()) {

        warning("'", lib, "' library not available, switching to mc.cores=1")
        1
    } else {
        mc.cores
    }
}

.xlapply <- function(X, FUN, ..., mc.cores = 1)
{   # Wrapper to choose between lapply and mclapply accordingly
    actual.cores <- .check.mc(mc.cores)

    if (actual.cores > 1) {
        mclapply(X=X, FUN=FUN, ...=..., mc.cores=actual.cores)
    } else {
        lapply(X=X, FUN=FUN, ...=...)
    }
}

# Simple function for returning the middle point of a RangedData or of a
# GRanges (normal mid doesn't work there)
.mid <- function(x)
    floor((start(x) + end(x)) / 2)
