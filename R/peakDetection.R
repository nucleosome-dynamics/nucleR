#' Detect peaks (local maximum) from values series
#'
#' This function allows a efficient recognition of the local maximums (peaks)
#' in a given numeric vector.
#'
#' It's recommended to smooth the input with \code{filterFFT} prior the
#' detection.
#'
#' @param data Input numeric values, or a list of them
#' @param threshold Threshold value from which the peaks will be selected. Can
#' be given as a percentage string (i.e., \code{"25\%"} will use the value in
#' the 1st quantile of \code{data}) or as an absolute coverage numeric value
#' (i.e., \code{20} will not look for peaks in regions without less than 20
#' reads (or reads per milion)).
#' @param width If a positive integer > 1 is given, the peaks are returned as a
#' range of the given width centered in the local maximum. Useful for
#' nucleosome calling from a coverage peak in the dyad.
#' @param score If TRUE, the results will be scored using \code{peakScoring}
#' function.
#' @param min.cov Minimum coverage that a peak needs in order to be considered
#' as a nucleosome call.
#' @param mc.cores The number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. Parallelization requires at least two
#' cores.
#'
#' @return The type of the return depends on the input parameters:
#'
#' \code{numeric} (or a list of them) if \code{width==1 & score==FALSE}
#' containing the position of the peaks.
#' 
#' \code{data.frame} (or list of them) if \code{width==1 & score==TRUE}
#' containing a 'peak' column with the position of the peak plus a 'score'
#' column with its score.
#' 
#' \code{IRanges} (or \code{IRangesList}) if \code{width>1 & score==FALSE}
#' containing the ranges of the peaks.
#' 
#' \code{RangedData} if \code{width>1 & score==TRUE} containing the ranges of
#' the peaks and the assigned score.
#' @note If \code{width} > 1, those ranges outside the range
#' \code{1:length(data)} will be skipped.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{filterFFT}}, \code{\link{peakScoring}}
#' @keywords manip
#' @rdname peakDetection
#'
#' @examples
#' # Generate a random peaks profile
#' reads <- syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#' cover <- coverage(reads)
#'
#' # Filter them
#' cover_fft <- filterFFT(cover)
#'
#' # Detect and plot peaks (up a bit the threshold for accounting synthetic
#' # data)
#' peaks <- peakDetection(cover_fft, threshold="40%", score=TRUE)
#' plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#'
#' #Now use ranges version, which accounts for fuzziness when scoring
#' peaks <- peakDetection(cover_fft, threshold="40%", score=TRUE, width=147)
#' plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#'
#' @export
#'
setGeneric(
    "peakDetection",
    function(data, threshold=0.25, width=1, score=TRUE, min.cov=2, mc.cores=1)
        standardGeneric("peakDetection")
)

#' @rdname peakDetection
setMethod(
    "peakDetection",
    signature(data="list"),
    function (data, threshold="25%", width=1, score=TRUE, min.cov=2,
            mc.cores=1) {

        res <- .xlapply(
            data,
            peakDetection,
            threshold = threshold,
            width     = width,
            score     = score,
            min.cov   = min.cov,
            mc.cores  = mc.cores
        )

        # Process the result, case with ranges
        if (width > 1) {
            if (score) {
                # res is a list of IRanges, conversion is direct
                res <- res[!sapply(res, is.null)]
                for (name in names(res)) {
                    names(res[[name]]) <- name
                }
                # Combine RangedData objects into single one
                res <- do.call(c, unname(res))

            } else {
                # res is a list of RangedData
                res <- unlist(res)
                if (length(res)) {
                    res <- IRanges::IRangesList(res)
                } else {
                    res <- IRanges::IRangesList()
                }
            }
        }

        return (res)
    }
)

#' @rdname peakDetection
setMethod(
    "peakDetection",
    signature(data="numeric"),
    function (data, threshold="25%", width=1, score=TRUE, min.cov=2,
            mc.cores=1) {

        if (width < 1) {
            stop("'width' attribute should be greater than 1")
        }

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshdol is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- stats::quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        ranges <- IRanges::IRanges(!is.na(data) & data > threshold)
        if (length(ranges) == 0) {
            return(NULL)
        }

        covers <- lapply(ranges, function(x) data[x])

        # For each range, look for changes of trend and keep the starting
        # position of trend change
        pea <- .xlapply(
            covers,
            function(x) {
                if (length(x) == 1) {
                    1
                } else {
                    IRanges::start(IRanges::IRanges(
                        x[2:length(x)] <
                        x[1:(length(x) - 1)]
                    ))
                }
            },
            mc.cores=mc.cores
        )

        # Some peaks can have only one trend, correct them
        unitrend <- which(sapply(pea, function(x) length(x) == 0))
        pea[unitrend] <- sapply(
            covers[unitrend],
            function(x) which(x == max(x))
        )[1]

        # Add start offset to peaks relative to the start of the range
        starts <- IRanges::start(ranges)
        res <- unlist(sapply(
            1:length(starts),
            function(i) pea[[i]] + starts[[i]]
        ))

        res <- res[res <= length(data)]
        # the FFT coverage at the peak should be bigger than a given number
        # (by default, 2)
        res <- res[data[res] > min.cov]

        # Extension
        if (width > 1) {
            ext <- floor(width / 2)
            starts <- res - ext
            # Odd/pair correction
            ends <- res + ifelse(width %% 2 == 0, ext - 1, ext)
            res <- IRanges::IRanges(start=starts, end=ends)
            # Remove out of bounds
            res <- res[IRanges::start(res) > 1 & IRanges::end(res) < length(data)]
        }

        if (score) {
            return (peakScoring(peaks=res, data=data, threshold=threshold))
        } else {
            return (res)
        }
    }
)
