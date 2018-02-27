#' Detect peaks (local maximum) from values series
#'
#' This function allows a efficient recognition of the local maximums (peaks)
#' in a given numeric vector.
#'
#' It's recommended to smooth the input with `filterFFT` prior the detection.
#'
#' @param data Input numeric values, or a list of them
#' @param chromosome Optionally specify the name of the chromosome for input
#'   data that doesn't specify it.
#' @param threshold Threshold value from which the peaks will be selected. Can
#'   be given as a percentage string (i.e., `"25\\%"` will use the value in the
#'   1st quantile of `data`) or as an absolute coverage numeric value (i.e.,
#'   `20` will not look for peaks in regions without less than 20 reads (or
#'   reads per milion)).
#' @param width If a positive integer > 1 is given, the peaks are returned as a
#'   range of the given width centered in the local maximum. Useful for
#'   nucleosome calling from a coverage peak in the dyad.
#' @param score If TRUE, the results will be scored using [peakScoring()]
#'   function.
#' @param min.cov Minimum coverage that a peak needs in order to be considered
#'   as a nucleosome call.
#' @param mc.cores The number of cores to use, i.e. at most how many child
#'   processes will be run simultaneously. Parallelization requires at least
#'   two cores.
#'
#' @return The type of the return depends on the input parameters:
#'
#'   * `numeric` (or a list of them) if `width==1 & score==FALSE` containing
#'     the position of the peaks.
#'   * `data.frame` (or list of them) if `width==1 & score==TRUE` containing a
#'     'peak' column with the position of the peak plus a 'score' column with
#'     its score.
#'   * `IRanges` (or `IRangesList`) if `width>1 & score==FALSE` containing the
#'     ranges of the peaks.
#'   * `GRanges` if `width>1 & score==TRUE` containing the ranges of the
#'      peaks and the assigned score.
#'
#' @note If `width` > 1, those ranges outside the range `1:length(data)` will
#'   be skipped.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso [filterFFT()], [peakScoring()]
#' @keywords manip
#' @rdname peakDetection
#'
#' @examples
#' # Generate a random peaks profile
#' reads <- syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#' cover <- coverage.rpm(reads)
#'
#' # Filter them
#' cover_fft <- filterFFT(cover)
#'
#' # Detect and plot peaks (up a bit the threshold for accounting synthetic
#' # data)
#' peaks <- peakDetection(cover_fft, threshold="40%", score=TRUE)
#' plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#'
#' # Now use ranges version, which accounts for fuzziness when scoring
#' peaks <- peakDetection(cover_fft, threshold="40%", score=TRUE, width=147)
#' plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#'
#' @export
#'
setGeneric(
    "peakDetection",
    function(data, threshold=0.25, chromosome=NULL, width=1, score=TRUE,
             min.cov=2, mc.cores=1)
        standardGeneric("peakDetection")
)

#' @rdname peakDetection
#' @importFrom IRanges IRangesList
#' @importMethodsFrom GenomeInfoDb "seqlevels<-"
setMethod(
    "peakDetection",
    signature(data="list"),
    function (data, threshold="25%", width=1, score=TRUE, min.cov=2,
            mc.cores=1) {

        res <- .xlapply(
            names(data),
            function (x)
                peakDetection(data[[x]],
                              chromosome = x,
                              threshold  = threshold,
                              width      = width,
                              score      = score,
                              min.cov    = min.cov,
                              mc.cores   = mc.cores)
        )
        names(res) <- names(data)

        # Process the result, case with ranges
        if (width > 1) {
            if (score) {
                # res is a list of GRanges
                res <- res[!sapply(res, is.null)]
                for (name in names(res)) {
                    seqlevels(res[[name]]) <- names(res)
                }
                res <- do.call(`c`, unname(res))

            } else {
                # res is a list of IRanges
                res <- unlist(res)
                if (length(res)) {
                    res <- IRangesList(res)
                } else {
                    res <- IRangesList()
                }
            }
        }

        return (res)
    }
)

#' @rdname peakDetection
#' @importFrom IRanges IRanges
#' @importFrom stats quantile
#' @importMethodsFrom BiocGenerics start end
setMethod(
    "peakDetection",
    signature(data="numeric"),
    function (data, threshold="25%", chromosome=NULL, width=1, score=TRUE,
              min.cov=2, mc.cores=1) {

        if (width < 1) {
            stop("'width' attribute should be greater than 1")
        }

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshdol is given as a string with percentage, convert it
            if (grepl("%$", threshold)) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        ranges <- IRanges(!is.na(data) & data > threshold)
        if (length(ranges) == 0) {
            return(NULL)
        }

        covers <- .lapplyIRange(
            ranges,
            function (x) data[.unlist_as_integer(x)]
        )

        # For each range, look for changes of trend and keep the starting
        # position of trend change
        pea <- .xlapply(
            covers,
            function(x) {
                if (length(x) == 1) {
                    1
                } else {
                    start(IRanges(
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
        starts <- start(ranges)
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
            res <- IRanges(start=starts, end=ends)
            # Remove out of bounds
            res <- res[start(res) > 1 & end(res) < length(data)]
        }

        if (score) {
            return (peakScoring(
                peaks      = res,
                chromosome = chromosome,
                data       = data,
                threshold  = threshold
            ))
        } else {
            return (res)
        }
    }
)
