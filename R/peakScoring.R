#' Peak scoring function
#'
#' Scores peaks detected with function \code{peakDetection} according the
#' height and the sharpness (width) of the peak. This function can be called
#' automatically from \code{peakDetection} if \code{score=TRUE}.
#'
#' This function scores each previously identified peak according its height
#' and sharpness.
#'
#' The height score (\code{score_h}) tells how large is a peak, higher means
#' more coverage or intensity, so better positioned nucleosome. This score is
#' obtained by checking the observed peak value in a Normal distribution with
#' the mean and sd of \code{data}. This value is between 0 and 1.
#'
#' The width score (\code{score_w}) is a mesure of how sharp is a peak. With a
#' NGS coverage in mind, a perfect phased (well-positioned) nucleosome is this
#' that starts and ends exactly in the same place many times. The shape of this
#' ideal peak will be a rectangular shape of the lenght of the read. A wider
#' top of a peak could indicate fuzzyness. The parameter \code{dyad.length}
#' tells how long should be the "flat" region of an ideal peak. The optimum
#' value for this parameter is the lenght of the read in single-ended data or
#' the \code{trim} value of the function \code{processReads}. For Tiling Array,
#' the default value should be fine.
#'
#' This score is obtained calculating the ratio between the mean of the
#' nucleosome scope (the one provided by range in the elements of \code{peaks})
#' and the \code{dyad.length} central bases. This value is normalized between 0
#' and 1.
#'
#' For punctual, single points peaks (provided by \code{numeric} vector or list
#' as \code{peaks} attribute) the score returned is the height score.
#'
#' For range \code{peaks} the weighted sum of the heigth and width scores is
#' used. This is: \code{((score_h * weigth.height) / sum.wei) + ((score_w *
#' weigth.width) / sum.wei)}. Note that you can query for only one score by
#' weting its weight to 1 and the other to 0.
#'
#' @param peaks The identified peaks resulting from \code{peakDetection}. Could
#' be a \code{numeric} vector with the position of the peaks, or a
#' \code{IRanges} object with the extended range of the peak. For both types,
#' list support is implemented as a \code{numeric} list or a \code{IRangesList}
#' @param data Data of nucleosome coverage or intensites.
#' @param threshold The non-default \code{threshold} previously used in
#' \code{peakDetection} function, if applicable. Can be given as a percentage
#' string (i.e., \code{"25\%"} will use the value in the 1st quantile of
#' \code{data}) or as an absolute coverage numeric value (i.e., \code{20} will
#' not look for peaks in regions without less than 20 reads (or reads per
#' million)).
#' @param dyad.length How many bases account in the nucleosome dyad for
#' sharpness description. If working with NGS data, works best with the reads
#' width value for single-ended data or the \code{trim} value given to the
#' \code{processReads} function.
#' @param weight.height,weight.width If the score is a range, the height and
#' the widht score (coverage and fuzzynes) can be defined with different
#' weigths with these parameters. See details.
#' @param mc.cores If input is a \code{list} or \code{IRangeList}, and multiple
#' cores support is available, the maximum number of cores for parallel
#' processing.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return In the case of \code{numeric} input, the value returned is a
#' \code{data.frame} containing a 'peak' and a 'score' column. If the input is
#' a \code{list}, the result will be a \code{list} of \code{data.frame}.
#'
#' If input is a \code{IRanges} or \code{IRangesList}, the result will be a
#' RangedData object with one or multiple spaces respectively and a 3 data
#' column with the mixed, width and heigh score.
#'
#' @author Oscar Flores \email{oflores@@mmb.cpb.ub.es}
#' @seealso \code{\link{peakDetection}}, \code{\link{processReads}},
#' @keywords manip
#' @rdname peakScoring
#'
#' @examples
#' # Generate a synthetic map
#' 
#' # Trimmed length nucleosome map
#' map <- syntheticNucMap(nuc.len=40, lin.len=130)
#' 
#' # Get the information of dyads and the coverage
#' peaks <- c(map$wp.starts, map$fz.starts)
#' cover <- filterFFT(coverage(map$syn.reads))
#' 
#' # Calculate the scores
#' scores <- peakScoring(peaks, cover)
#' plotPeaks(scores$peak, cover, scores=scores$score, start=5000, end=10000)
#'
#' @export
setGeneric(
    "peakScoring",
    function(peaks, data, threshold=0.25, ...)
        standardGeneric("peakScoring")
)

#' @rdname peakScoring
setMethod(
    "peakScoring",
    signature(peaks="list"),
    function (peaks, data, threshold="25%", mc.cores=1)
        # Return the list directly
        .xlapply(
            peaks,
            peakScoring,
            data      = data,
            threshold = threshold,
            mc.cores  = mc.cores
        )
)

#' @rdname peakScoring
setMethod(
    "peakScoring",
    signature(peaks="IRangesList"),
    function (peaks, data, threshold="25%", weight.width=1, weight.height=1,
            dyad.length=38, mc.cores=1) {

        res <- .xlapply(
            names(peaks),
            function(x)
                peakScoring(
                    peaks         = peaks[[x]],
                    data          = data[[x]],
                    threshold     = threshold,
                    dyad.length   = dyad.length,
                    weight.width  = weight.width,
                    weight.height = weight.height
                ),
            mc.cores=mc.cores
        )
        names(res) <- names(peaks)

        # Result should be returned as a single RangedData object
        # Put the correct space in RangedData
        for (name in names(res)) {
            if (nrow(res[[name]])) {
                names(res[[name]]) <- name
            }
        }

        # Combine RangedData objects into single one
        do.call(c, unname(res))
    }
)

#' @rdname peakScoring
setMethod(
    "peakScoring",
    signature(peaks="numeric"),
    function (peaks, data, threshold="25%") {

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshold is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- stats::quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        mean <- mean(data[!is.na(data) & data > threshold], na.rm=TRUE)
        sd <- sd(data[!is.na(data) & data > threshold], na.rm=TRUE)

        res <- stats::pnorm(data[peaks], mean=mean, sd=sd, lower.tail=TRUE)
        return (data.frame(peak=peaks, score=res))
    }
)

#' @rdname peakScoring
setMethod(
    "peakScoring",
    signature(peaks="IRanges"),
    function (peaks, data, threshold="25%", weight.width=1, weight.height=1,
            dyad.length=38) {

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

        # Hack for TA, that could have negative values in range calculations
        if (min(data, na.rm=TRUE) < 0) {
            data <- data + abs(min(data, na.rm=TRUE))
        }

        mean <- mean(data[data > threshold], na.rm=TRUE)
        sd <- sd(data[data > threshold], na.rm=TRUE)

        # Calculate dyad range
        dyad.middl <- IRanges::start(peaks) + floor(IRanges::width(peaks) / 2)
        dyad.start <- dyad.middl - floor((dyad.length / 2))
        dyad.end <- dyad.start + (dyad.length - 1)
        dyads <- IRanges::IRanges(start=dyad.start, end=dyad.end)

        sums.range <- lapply(peaks, function(x) mean(data[x], na.rm=TRUE))
        sums.dyad  <- lapply(dyads, function(x) mean(data[x], na.rm=TRUE))

        # Score the heigh of the peak
        scor.heigh <- stats::pnorm(data[dyad.middl], mean=mean, sd=sd, lower.tail=TRUE)

        # Score the width (dispersion) of the peak
        scor.width <- unlist(sums.dyad) / unlist(sums.range)
        scor.width <- scor.width / max(scor.width)

        # Final score
        sum.wei <- weight.width + weight.height
        scor.final <- ((scor.heigh * weight.height) / sum.wei) +
            ((scor.width * weight.width) / sum.wei)

        # Return everything or just merged score
        return (IRanges::RangedData(
            peaks,
            score   = scor.final,
            score_w = scor.width,
            score_h = scor.heigh
        ))
    }
)
