#' Peak scoring function
#'
#' Scores peaks detected with function `peakDetection` according the height and
#' the sharpness (width) of the peak. This function can be called automatically
#' from `peakDetection` if `score=TRUE`.
#'
#' This function scores each previously identified peak according its height
#' and sharpness.
#'
#' The height score (`score_h`) tells how large is a peak, higher means more
#' coverage or intensity, so better positioned nucleosome. This score is
#' obtained by checking the observed peak value in a Normal distribution with
#' the mean and sd of `data`. This value is between 0 and 1.
#'
#' The width score (`score_w`) is a mesure of how sharp is a peak. With a NGS
#' coverage in mind, a perfect phased (well-positioned) nucleosome is this that
#' starts and ends exactly in the same place many times. The shape of this
#' ideal peak will be a rectangular shape of the lenght of the read. A wider
#' top of a peak could indicate fuzzyness. The parameter `dyad.length` tells
#' how long should be the "flat" region of an ideal peak. The optimum value for
#' this parameter is the lenght of the read in single-ended data or the `trim`
#' value of the function `processReads`. For Tiling Arrays, the default value
#' should be fine.
#'
#' This score is obtained calculating the ratio between the mean of the
#' nucleosome scope (the one provided by range in the elements of `peaks`) and
#' the `dyad.length` central bases. This value is normalized between 0 and 1.
#'
#' For punctual, single points peaks (provided by `numeric` vector or list as
#' `peaks` attribute) the score returned is the height score.
#'
#' For range `peaks` the weighted sum of the heigth and width scores is used.
#' This is: `((score_h * weigth.height) / sum.wei) + ((score_w * weigth.width)
#' / sum.wei)`. Note that you can query for only one score by setting its
#' weight to 1 and the other to 0.
#'
#' @param peaks The identified peaks resulting from `peakDetection`. Could be a
#'   `numeric` vector with the position of the peaks, or a `IRanges` object
#'   with the extended range of the peak. For both types, list support is
#'   implemented as a `numeric` list or a `IRangesList`
#' @param data Data of nucleosome coverage or intensites.
#' @param chromosome Optionally specify the name of the chromosome for input
#'   data that doesn't specify it.
#' @param threshold The non-default `threshold` previously used in
#'   `peakDetection` function, if applicable. Can be given as a percentage
#'   string (i.e., `"25\\%"` will use the value in the 1st quantile of `data`)
#'   or as an absolute coverage numeric value (i.e., `20` will not look for
#'   peaks in regions without less than 20 reads (or reads per million)).
#' @param dyad.length How many bases account in the nucleosome dyad for
#'   sharpness description. If working with NGS data, works best with the reads
#'   width value for single-ended data or the `trim` value given to the
#'   `processReads` function.
#' @param weight.height,weight.width If the score is a range, the height and
#'   the widht score (coverage and fuzzynes) can be defined with different
#'   weigths with these parameters. See details.
#' @param mc.cores If input is a `list` or `IRangeList`, and multiple cores
#'   support is available, the maximum number of cores for parallel processing.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return In the case of `numeric` input, the value returned is a `data.frame`
#' containing a 'peak' and a 'score' column. If the input is a `list`, the
#' result will be a `list` of `data.frame`.
#'
#' If input is a `IRanges` or `IRangesList`, the result will be a `data.frame`
#' or [GenomicRanges::GRanges] object with one or multiple spaces respectively
#' and a 3 data column with the mixed, width and heigh score.
#'
#' @author Oscar Flores \email{oflores@@mmb.cpb.ub.es}
#' @seealso [peakDetection()], [processReads()],
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
#' cover <- filterFFT(coverage.rpm(map$syn.reads))
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
#' @importMethodsFrom GenomeInfoDb "seqlevels<-"
setMethod(
    "peakScoring",
    signature(peaks="IRangesList"),
    function (peaks, data, threshold="25%", weight.width=1, weight.height=1,
              dyad.length=38, mc.cores=1) {
        chrs <- names(peaks)
        res <- .xlapply(
            chrs,
            function(x)
                peakScoring(
                    peaks         = peaks[[x]],
                    data          = data[[x]],
                    chromosome    = x,
                    threshold     = threshold,
                    dyad.length   = dyad.length,
                    weight.width  = weight.width,
                    weight.height = weight.height
                ),
            mc.cores=mc.cores
        )
        names(res) <- chrs
        for (i in chrs) {
            seqlevels(res[[i]]) <- chrs
        }
        do.call(`c`, unname(res))
    }
)

#' @rdname peakScoring
#' @importFrom stats pnorm quantile
setMethod(
    "peakScoring",
    signature(peaks="numeric"),
    function (peaks, data, chromosome=NULL, threshold="25%") {

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshold is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        mean <- mean(data[!is.na(data) & data > threshold], na.rm=TRUE)
        sd <- sd(data[!is.na(data) & data > threshold], na.rm=TRUE)

        res <- pnorm(data[peaks], mean=mean, sd=sd, lower.tail=TRUE)
        if (!is.null(chromosome)) {
            return (data.frame(peak=peaks, score=res, chromosome=chromosome))
        } else {
            return (data.frame(peak=peaks, score=res))
        }
    }
)

#' @rdname peakScoring
#' @importFrom stats pnorm
#' @importFrom IRanges IRanges
#' @importFrom stats quantile
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom IRanges start end width
setMethod(
    "peakScoring",
    signature(peaks="IRanges"),
    function (peaks, data, chromosome=NULL, threshold="25%", weight.width=1,
              weight.height=1, dyad.length=38) {

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshdol is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
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
        dyad.middl <- start(peaks) + floor(width(peaks) / 2)
        dyad.start <- dyad.middl - floor((dyad.length / 2))
        dyad.end <- dyad.start + (dyad.length - 1)
        dyads <- IRanges(start=dyad.start, end=dyad.end)

        sums.range <- .lapplyIRange(
            peaks,
            function (x)
                mean(data[.unlist_as_integer(x)], na.rm=TRUE)
        )
        sums.dyad <- .lapplyIRange(
            dyads,
            function (x)
                mean(data[.unlist_as_integer(x)], na.rm=TRUE)
        )

        # Score the heigh of the peak
        scor.heigh <- pnorm(
            data[dyad.middl], mean=mean, sd=sd, lower.tail=TRUE
        )

        # Score the width (dispersion) of the peak
        scor.width <- unlist(sums.dyad) / unlist(sums.range)
        scor.width <- scor.width / max(scor.width)

        # Final score
        sum.wei <- weight.width + weight.height
        scor.final <- ((scor.heigh * weight.height) / sum.wei) +
            ((scor.width * weight.width) / sum.wei)

        if (is.null(chromosome)) {
            chromosome <- "*"
        }

        return (GRanges(
            seqnames = chromosome,
            ranges   = peaks,
            score    = scor.final,
            score_w  = scor.width,
            score_h  = scor.heigh
        ))
    }
)
