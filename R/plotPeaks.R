#' Nucleosome calling plot function
#'
#' Helper function for a quick and convenient overview of nucleosome calling
#' data.
#'
#' This function is intended to plot data previously processed with `nucleR`
#' pipeline. It shows a coverage/intensity profile toghether with the
#' identified peaks. If available, score of each peak is also shown.
#'
#' @param peaks `numeric`, `data.frame`, `IRanges` or `GRanges` object
#'   containing the detected peaks information. See help of [peakDetection()]
#'   or [peakScoring()] for more details.
#' @param data Coverage or Tiling Array intensities
#' @param threshold Threshold applied in `peakDetection`
#' @param scores If `peaks` is a `data.frame` or a `GRanges` it's obtained from
#'  'score' column, otherwise, `scores` can be given here as a `numeric`
#'   vector.
#' @param start,end Start and end points defining a subset in the range of
#'   `data`. This is a convenient way to plot only a small region of data,
#'   without dealing with subsetting of range or score objects.
#' @param dyn.pos If peaks are ranges, should they be positioned dynamicaly on
#'   top of the peaks or staticaly at `threshold` baseline. Spacing of
#'   overlapping ranges is automatically applied if `FALSE`.
#' @param xlab,ylab,type,col.points Default values with general properties of
#'   the plot
#' @param thr.lty,thr.lwd,thr.col Default values with general properties for
#'   threshold representation
#' @param rect.thick,rect.lwd,rect.border Default values for
#'   [ggplot2::geom_rect()] representation of ranges. `rect.thick` indicates
#'   the thickness of the rectangles.
#' @param scor.col,scor.nudge,scor.cex,scor.digits Default values for
#'   [ggplot2::geom_text()] representation for score numbers, if available.
#' @param indiv.scores Show or hide individual scores for width and height in
#'   brakets besides the mixed score.
#' @param \dots Arguments to be passed to other methods.
#'
#' @return (none)
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
#' @seealso [peakDetection()], [peakScoring()], [ggplot2::ggplot()],
#' @keywords hplot
#' @rdname plotPeaks
#'
#' @examples
#'
#' # Generate a random peaks profile
#' reads <- syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#' cover <- coverage.rpm(reads)
#'
#' # Filter them
#' cover_fft <- filterFFT(cover)
#'
#' # Detect peaks
#' peaks <- peakDetection(cover_fft, threshold="40%", score=TRUE, width=140)
#'
#' # Plot peaks and coverage profile (show only a window)
#' plotPeaks(peaks, cover_fft, threshold="40%", start=1000, end=6000)
#'
#' @export
#'
setGeneric(
    "plotPeaks",
    function (peaks, data, ...)
        standardGeneric("plotPeaks")
)

#' @rdname plotPeaks
#' @importFrom ggplot2 ggplot xlim xlab ylab theme
setMethod(
    "plotPeaks",
    signature(peaks="numeric"),
    function (peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
              xlab="position", ylab="coverage", type=1, col.points="red",
              thr.lty=1, thr.lwd=1, thr.col="darkred", scor.col=col.points,
              scor.cex=2.5, scor.digits=2, scor.nudge=2000)
    {
        threshold <- .getThreshold(threshold, data)
        covdf <- .makeDataDf(data, start, end)
        peakdf <- .makePeakDf(peaks, scores, data, start, end)

        ggplot() +
            .plotCov(covdf, linetype=type) +
            .plotScoreOrDots(df=peakdf,
                             col.points=col.points,
                             scor.col=scor.col,
                             scor.digits=scor.digits,
                             scor.nudge=scor.nudge,
                             scor.cex=scor.cex) +
            .plotThresh(threshold, color=thr.col, lty=thr.lty, lwd=thr.lwd) +
            xlim(start, end) +
            xlab(xlab) +
            ylab(ylab) +
            theme(legend.position="none")
    }
)

#' @rdname plotPeaks
setMethod(
    "plotPeaks",
    signature(peaks="data.frame"),
    function (peaks, data, ...)
        plotPeaks(peaks=peaks$peak, data=data, scores=peaks$score, ...)
)

#' @rdname plotPeaks
#' @importMethodsFrom S4Vectors values runLength
#' @importMethodsFrom IRanges ranges
setMethod(
    "plotPeaks",
    signature(peaks="GRanges"),
    function (peaks, data, ...)
    {
        if (sum(runLength(seqnames(peaks)) > 0) > 1) {
            stop("Only GRanges with a single seqname are supported")
        }
        scoreMatrix <- as.data.frame(values(peaks))
        if (ncol(scoreMatrix) == 0) {
            scoreMatrix <- NULL
        }
        plotPeaks(peaks=ranges(peaks), data=data, scores=scoreMatrix, ...)
    }
)

#' @rdname plotPeaks
#' @importFrom IRanges IRanges
#' @importMethodsFrom S4Vectors space values
setMethod(
    "plotPeaks",
    signature(peaks="RangedData"),
    function (peaks, data, ...)
    {
        if (length(unique(space(peaks))) > 1) {
            stop("Only uni-spatial RangedData is supported")
        }
        scoreMatrix <- as.data.frame(values(peaks)[[1]])
        if (ncol(scoreMatrix) == 0) {
            scoreMatrix <- NULL
        }
        plotPeaks(peaks=ranges(peaks)[[1]], data=data, scores=scoreMatrix, ...)
    }
)

#' @rdname plotPeaks
#' @importFrom ggplot2 ggplot scale_alpha_manual xlim xlab ylab theme
setMethod(
    "plotPeaks",
    signature(peaks="IRanges"),
    function (peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
              dyn.pos=TRUE, xlab="position", ylab="coverage", type=1,
              col.points="red", thr.lty=1, thr.lwd=1, thr.col="darkred",
              rect.thick=2, rect.lwd=0.5, rect.border="black",
              scor.col=col.points, scor.cex=2.5, scor.digits=2,
              indiv.scores=FALSE, scor.nudge=2000)
    {
        threshold <- .getThreshold(threshold, data)

        covdf <- .makeDataDf(data, start, end)
        df <- .initPeaksDf(peaks, scores, start, end)

        pc <- .getPc(covdf)
        df[, "bottom"] <- .getYBottom(df,
                                      covdf,
                                      pc,
                                      dyn.pos=dyn.pos,
                                      threshold=threshold)

        df[, "ymin"] <- df[, "bottom"] + pc
        df[, "ymax"] <- df[, "bottom"] + pc * rect.thick
        if ("nmerge" %in% names(df)) {
            df[, "type"] <- ifelse(df[, "nmerge"] > 1, "fuz", "wp")
        }

        ggplot() +
            .plotCov(covdf, linetype=type) +
            .plotPeakRects(df,
                           color = rect.border,
                           fill  = scor.col,
                           size  = rect.lwd) +
            scale_alpha_manual(values=c(wp=1, fuz=0.2)) +
            .plotScores(df,
                        scor.digits  = scor.digits,
                        indiv.scores = indiv.scores,
                        pc           = pc,
                        size         = scor.cex,
                        color        = scor.col,
                        nudge_y      = scor.nudge) +
            .plotThresh(threshold, color=thr.col, lty=thr.lty, lwd=thr.lwd) +
            xlim(start, end) +
            xlab(xlab) +
            ylab(ylab) +
            theme(legend.position="none")
    }
)

.makeDataDf <- function (data, start, end)
{
    covdf <- data.frame(x=seq_along(data), y=data)
    if (!missing(start) && !is.null(start)) {
        covdf <- covdf[covdf[, "x"] > start, ]
    }
    if (!missing(end) && !is.null(end)) {
        covdf <- covdf[covdf[, "x"] < end, ]
    }
    covdf
}

.initPeaksDf <- function (peaks, scores, start, end)
{
    df <- as.data.frame(peaks)
    if (!missing(scores) && !is.null(scores)) {
        df <- cbind(df, scores)
    }
    df$mid <- .mid(peaks)
    if (!missing(start) && !is.null(start)) {
        df <- df[df[, "start"] > start, ]
    }
    if (!missing(end) && !is.null(end)) {
        df <- df[df[, "end"] < end, ]
    }
    df
}

#' @importFrom ggplot2 geom_line aes_string
.plotCov <- function (covdf, ...)
    geom_line(data=covdf, mapping=aes_string(x="x", y="y"), ...)

.getPc <- function (data)
    0.01 * max(data[, "y"])

#' @importFrom stats quantile
#' @importMethodsFrom IRanges disjointBins
.getYBottom <- function (df, data, pc, dyn.pos=TRUE, threshold=0)
{
    if (dyn.pos) {
        win_m <- round(min(df$width) / 2) / 2  # This takes a widow |half|
        return(sapply(
            df$mid,
            function(x)  {
                i <- data[, "x"] %in% ((x-win_m):(x+win_m))
                max(data[i, "y"])
            }
        ))
    } else {
        ran <- IRanges(start=df[, "start"], end=df[, "end"]+10)
        bins <- disjointBins(ran)
        bottom <- quantile(data[, "y"], threshold, na.rm=TRUE)
        return((bins-1) * pc*2 + bottom)
    }
}

#' @importFrom ggplot2 geom_rect aes_string
.plotPeakRects <- function (df, ...)
{
    if ("type" %in% names(df)) {
        geom_rect(data=df,
                  mapping=aes_string(xmin  = "start",
                                     xmax  = "end",
                                     ymin  = "ymin",
                                     ymax  = "ymax",
                                     alpha = "type"),
                  ...)
    } else {
        geom_rect(data=df,
                  mapping=aes_string(xmin = "start",
                                     xmax = "end",
                                     ymin = "ymin",
                                     ymax = "ymax"),
                  ...)
    }
}

.getScoreTxt <- function (df, scor.digits, indiv.scores)
{
    .f <- function(x) format(x, digits=scor.digits)
    info.avail <-  "score_w" %in% names(df) & "score_h" %in% names(df)
    if (info.avail && indiv.scores) {
        paste0(
            .f(df[, "score"]), " (",
            .f(df[, "score_h"]), "h | ",
            .f(df[, "score_w"]), "w)"
        )
    } else {
        .f(df[, "score"])
    }
}

#' @importFrom ggplot2 geom_text geom_blank aes_string
.plotScores <- function (df, scor.digits, indiv.scores, pc, ...)
{
    if ("score" %in% names(df)) {
        df[, "txt"] <- .getScoreTxt(df, scor.digits, indiv.scores)
        df[, "txt_pos"] <- df[, "ymax"] + pc*3
        geom_text(data=df,
                  mapping=aes_string(x="mid", y="txt_pos", label="txt"),
                  ...)
    } else {
        geom_blank()
    }
}

#' @importFrom ggplot2 geom_hline geom_blank
.plotThresh <- function (threshold, ...)
{
    if (threshold != 0) {
        geom_hline(yintercept=threshold, ...)
    } else {
        geom_blank()
    }
}

.makePeakDf <- function (peaks, scores, data, start, end)
{
    peakdf <- data.frame(x=peaks, y=data[peaks])
    if (!is.null(scores)) {
        peakdf[, "scores"] <- scores
    }
    if (!missing(start) && !is.null(start)) {
        peakdf <- peakdf[peakdf[, "x"] > start, ]
    }
    if (!missing(end) && !is.null(end)) {
        peakdf <- peakdf[peakdf[, "x"] < end, ]
    }
    peakdf
}

#' @importFrom ggplot2 geom_text geom_point aes_string
.plotScoreOrDots <- function (df, col.points, scor.col, scor.digits,
                              scor.nudge, scor.cex, ...)
{
    if ("scores" %in% names(df)) {
        df[, "scores"] <- format(df[, "scores"], digits=scor.digits)
        geom_text(data    = df,
                  mapping = aes_string(x="x", y="y", label="scores"),
                  nudge_y = scor.nudge,
                  color   = scor.col,
                  size    = scor.cex,
                  ...)
    } else {
        geom_point(data    = df,
                   mapping = aes_string(x="x", y="y"),
                   shape   = 1,
                   color   = col.points,
                   size    = 3,
                   ...)
    }
}
