#' Nucleosome calling plot function
#'
#' Helper function for a quick and convenient overview of nucleosome calling
#' data.
#'
#' This function is intended to plot data previously processed with
#' \code{nucleR} pipeline. It shows a coverage/intensity profile toghether with
#' the identified peaks. If available, score of each peak is also shown.
#'
#' @param peaks \code{numeric}, \code{data.frame}, \code{IRanges} or
#' \code{RangedData} object containing the detected peaks information. See help
#' of \code{\link{peakDetection}} or \code{\link{peakScoring}} for more
#' details.
#' @param data Coverage or Tiling Array intensities
#' @param threshold Threshold applied in \code{peakDetection}
#' @param scores If \code{peaks} is a \code{data.frame} or a \code{RangedData}
#' it's obtained from 'score' column, otherwise, \code{scores} can be given
#' here as a \code{numeric} vector.
#' @param start,end Start and end points defining a subset in the range of
#' \code{data}. This is a convenient way to plot only a small region of data,
#' without dealing with subsetting of range or score objects.
#' @param dyn.pos If peaks are ranges, should they be positioned dynamicaly on
#' top of the peaks or staticaly at \code{threshold} baseline. Spacing of
#' overlapping ranges is automatically applied if \code{FALSE}.
#' @param xlab,type,col.points Default values to be passed to \code{plot} and
#' \code{points}
#' @param thr.lty,thr.lwd,thr.col Default values to be passed to \code{abline}
#' for threshold representation
#' @param rect.thick,rect.lwd,rect.border Default values for \code{rect}
#' representation or ranges. \code{rect.thick} indicates the thickness (in
#' percentage relative to y-axis range) of the rectangles.
#' @param scor.col,scor.font,scor.adj,scor.cex,scor.digits Default values for
#' \code{text} representation for score numbers, if available.
#' @param indiv.scores Show or hide individual scores for width and height in
#' brakets besides the mixed score.
#' @param \dots Other parameters passed to \code{\link{plot}} function
#'
#' @return (none)
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{peakDetection}}, \code{\link{peakScoring}},
#' \code{\link{plot}},
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
    function(peaks, data, ...)
        standardGeneric("plotPeaks")
)

#' @rdname plotPeaks
#' @importFrom graphics plot points abline text
#' @importFrom stats quantile
setMethod(
    "plotPeaks",
    signature(peaks="numeric"),
    function (peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
            xlab="position", type="l", col.points="red", thr.lty=1, thr.lwd=1,
            thr.col="darkred", scor.col=col.points, scor.font=2,
            scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2, ...) {

        # Calculate the ranges in threshold and get the coverage
        # If threshold is given as a string with percentage, convert it
        if (!is.numeric(threshold)) {
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%","", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        data <- data[start:end]
        names(data) <- start:end

        subset <- which(peaks >= start & peaks <= end)
        peaks <- peaks[subset]

        if (!is.null(scores)) {
            scores <- scores[subset]
        }

        plot(start:end, data, type=type, xlab=xlab, ...)
        if (is.null(scores)) {
            points(peaks, data[as.character(peaks)], col=col.points)
        }
        if (threshold != 0) {
            abline(h=threshold, col=thr.col, lty=thr.lty, lwd=thr.lwd)
        }

        if (!is.null(scores))
            text(
                peaks, data[as.character(peaks)],
                labels=format(scores, digits=scor.digits),
                cex=scor.cex, adj=scor.adj, col=scor.col,
                font=scor.font
            )
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
#' @importFrom IRanges IRanges
#' @importMethodsFrom S4Vectors space values
setMethod(
    "plotPeaks",
    signature(peaks="RangedData"),
    function(peaks, data, ...) {
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
#' @importFrom graphics plot rect abline text
#' @importFrom stats quantile
#' @importMethodsFrom IRanges start end width disjointBins
setMethod(
    "plotPeaks",
    signature(peaks="IRanges"),
    function (peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
            dyn.pos=TRUE, xlab="position", type="l", col.points="red",
            thr.lty=1, thr.lwd=1, thr.col="darkred", rect.thick=2,
            rect.lwd=1, rect.border="black", scor.col=col.points,
            scor.font=2, scor.adj=c(0.5,0), scor.cex=0.75, scor.digits=2,
            indiv.scores=TRUE, ...) {

        if (!is.null(scores) & is.numeric(scores)) {
            scores <- data.frame(score=scores)
        }

        # Calculate the ranges in threshold and get the coverage
        # If threshold is given as a string with percentage, convert it
        if (!is.numeric(threshold)) {
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%","", threshold)) / 100,
                    na.rm=TRUE)
            }
        }

        # Subset data
        subset <- which(end(peaks) >= start & start(peaks) <= end)
        data_full <- data
        data <- data[start:end]
        names(data) <- start:end
        peaks <- peaks[subset]

        if (!is.null(scores)) {
            scores <- scores[subset,]
        }

        # One percentile for vertical measures
        pc <- (0.01 * max(data))

        # Dynamic positioning (on top of the peaks) or all the same
        win_m <- round(min(width(peaks)) / 2) / 2  # This takes a widow |half|
        midpoints <- .mid(peaks)

        if (dyn.pos) {
            ybottom <- sapply(
                midpoints,
                function(x)
                    max(data_full[(x - win_m):(x + win_m)])
            )
        } else {
            ybottom <- quantile(data, threshold, na.rm=TRUE)
        }

        # Overlap correction
        if (!dyn.pos) {
            bins <- disjointBins(IRanges(start(peaks), end(peaks) + 10))
            ybottom <- ((bins - 1) * pc * 2) + ybottom
        }

        # Plot coverage
        plot(start:end, data, type=type, xlab=xlab, ...)

        # We have information about merged calls
        if (!is.null(scores) & "nmerge" %in% names(scores)) {
            # Plot WP nucleosomes
            p1 <- peaks[scores$nmerge == 1, ]
            ybot <- ybottom[scores$nmerge == 1]

            if (length(p1) > 0) {
                rect(
                    start(p1),
                    ybot + pc,
                    end(p1),
                    ybot + pc * rect.thick,
                    lwd=rect.lwd,
                    col=col.points,
                    border=rect.border
                )
            }

            # Plot fuzzy nucleosomes
            p2 <- peaks[scores$nmerge > 1, ]
            ybot <- ybottom[scores$nmerge > 1]
            if (length(p2) > 0) {
                rect(
                    start(p2),
                    ybot + pc,
                    end(p2),
                    ybot + pc * rect.thick,
                    lwd=rect.lwd,
                    col=col.points,
                    border=rect.border,
                    density=30
                )
            }

        } else {  # All are normal calls
            rect(
                start(peaks),
                ybottom + pc,
                end(peaks),
                ybottom + pc * rect.thick,
                lwd=rect.lwd,
                col=col.points,
                border=rect.border
            )
        }

        if (threshold != 0) {
            abline(h=threshold, col=thr.col, lty=thr.lty, lwd=thr.lwd)
        }

        # Add text
        if (!is.null(scores)) {
            .f <- function(x) format(x, digits=scor.digits)

            # Composite or simple score
            if ("score_w" %in% names(scores) &
                "score_h" %in% names(scores) &
                indiv.scores) {
                scores <- paste(
                    .f(scores$score), " (",
                    .f(scores$score_h), "h | ",
                    .f(scores$score_w), "w)",
                    sep=""
                )
            } else {
                scores <- .f(scores$score)
            }

            text(
                midpoints,
                ybottom + pc * rect.thick + pc * 3,
                labels=scores,
                cex=scor.cex,
                adj=scor.adj,
                col=scor.col,
                font=scor.font
            )
        }
    }
)
