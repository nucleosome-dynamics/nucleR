#' Automatic merging of overlapped nucleosome calls
#'
#' This function joints close nucleosome calls into one larger, fuzzy
#' nucleosome.
#'
#' This functions looks for overlapped calls and join those with more than
#' \code{min.overlap} bases overlapped. More than two reads can be joined in
#' one single call if all of them are overlapped at least that distance with
#' almost another read in the range.
#'
#' Joining is performed in chain, so if nucleosome call A is close to B and B
#' is close to C, the final call will comprise the range A-B-C. The resulting
#' scores (mixed, width, height) of the final joined call will be the average
#' value of the individual scores.
#'
#' The parameter \code{discard.low} allows to ignore the small peaks that could
#' be merged with larger ones, originating large calls. In the case that all of
#' the overlapped reads in a given position have \code{score_h} less than
#' \code{discard.low}, all of them will be selected instead of deleting that
#' call.
#'
#' @param calls \code{RangedData} with scored and ranged nucleosome calls from
#' \code{peakScoring} or \code{peakDetection(..., score=TRUE)}.
#' @param min.overlap Minimum overlap between two reads for merge them
#' @param discard.low Discard low covered calls (i.e. calls with \code{score_h
#' < discard.low} will be discarded)
#' @param mc.cores Number of cores available to parallel data processing.
#' @param verbose Show progress info?
#' @return \code{RangedData} with merged calls and the additional data column
#' \code{nmerge}, with the count of how many original ranges are merged in the
#' resulting range.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{peakScoring}}
#' @keywords manip
#'
#' @examples
#' # Generate a synthetic coverage map (assuming reads of 40bp and fragments
#' # of 130)
#' map <- syntheticNucMap(
#'     wp.num=20, fuz.num=20,  nuc.len=40, lin.len=130, rnd.seed=1
#' )
#' cover <- filterFFT(coverage(map$syn.reads))
#'
#' # Find peaks over FFT filtered coverage
#' calls <- peakDetection(filterFFT(
#'     cover, pcKeepComp=0.02), width=130, score=TRUE
#' )
#'
#' # Merge overlapped calls
#' merged_calls = mergeCalls(calls)
#'
#' plotPeaks(merged_calls, cover)
#'
#' @export mergeCalls
#'
mergeCalls <- function (calls, min.overlap = 50, discard.low = 0.2,
                        mc.cores = 1, verbose = TRUE)
{
    res <- lapply(
        calls,
        .mergeSpace,
        min.overlap = min.overlap,
        discard.low = discard.low,
        mc.cores    = mc.cores,
        verbose     = verbose
    )
    return(do.call(c, unname(res)))
}

.mergeSpace <- function(calls, min.overlap, discard.low, mc.cores, verbose)
{
    if (verbose) {
        message("* Starting space: ", names(calls))
    }

    mc.cores <- .check.mc(mc.cores)

    # Set only one state
    calls$nmerge <- 1

    # Efficient call to find overlaps
    if (verbose) {
        message(" - Finding overlapped reads")
    }
    ovlps <- IRanges::findOverlaps(
        calls,
        minoverlap = min.overlap,
        type       = "any",
        select     = "all"
    )

    # Select those reads wich are overlapped (by construction with the n+1 read)
    hits <- S4Vectors::queryHits(ovlps[[1]])
    # This is the rownumber of ALL the overlapped reads
    selhits <- unique(sort(c(hits, hits + 1)))

    #No overlaped reads
    if (length(selhits) == 0) {
        return(calls)
    }

    # Make a list of the id of those reads wich are overlapped and with
    # how many following reads they are overlapped
    red <- IRanges::reduce(IRanges::IRanges(start=hits, width=1))

    # Make list of "grouping" calls
    # So [[1]] = 23, 24   means that the first "merged" call is the overlap
    # of rows (calls) 23 and 24
    if (verbose) {
        message(" - Constructing merge list")
    }

    xs <- mapply(
        function (x, y) seq.int(from=x, length.out=y),
        IRanges::start(red),
        IRanges::width(red) + 1,
        SIMPLIFY=FALSE
    )

    # This saves a lot of time later, just create vectors
    dfcalls <- as.data.frame(calls)
    df_start <- dfcalls$start
    df_end <- dfcalls$end
    df_score <- dfcalls$score
    df_score_w <- dfcalls$score_w
    df_score_h <- dfcalls$score_h

    if (verbose) {
        message(" - Merging calls")
    }
    # Join function
    .join <- function(xi) {
        # This is the heigth selection, to avoid low nucleosomes be merged
        # with big ones
        x <- xi[which(df_score_h[xi] > discard.low)]
        if (length(x) == 0) {
            x <- xi
        }

        start <- min(df_start[x])
        end <- max(df_end[x])
        score <- mean(df_score[x])
        score_w <- mean(df_score_w[x])
        score_h <- mean(df_score_h[x])
        nmerge <- length(x)
        return (c(start, end, score, score_w, score_h, nmerge))
    }

    res <- .xlapply(xs, .join, mc.cores=mc.cores)

    if (verbose) {
        message (" - Formatting results")
    }
    # Join the results as a dataframe
    resdf <- data.frame(matrix(unlist(res), ncol=6, byrow=TRUE))
    names(resdf) <- c("start", "end", "score", "score_w", "score_h", "nmerge")

    # Add other derivate data and order according what expects
    # RangedData(...) call
    resdf$space <- names(ovlps)
    resdf$width <- resdf$end - resdf$start
    order <- c(
        "space",
        "start",
        "end",
        "width",
        "score",
        "score_w",
        "score_h",
        "nmerge"
    )
    fuz <- resdf[, order]

    # Select WP nucleosomes
    wp <- as.data.frame(calls[-selhits, ])

    # Join and order (the last maybe is not needed, but is nice)
    all <- rbind(wp, fuz)
    all <- all[order(all$start), ]

    # Return all of them
    if (verbose) {
        message(
            " - Done (",
            nrow(wp),
            " non-overlapped | ",
            nrow(fuz),
            " merged calls)"
        )
    }
    return(IRanges::RangedData(all))
}
