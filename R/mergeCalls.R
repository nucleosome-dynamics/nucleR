#' Automatic merging of overlapped nucleosome calls
#'
#' This function joints close nucleosome calls into one larger, fuzzy
#' nucleosome.
#'
#' This functions looks for overlapped calls and join those with more than
#' `min.overlap` bases overlapped. More than two reads can be joined in
#' one single call if all of them are overlapped at least that distance with
#' almost another read in the range.
#'
#' Joining is performed in chain, so if nucleosome call A is close to B and B
#' is close to C, the final call will comprise the range A-B-C. The resulting
#' scores (mixed, width, height) of the final joined call will be the average
#' value of the individual scores.
#'
#' The parameter `discard.low` allows to ignore the small peaks that could be
#' merged with larger ones, originating large calls. In the case that all of
#' the overlapped reads in a given position have `score_h` less than
#' `discard.low`, all of them will be selected instead of deleting that call.
#'
#' @param calls [GenomicRanges::GRanges] with scored and ranged nucleosome
#'   calls from `peakScoring` or `peakDetection(..., score=TRUE)`.
#' @param min.overlap Minimum overlap between two reads for merge them
#' @param discard.low Discard low covered calls (i.e. calls with `score_h
#'   < discard.low` will be discarded)
#' @param mc.cores Number of cores available to parallel data processing.
#' @param verbose Show progress info?
#'
#' @return [GenomicRanges::GRanges] with merged calls and the additional data
#'   column `nmerge`, with the count of how many original ranges are merged in
#'   the resulting range.
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso [peakScoring()]
#' @keywords manip
#'
#' @examples
#' # Generate a synthetic coverage map (assuming reads of 40bp and fragments
#' # of 130)
#' map <- syntheticNucMap(
#'     wp.num=20, fuz.num=20,  nuc.len=40, lin.len=130, rnd.seed=1
#' )
#' cover <- filterFFT(coverage.rpm(map$syn.reads))
#'
#' # Find peaks over FFT filtered coverage
#' calls <- peakDetection(filterFFT(
#'     cover, pcKeepComp=0.02), width=130, score=TRUE
#' )
#'
#' # Merge overlapped calls
#' merged_calls <- mergeCalls(calls)
#'
#' plotPeaks(merged_calls, cover)
#'
#' @export mergeCalls
#' @importMethodsFrom GenomeInfoDb seqnames
#'
mergeCalls <- function (calls, min.overlap=50, discard.low=0.2, mc.cores=1,
                        verbose=TRUE)
{
    res <- .xlapply(
        split(calls, seqnames(calls)),
        .mergeSpace,
        min.overlap = min.overlap,
        discard.low = discard.low,
        verbose     = verbose,
        mc.cores    = mc.cores
    )
    return(do.call(c, unname(res)))
}

#' Space merger
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors queryHits
#' @importFrom plyr ldply
#' @importMethodsFrom GenomeInfoDb seqnames
#' @importMethodsFrom IRanges start width findOverlaps reduce
#' @importMethodsFrom S4Vectors runValue
.mergeSpace <- function (calls, min.overlap, discard.low, verbose)
{
    if (verbose) {
        message("* Starting space: ", runValue(seqnames(calls)))
    }

    if (verbose) {
        message(" - Finding overlapped reads")
    }
    ovlps <- findOverlaps(
        calls,
        minoverlap      = min.overlap,
        type            = "any",
        select          = "all",
        drop.self       = TRUE,
        drop.redundant  = TRUE
    )

    # Select those reads wich are overlapped (by construction with the n+1
    # read)
    hits <- queryHits(ovlps)
    # This is the rownumber of ALL the overlapped reads
    selhits <- unique(sort(c(hits, hits+1)))

    # No overlaped reads
    if (length(selhits) == 0) {
        return (calls)
    }

    # Make a list of the id of those reads wich are overlapped and with
    # how many following reads they are overlapped
    red <- reduce(IRanges(start=hits, width=1))

    # Make list of "grouping" calls
    # So [[1]] = 23, 24   means that the first "merged" call is the overlap
    # of rows (calls) 23 and 24
    if (verbose) {
        message(" - Constructing merge list")
    }

    xs <- mapply(
        function (x, y) seq.int(from=x, length.out=y),
        start(red),
        width(red) + 1,
        SIMPLIFY=FALSE
    )

    if (verbose) {
        message(" - Merging calls")
    }
    # Join function
    .join <- function (xi, dfcalls, discard.low) {
        # This is the heigth selection, to avoid low nucleosomes be merged
        # with big ones
        x <- xi[which(dfcalls[xi, "score_h"] > discard.low)]
        if (length(x) == 0) {
            x <- xi
        }
        return(data.frame(
            start   = min(dfcalls[x, "start"]),
            end     = max(dfcalls[x, "end"]),
            score   = mean(dfcalls[x, "score"]),
            score_w = mean(dfcalls[x, "score_w"]),
            score_h = mean(dfcalls[x, "score_h"]),
            nmerge  = length(x)
        ))
    }
    # calculating it on a data.frame is faster than on a GRanges
    dfcalls <- as.data.frame(calls)
    resdf <- ldply(xs, .join, dfcalls, discard.low)
    resdf$seqnames <- runValue(seqnames(calls))
    fuz <- makeGRangesFromDataFrame(resdf, keep.extra.columns=TRUE)

    # Select WP nucleosomes
    wp <- calls[-selhits, ]
    wp$nmerge <- 1

    # Join and order (the last maybe is not needed, but it is nicer)
    all <- sort(c(wp, fuz))

    # Return all of them
    if (verbose) {
        message(
            " - Done (",
            length(wp),
            " non-overlapped | ",
            length(fuz),
            " merged calls)"
        )
    }
    return (all)
}
