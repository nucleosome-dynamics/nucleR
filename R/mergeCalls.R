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
#' @rdname mergeCalls
#' @export mergeCalls
#'
setGeneric(
    "mergeCalls",
    function (calls, min.overlap=50, discard.low=0.2, mc.cores=1, verbose=TRUE)
        standardGeneric("mergeCalls")
)

#' @rdname mergeCalls
#' @importMethodsFrom GenomeInfoDb seqnames
setMethod(
    "mergeCalls",
    signature(calls="GRanges"),
    function (calls, min.overlap=50, discard.low=0.2, mc.cores=1, verbose=TRUE)
    {
       res <- .xlapply(
           split(calls, seqnames(calls)),
           .mergeChrom,
           min.overlap = min.overlap,
           discard.low = discard.low,
           verbose     = verbose,
           mc.cores    = mc.cores
       )
       return(do.call(c, unname(res)))
   }
)

#' @rdname mergeCalls
setMethod(
    "mergeCalls",
    signature(calls="RangedData"),
    function (calls, min.overlap=50, discard.low=0.2, mc.cores=1, verbose=TRUE)
    {
       res <- .xlapply(
           calls,
           .mergeChrom,
           min.overlap = min.overlap,
           discard.low = discard.low,
           verbose     = verbose,
           mc.cores    = mc.cores
       )
       return(do.call(c, unname(res)))
   }
)

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom plyr ddply
#' @importMethodsFrom BiocGenerics as.data.frame
.mergeChrom <- function (calls, min.overlap, discard.low, verbose)
{
    if (verbose) {
        message("* Starting space: ", .whichChr(calls))
        message(" - Finding overlapped reads")
    }

    # we'll use data.frames to merge them more easily; remove unused
    # columns so we can forget about them
    merge.group <- .getMergeGroup(calls, min.overlap)
    dfcalls <- BiocGenerics::as.data.frame(calls)
    dfcalls[, "strand"] <- NULL
    dfcalls[, "width"] <- NULL
    dfcalls[, "merge.group"] <- merge.group

    if (verbose) {
        message(" - Merging calls")
    }
    resdf <- ddply(dfcalls, "merge.group", .joinNucs, discard.low)
    resdf[, "merge.group"] <- NULL  # we don't need this column anymore
    res <- makeGRangesFromDataFrame(resdf, keep.extra.columns=TRUE)

    if (verbose) {
        message(
            " - Done (",
            sum(merge.group == 0),
            " non-overlapped | ",
            sum(merge.group != 0),
            " merged calls)"
        )
    }

    return (sort(res))
}

#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importMethodsFrom IRanges start end findOverlaps reduce
.getMergeGroup <- function (calls, min.overlap)
{
    ovlps <- findOverlaps(
        calls,
        minoverlap     = min.overlap,
        type           = "any",
        select         = "all",
        drop.self      = TRUE,
        drop.redundant = TRUE
    )

    # Select those reads wich are overlapped (by construction with the n+1
    # read)
    hits <- queryHits(ovlps)

    # Make a list of the id of those reads wich are overlapped and with
    # how many following reads they are overlapped
    red <- reduce(IRanges(start=hits, width=1))

    # we give a value of 0 to nucleosomes that are not to be merged and a
    # positive number to the ones that will be merged, giving the same
    # number to the ones that are to be merged together
    merge.group <- rep(0, .countRows(calls))
    for (i in seq_along(red)) {
        from <- start(red)[i]
        to <- end(red)[i] + 1
        merge.group[from:to] <- i
    }

    return (merge.group)
}

.joinNucs <- function (x, discard.low)
{
    seqnames <- x[1, "seqnames"]
    space <- x[1, "space"]
    if (!is.null(seqnames)) {
        chr <- seqnames
    } else if (!is.null(space)) {  # in case we started with a RangedData
        chr <- space
        x[, "seqnames"] <- chr
        x[, "space"] <- NULL
    }

    if (x[1, "merge.group"] == 0) {  # nucleosmes that won't be merged
        x[, "nmerge"] <- 1
        return(x)
    } else {  # nucleosome that will be merged
        x <- x[x[, "score_h"] > discard.low, ]
        return(data.frame(
            seqnames = x[1, "seqnames"],
            start    = min(x[, "start"]),
            end      = max(x[, "end"]),
            score    = mean(x[, "score"]),
            score_w  = mean(x[, "score_w"]),
            score_h  = mean(x[, "score_h"]),
            nmerge   = nrow(x)
        ))
    }
}
