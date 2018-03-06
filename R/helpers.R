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

#' mclapply warapper
#'
#' Wrapper to choose between lapply and mclapply accordingly
#'
#' @importFrom parallel mclapply
.xlapply <- function(X, FUN, ..., mc.cores = 1)
{
    actual.cores <- .check.mc(mc.cores)

    if (actual.cores > 1) {
        mclapply(X=X, FUN=FUN, ...=..., mc.cores=actual.cores)
    } else {
        lapply(X=X, FUN=FUN, ...=...)
    }
}

#' Find midpoints
#'
#' Simple function for returning the middle point of a of a GRanges or a
#' data.frame (normal mid doesn't work there)
setGeneric(
    ".mid",
    function (x)
        setGeneric(".mid")
)

#' @importMethodsFrom BiocGenerics start end
setMethod(
    ".mid",
    signature(x="GRanges"),
    function (x)
        floor((start(x) + end(x)) / 2)
)

setMethod(
    ".mid",
    signature(x="data.frame"),
    function (x)
        floor((x$start + x$end) / 2)
)

.lapplyIRange <- function (x, fun, ...)
    lapply(seq_along(x), function (i) fun(x[i], ...))

setGeneric(
    ".whichChr",
    function (x)
        setGeneric(".whichChr")
)

#' @importMethodsFrom GenomeInfoDb seqnames
#' @importMethodsFrom S4Vectors runValue
setMethod(
    ".whichChr",
    signature(x="GRanges"),
    function (x)
        runValue(seqnames(x))
)

setMethod(
    ".whichChr",
    signature(x="RangedData"),
    function (x)
        names(x)
)

setGeneric(
    ".countRows",
    function (x)
        setGeneric(".countRows")
)

setMethod(
    ".countRows",
    signature(x="GRanges"),
    function (x)
        length(x)
)

setMethod(
    ".countRows",
    signature(x="RangedData"),
    function (x)
        nrow(x)
)

#' Threshold getter
#'
#' If threshold is given as a string with percentage, convert it
#'
#' @param threshold threshold given as an absolute value or as a string
#'   percentage
#' @param data vector with values from which to derive the threshold if it's
#'   relative
#'
#' @return a numeric vector
#'
#' @importFrom stats quantile
#
.getThreshold <- function (threshold, data)
{
    if (!is.numeric(threshold) && grepl("%$", threshold)) {
        percent <- as.numeric(sub("%", "", threshold))
        quantile(data, percent/100, na.rm=TRUE)
    } else {
        threshold
    }
}

#' Unlist an IRanges object into a vector
#'
#' Internal function from the IRanges package, lifted and slightly modified to
#' prevent a NOTE warning about the use of a non-exported function
#'
#' @author  H. Pagès, P. Aboyoun and M. Lawrence
#' @importFrom utils getFromNamespace
#' @importMethodsFrom IRanges pos
#'
.unlist_as_integer <- function (x)
{
    stopifnot(is(x, "Ranges"))
    if (is(x, "Pos")) {
        return(pos(x))
    } else {
        fancy_mseq <- getFromNamespace("fancy_mseq", "S4Vectors")
        return(fancy_mseq(width(x), offset=start(x)-1L))
    }
}

#' Vectorized version of `all`
#'
#' Helper function that behaves as a vectorized version of the function `all`
#'
#' @param ... arbitraty amount of `logical` vectors, expected to have the same
#'   length
#' @return `logical` vector
#'
.vectorizedAll <- function(...)
    Reduce(`&`, list(...))

#' File loader
#
#' Higher order function to import BAM or Bowtie files.
#' Deals with wether type is `single` or `paired` and with the number of input
#' files
#'
#' @importFrom GenomicRanges GRanges GRangesList
.loadFiles <- function (singleLoad, pairedLoad)
{
    function (files, type="paired") {
        if (type == "single") {
            f <- .loadSingleBam
        } else if (type == "paired") {
            f <- .loadPairedBam
        } else {
            stop("type must be `single` or `paired`")
        }

        len <- length(files)
        if (len == 0) {
            GRanges()
        } else if (len == 1) {
            f(files[[1]])
        } else {
            GRangesList(lapply(files, f))
        }
    }
}
