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
#' Simple function for returning the middle point of a of a GRanges (normal mid
#' doesn't work there)
#' @importMethodsFrom IRanges start end
.mid <- function(x)
    floor((start(x) + end(x)) / 2)

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
        percent <- as.numeric(sub("%","", threshold))
        quantile(data, percent/100, na.rm=TRUE)
    } else {
        threshold
    }
}
