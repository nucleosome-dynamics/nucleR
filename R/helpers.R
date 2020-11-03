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
#' @importMethodsFrom BiocGenerics start end
.mid <- function(x)
    floor((start(x) + end(x)) / 2)

.lapplyIRange <- function (x, fun, ...)
    lapply(seq_along(x), function (i) fun(x[i], ...))

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
#' Wrapper to internal function from the IRanges package. Avoids use of
#' \code{:::} and thus prevents a NOTE warning about the use of a non-exported
#' function
#'
#' @author  H. Pages, P. Aboyoun and M. Lawrence
#' @importFrom utils getFromNamespace
#' @importMethodsFrom IRanges pos
#'
.unlist_as_integer <- function (x)
    getFromNamespace("unlist_as_integer", "IRanges")(x)

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
