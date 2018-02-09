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
#' Simple function for returning the middle point of a RangedData or of a
#' GRanges (normal mid doesn't work there)
#' @importMethodsFrom IRanges start end
.mid <- function(x)
    floor((start(x) + end(x)) / 2)
