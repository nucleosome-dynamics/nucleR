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

.xlapply <- function(X, FUN, ..., mc.cores = 1)
{   # Wrapper to choose between lapply and mclapply accordingly
    actual.cores <- .check.mc(mc.cores)

    if (actual.cores > 1) {
        mclapply(X=X, FUN=FUN, ...=..., mc.cores=actual.cores)
    } else {
        lapply(X=X, FUN=FUN, ...=...)
    }
}

# Simple function for returning the middle point of a RangedData or of a
# GRanges (normal mid doesn't work there)
.mid <- function(x)
    floor((start(x) + end(x)) / 2)
