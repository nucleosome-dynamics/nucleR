#' Auto detection of a fitted \code{pcKeepComp} param for filterFFT function
#'
#' This function tries to obtain the minimum number of components needed in a
#' FFT filter to achieve or get as close as possible to a given correlation
#' value. Usually you don't need to call directly this function, is used in
#' \code{filterFFT} by default.
#' 
#' This function predicts a suitable \code{pcKeepComp} value for
#' \code{filterFFT} function. This is the recommended amount of components (in
#' percentage) to keep in the \code{filterFFT} function to obtain a correlation
#' of (or near of) \code{cor.target}.
#'
#' The search starts from two given values \code{pc.min, pc.max} and uses
#' linial interpolation to quickly reach a value that gives a corelation
#' between the filtered and the original near \code{cor.target} within the
#' specified tolerance \code{cor.tol}.
#'
#' To allow a quick detection without an exhaustive search, this function uses
#' a subset of the data by randomly sampling those regions with meaningful
#' coverage values (i,e, different from 0 or NA) larger than
#' \code{smpl.min.size}. If it's not possible to obtain \code{smpl.max.size}
#' from this region (this could be due to flanking 0's, for example) at least
#' \code{smpl.min.size} will be used to check correlation. Mean correlation
#' between all sampled regions is used to test the performance of the
#' pcKeepComp parameter.
#'
#' If the number of meaningful bases in \code{data} is less than
#' \code{smpl.min.size * (smpl.num/2)} all the \code{data} vector will be used
#' instead of using sampling.
#'
#' @param data Numeric vector to be filtered
#' @param pc.min,pc.max Range of allowed values for pcKeepComp (minimum and
#' maximum), in the range 0:1.
#' @param max.iter Maximum number of iterations
#' @param verbose Extra information (debug)
#' @param cor.target Target correlation between the filtered and the original
#' profiles. A value around 0.99 is recommeded for Next Generation Sequencing
#' data and around 0.7 for Tiling Arrays.
#' @param cor.tol Tolerance allowed between the obtained correlation an the
#' target one.
#' @param smpl.num If \code{data} is a large vector, some samples from the
#' vector will be used instead the whole dataset. This parameters tells the
#' number of samples to pick.
#' @param smpl.min.size,smpl.max.size Minimum and maximum size of the samples.
#' This is used for selection and sub-selection of ranges with meaningful
#' values (i,e, different from 0 and NA). Power of 2 values are recommended,
#' despite non-mandatory.
#' @param \dots Parameters to be pass to \code{autoPcKeepComp}
#'
#' @return Fitted \code{pcKeepComp} value
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}, David Rosell
#' \email{david.rosell@@irbbarcelona.org}
#' @keywords attribute
#'
#' @examples
#' # Load dataset
#' data(nucleosome_htseq)
#' data <- as.vector(coverage.rpm(nucleosome_htseq)[[1]])
#'
#' # Get recommended pcKeepComp value
#' pckeepcomp <- pcKeepCompDetect(data, cor.target=0.99)
#' print(pckeepcomp)
#'
#' # Call filterFFT
#' f1 <- filterFFT(data, pcKeepComp=pckeepcomp)
#'
#' # Also this can be called directly
#' f2 <- filterFFT(data, pcKeepComp="auto", cor.target=0.99)
#'
#' # Plot
#' plot(data[1:2000], col="black", type="l", lwd=2)
#' lines(f1[1:2000], col="red", lwd=2)
#' lines(f2[1:2000], col="blue", lwd=2, lty=2)
#' legend(
#'     "bottom", c("original", "two calls", "one call"),
#'     col=c("black", "red", "blue"), lty=c(1,1,2), horiz=TRUE, bty="n"
#' )
#'
#' @export pcKeepCompDetect
#'
#' @importFrom stats runif
#' @importFrom IRanges IRanges
#' @importMethodsFrom stats cor
#' @importMethodsFrom IRanges width start
pcKeepCompDetect <- function (data, pc.min = 0.01, pc.max = 0.1, max.iter = 20,
                            verbose = FALSE, cor.target = 0.98, cor.tol = 1e-3,
                            smpl.num = 25, smpl.min.size = 2 ^ 10,
                            smpl.max.size = 2 ^ 14)
{
    # Sample data
    samp <- .sampleData(
        data,
        smpl.num,
        smpl.min.size,
        smpl.max.size,
        verbose=verbose
    )

    # Original params
    pc.min.input <- pc.min
    pc.max.input <- pc.max

    # Calculates the correlation between fft and normal for a given pcKeepComp
    .meancor <- function(pc) {
        if (pc == 1) {
            return(1)
        }

        if (pc == 0) {
            return(0)
        }
        cors <- sapply(samp, function(x) cor(filterFFT(x, pcKeepComp=pc), x))
        return(mean(cors, na.rm=TRUE))
    }

    # Initial evaluation
    eval.min <- .meancor(pc.min)
    eval.max <- .meancor(pc.max)
    desv.min <- eval.min - cor.target
    desv.max <- eval.max - cor.target

    if (abs(desv.min) < abs(desv.max)) {
        pc.best <- pc.min
        dv.best <- desv.min
    } else {
        pc.best <- pc.max
        dv.best <- desv.max
    }

    # Iterate
    i <- 0
    while (
        (abs(dv.best) > cor.tol) &
        (i < max.iter) &
        (pc.best > pc.min.input) &
        (pc.best < pc.max.input)
    ) {
        i <- i + 1
        if (verbose) {
            cat(
                "Iteration",
                i,
                ".  pcKeepComp=",
                pc.best,
                ".  Correlation=",
                dv.best + cor.target,
                "\n"
            )
        }

        # Lineal interpolation of optimal value
        pc.new <- pc.min -
            (desv.min * ((pc.max - pc.min) / (desv.max - desv.min)))
        if (pc.new > 1) {
            pc.new <- 1
        } else if (pc.new < 0) {
            pc.new <- 0
        }

        #New evaluation
        eval.new <- .meancor(pc.new)
        desv.new <- eval.new - cor.target

        #Update limits
        if (desv.min < 0 & desv.max > 0) {
            if (desv.new > 0) {
                pc.max <- pc.new
                desv.max <- desv.new
            } else {
                pc.min <- pc.new
                desv.min <- desv.new
            }

        } else if (desv.min > 0 & desv.max > 0) {
            pc.max <- pc.min
            desv.max <- desv.min
            pc.min <- pc.new
            desv.min <- desv.new

        } else if (desv.min < 0 & desv.max < 0) {
            pc.min <- pc.max
            desv.min <- desv.new
            pc.max <- pc.new
            desv.max <- desv.new
        }

        # New best value update
        if (abs(desv.new) < abs(dv.best)) {
            pc.best <- pc.new
            dv.best <- desv.new
        }
    }

    if (pc.best > pc.max.input) {
        pc.best <- pc.max.input
    }
    if (pc.best < pc.min.input) {
        pc.best <- pc.min.input
    }

    return(pc.best)
}

.sampleData <- function (data, smpl.num, smpl.min.size, smpl.max.size,
                        verbose = FALSE)
{
    # Check coherency
    if (smpl.min.size > smpl.max.size) {
        warning("smpl.min.size > smpl.max.size, using only smpl.min.size")
        smpl.max.size <- smpl.min.size
    }

    res <- list()

    # For short sequences, use all the data
    if (length(data) < smpl.min.size * (smpl.num / 2)) {
        if(verbose) {
            message("Short fragment. Using all data")
        }
        res[[1]] <- data

        # For long sequence, use sampling
    } else {

        if(verbose) {
            message("Long sequence. Trying sampling...")
        }

        # Select ranges <> from 0 or NA and longer than minimum size
        rang <- IRanges(data != 0 & !is.na(data))
        rang <- rang[width(rang) > smpl.min.size]

        tota <- sum(width(rang))
        # If the overall useful bases don't satisfy the criteria, return all
        if (tota < smpl.min.size * (smpl.num / 2)) {
            if (verbose) {
                message(
                    "No enough covered bases to apply sampling. ",
                    "Using all data"
                )
            }
            res[[1]] <- data
        } else {
            if(verbose) {
                message(" Selecting regions for sampling")
            }
            # Sampling with weighted probabilities according range's width
            reps <- sample(
                1:length(rang),
                size=smpl.num,
                replace=TRUE,
                prob=width(rang) / tota
            )

            # Random selection of start points
            # In short, if a range has a width higher than smpl.max.size we
            # will pickup a random start position inside it that still it's in
            # the limits. The width of the final ranges will be the minimum of
            # the range's size and smpl.max.size
            wids <- width(rang[reps])
            marg <- unlist(sapply(wids - smpl.max.size, max, 0))
            rnof <- unlist(sapply(
                marg[marg > 0],
                function(x) floor(runif(n=1, max=x, min=1))
            ))
            marg[marg > 0] <- rnof

            fina <- IRanges(
                start=start(rang[reps]) + marg,
                width=sapply(wids, min, smpl.max.size)
            )

            res <- lapply(fina, function(x) data[x])
        }
    }
    if (verbose) {
        message("Returning ", length(res), " regions")
    }
    return(res)
}
