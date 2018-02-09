#' Clean noise and smoothing for genomic data using Fourier-analysis
#'
#' Remove noise from genomic data smoothing and cleaning the observed signal.
#' This function doesn't alter the shape or the values of the signal as much as
#' the traditional method of sliding window average does, providing a great
#' correlation within the original and filtered data (>0.99).
#'
#' Fourier-analysis principal components selection is widely used in signal
#' processing theory for an unbiased cleaning of a signal over the time.
#' 
#' Other procedures, as the traditional sliding window average, can change too
#' much the shape of the results in function of the size of the window, and
#' moreover they don't only smooth the noise without removing it.
#'
#' With a Fourier Transform of the original signal, the input signal is
#' descomposed in diferent wavelets and described as a combination of them.
#' Long frequencies can be explained as a function of two ore more periodical
#' shorter frequecies. This is the reason why long, unperiodic sequences are
#' usually identified as noise, and therefore is desireable to remove them from
#' the signal we have to process.
#'
#' This procedure here is applied to genomic data, providing a novel method to
#' obtain perfectly clean values wich allow an efficient detection of the peaks
#' which can be used for a direct nucleosome position recognition.
#'
#' This function select a certain number of components in the original power
#' spectrum (the result of the Fast Fourier Transform which can be seen with
#' \code{showPowerSpec=TRUE}) and sets the rest of them to 0 (component
#' knock-out).
#'
#' The amout of components to keep (given as a percentage of the input lenght)
#' can be set by the \code{pcKeepComp}. This will select the first components
#' of the signal, knock-outing the rest. If this value is close to 1, more
#' components will be selected and then more noise will be allowed in the
#' output. For an effective filtering which removes the noise keeping almost
#' all relevant peaks, a value between 0.01 and 0.05 is usually sufficient.
#' Lower values can cause merging of adjacent minor peaks.
#'
#' This library also allows the automatic detection of a fitted value for
#' \code{pcKeepComp}. By default, if uses the \code{pcKeepCompDetect} function,
#' which looks which is the minimum percentage of components than can reproduce
#' the original signal with a corelation between the filtered and the original
#' one of 0.99. See the help page of \code{pcKeepCompDetect} for further
#' details and reference of available parameters.
#'
#' One of the most powerful features of \code{nucleR} is the efficient
#' implementation of the FFT to genomic data. This is achived trought few
#' tweaks that allow an optimum performance of the Fourier Transform. This
#' includes a by-range filtering, an automatic detection of uncovered regions,
#' windowed execution of the filter and padding of the data till nearest power
#' of 2 (this ensures an optimum case for FFT due the high factorization of
#' components). Internal testing showed up that in specific datasets, these
#' optimizations lead to a dramatic improvement of many orders of magnitude
#' (from 3 days to few seconds) while keeping the correlation between the
#' native \code{fft} call and our \code{filterFFT} higher than 0.99. So, the
#' use of these optimizations is highly recomended.
#'
#' If for some reason you want to apply the function without any kind of
#' optimizations you can specify the parameter \code{useOptim=FALSE} to bypass
#' them and get the pure knockout inverse from native FFT call. All other
#' parameters can be still applyied in this case.
#'
#' @param data Coverage or intensities values representing the results of the
#' NGS of TA experiment. This attribute could be a individual vector
#' representing a chromosome (\code{Rle} or \code{numeric} object) or a list of
#' them.
#' @param pcKeepComp Number of components to select, in percentage respect
#' total length of the sample. Allowed values are numeric (in range 0:1) for
#' manual setting or "auto" for automatic detection. See details.
#' @param showPowerSpec Plot the Power Spectrum of the Fast Fourier Transform
#' to visually identify the selected components (see details).
#' @param useOptim This function implements tweaks to a standard fft call to
#' improve (dramatically) the performance in large genomic data. These
#' optimizations can be bypassed by setting this parameter to \code{FALSE}.
#' @param mc.cores If multiple cores are available, maximum number of them to
#' use for parallel processing of \code{data} elements (only useful if
#' \code{data} is a list of elements)
#' @param \dots Other parameters to be passed to \code{pcKeepCompDetect}
#' function
#'
#' @return Numeric vector with cleaned/smoothed values
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}, David Rosell
#' \email{david.rossell@@irbbarcelona.org}
#' @references Smith, Steven W. (1999), The Scientist and Engineer's Guide to
#' Digital Signal Processing (Second ed.), San Diego, Calif.: California
#' Technical Publishing, ISBN 0-9660176-3-3 (availabe online:
#' http://www.dspguide.com/pdfbook.htm)
#' @keywords manip
#' @rdname filterFFT
#'
#' @examples
#' # Load example data, raw hybridization values for Tiling Array
#' raw_data <- get(data(nucleosome_tiling))
#'
#' #Filter data
#' fft_data <- filterFFT(raw_data, pcKeepComp=0.01)
#'
#' # See both profiles
#' par(mfrow=c(2,1), mar=c(3, 4, 1, 1))
#' plot(raw_data, type="l", xlab="position", ylab="Raw intensities")
#' plot(fft_data, type="l", xlab="position", ylab="Filtered intensities")
#'
#' #The power spectrum shows a visual representation of the components
#' fft_data <- filterFFT(raw_data, pcKeepComp=0.01, showPowerSpec=TRUE)
#'
#' @export
#'
setGeneric(
    "filterFFT",
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)
        standardGeneric("filterFFT")
)

#' @rdname filterFFT
setMethod(
    "filterFFT",
    signature(data="SimpleRleList"),
    function (data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE,
            mc.cores=1, ...) {
        return (filterFFT(
            lapply(data, as.vector),
            pcKeepComp,
            showPowerSpec,
            useOptim,
            mc.cores,
            ...
        ))
    }
)

#' @rdname filterFFT
setMethod(
    "filterFFT",
    signature(data="Rle"),
    function (data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)
        filterFFT(as.vector(data), pcKeepComp, showPowerSpec, useOptim, ...)
)

#' @rdname filterFFT
setMethod(
    "filterFFT",
    signature(data="list"),
    function (data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE,
            mc.cores=1, ...) {

        if (length(data) > 1 & showPowerSpec) {
            stop("showPowerSpec only can be applyied to lists of length = 1")
        }

        .xlapply(
            data,
            filterFFT,
            pcKeepComp,
            showPowerSpec,
            useOptim,
            mc.cores=mc.cores,
            ...
        )
    }
)

#' @rdname filterFFT
#' @importFrom graphics plot abline
#' @importFrom IRanges IRanges
#' @importMethodsFrom S4Vectors Rle
#' @importMethodsFrom IRanges start end width end<-
setMethod(
    "filterFFT",
    signature(data="numeric"),
    function (data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE,
            ...) {
        # Check the pcKeepComp
        if (is.numeric(pcKeepComp)) {
            if (pcKeepComp > 1 | pcKeepComp < 0) {
                stop("numeric pcKeepComp must be in range 0:1")
            }
        } else {
            pcKeepComp <- pcKeepCompDetect(data, ...)
        }

        if (showPowerSpec) {
            data[is.na(data)] <- 0
            upper <- min(250000, length(data))
            if (length(data) > 250000) {
                warning(paste0(
                    "warning: only first 250,000bp are represented in the ",
                    "power spectrum"
                ))
            }
            temp <- fft(data[1:upper])
            keep <- round(upper * pcKeepComp)
            plot(
                Re(temp[2:round(length(temp) / 2)]), type="l",
                xlab="Components", ylab="Power",
                sub="Selected components threshold marked as red line"
            )
            abline(v=keep, col="red", lwd=1, lty=2)
        }

        # Bypass all optimizations, very much slower
        if (!useOptim) {
            return (.filterFFT(data, pcKeepComp))
        }

        # Partition the data for available values
        ranges <- IRanges(!is.na(data))
        ranges <- ranges[width(ranges) > 50]  # Discard short regions
        defVal <- NA

        # Split, instead of NAs, the large 0 rows
        if (width(ranges[1]) == length(data)) {
            r <- Rle(data == 0)
            r@values[r@values == TRUE & r@lengths > 500] <- "SPLIT"
            ranges <- IRanges(as.vector(r) != "SPLIT")
            ranges <- ranges[width(ranges) > 100]
            defVal <- 0
        }

        .pad2power2 <- function(x) {
            pad <- rep(0, 2 ^ ceiling(log2(length(x))) - length(x))
            return (c(x, pad))
        }

        # Define FFT by regions for avoid large amount of memory and drop
        # in performance
        .fftRegion <- function(data2, pcKeepComp) {
            # Make windows overlapped 5000bp
            ran <- IRanges(
                start=seq(from=1, to=length(data2), by=(2 ^ 19) - 5000),
                width=2^19
            )
            # Last range has the lenght till the end of data
            end(ran[length(ran)]) <- length(data2)

            # If this gives only one range, return the fft directly
            if (length(ran) == 1) {
                return (.filterFFT(
                    .pad2power2(data2),
                    pcKeepComp
                )[1:end(ran[1])])
            }

            # If last range is shorter than 5000, extend the last-1
            # (we know now we have >1 ranges)
            if (width(ran)[length(ran)] < 5000) {
                end(ran[length(ran) - 1]) <- end(ran[length(ran)])  # Extend
                ran <- ran[1:(length(ran) - 1)]  # Remove last range
            }

            #Iterative case
            res <- NULL
            while (length(ran) > 1) {  # While not last range
                tmp <- .filterFFT(data2[start(ran[1]):end(ran[1])], pcKeepComp)
                if(is.null(res)){
                    res <- tmp
                } else {
                    res <- c(
                        res[1:(length(res) - 2500)],
                        tmp[2501:length(tmp)]
                    )
                }
                ran <- ran[2:length(ran)]
            }

            # Last range
            tmp <- .filterFFT(
                .pad2power2(data2[start(ran[1]):end(ran[1])]),
                pcKeepComp
            )[1:width(ran[1])]

            # If res not set, set it; otherwise append (by construction
            # is larger than 5000bp)
            if (is.null(res)) {
                res <- tmp
            } else {
                res <- c(res[1:(length(res) - 2500)], tmp[2501:length(tmp)])
            }

            return (res)
        }

        fft_ranges <- lapply(
            ranges,
            function(x) .fftRegion(data[x], pcKeepComp)
        )
        # Create a vector of default values and fill it
        res <- rep(defVal, length(data))
        for (i in seq_along(ranges)) {
            res[ranges[[i]]] <- fft_ranges[[i]]
        }

        # Set to default values the positions that have them in the input
        # (remove strange periodicities from large uncovered regions)
        if (is.na(defVal)) {
            rtmp <- IRanges(is.na(data))
            rtmp <- rtmp[width(rtmp) > 15]
            if (length(rtmp) > 0) {
                for (i in 1:length(rtmp)) {
                    res[rtmp[[i]]] <- NA
                }
            }
        } else if (defVal == 0) {
            rtmp <- IRanges(data == 0)
            rtmp <- rtmp[width(rtmp) > 15]
            if (length(rtmp) > 0) {
                for (i in 1:length(rtmp)) {
                    res[rtmp[[i]]] <- 0
                }
            }
            res[res < 0] <- 0
        }

        return (res)
    }
)

#' @importFrom stats fft
.filterFFT <- function(data, pcKeepComp)
{
    data[is.na(data)] <- 0
    temp <- fft(data)
    keep <- round(length(temp) * pcKeepComp)
    temp[keep:(length(temp) - keep)] <- 0
    return (Re(fft(temp, inverse=TRUE)) / length(temp))
}
