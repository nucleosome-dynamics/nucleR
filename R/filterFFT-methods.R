setMethod(
    "filterFFT",
    signature(data="SimpleRleList"),
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE,
             mc.cores=1, ...) {
        temp <- lapply(data, as.vector)
        return(filterFFT(temp, pcKeepComp, showPowerSpec, useOptim, mc.cores,
                         ...))
    }
)

setMethod(
    "filterFFT",
    signature(data="Rle"),
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...)
        filterFFT(as.vector(data), pcKeepComp, showPowerSpec, useOptim, ...)
)

setMethod(
    "filterFFT",
    signature(data="list"),
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE,
             mc.cores=1, ...) {

        if (length(data) > 1 & showPowerSpec) {
            stop("showPowerSpec only can be applyied to lists of length = 1")
        }

        .xlapply(
            data,
            filterFFT,
            pcKeepComp, showPowerSpec, useOptim,
            mc.cores=mc.cores,
            ...
        )
    }
)

# This is the filter itself
.filterFFT <- function(data, pcKeepComp)
{
    data[is.na(data)] <- 0
    temp <- fft(data)
    keep <- round(length(temp) * pcKeepComp)
    temp[keep:(length(temp) - keep)] <- 0
    return(Re(fft(temp, inverse=TRUE)) / length(temp))
}

setMethod(
    "filterFFT",
    signature(data="numeric"),
    function(data, pcKeepComp="auto", showPowerSpec=FALSE, useOptim=TRUE, ...) {
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
                warning("warning: only first 250,000bp are represented in the power spectrum")
            }
            temp <- fft(data[1:upper])
            keep <- round(upper * pcKeepComp)
            plot(Re(temp[2:round(length(temp) / 2)]), type="l",
                 xlab="Components", ylab="Power",
                 sub="Selected components threshold marked as red line")
            abline(v=keep, col="red", lwd=1, lty=2)
        }

        # Bypass all optimizations, very much slower
        if (!useOptim) {
            return(.filterFFT(data, pcKeepComp))
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
            return(c(x, pad))
        }

        # Define FFT by regions for avoid large amount of memory and drop
        # in performance
        .fftRegion <- function(data2, pcKeepComp) {
            # Make windows overlapped 5000bp
            ran <- IRanges(start=seq(from=1,
                                     to=length(data2),
                                     by=(2 ^ 19) - 5000),
                           width=2 ^ 19)
            # Last range has the lenght till the end of data
            end(ran[length(ran)]) <- length(data2)

            # If this gives only one range, return the fft directly
            if (length(ran) == 1) {
                return(.filterFFT(.pad2power2(data2),
                                  pcKeepComp)[1:end(ran[1])])
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
                    res <- c(res[1:(length(res) - 2500)], tmp[2501:length(tmp)])
                }
                ran <- ran[2:length(ran)]
            }

            # Last range
            tmp <- .filterFFT(.pad2power2(data2[start(ran[1]):end(ran[1])]),
                              pcKeepComp)[1:width(ran[1])]

            # If res not set, set it; otherwise append (by construction
            # is larger than 5000bp)
            if (is.null(res)) {
                res <- tmp
            } else {
                res <- c(res[1:(length(res) - 2500)], tmp[2501:length(tmp)])
            }

            return(res)
        }

        fft_ranges <- lapply(
            ranges,
            function(x)
                .fftRegion(data[x], pcKeepComp)
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
                for(i in 1:length(rtmp)) {
                    res[rtmp[[i]]] <- NA
                }
            }
        } else if (defVal == 0) {
            rtmp <- IRanges(data == 0)
            rtmp <- rtmp[width(rtmp) > 15]
            if (length(rtmp) > 0) {
                for(i in 1:length(rtmp)) {
                    res[rtmp[[i]]] <- 0
                }
            }
            res[res < 0] <- 0
        }

        return(res)
    }
)
