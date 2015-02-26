setMethod(
    "peakDetection",
    signature(data="list"),
    function(data, threshold="25%", width=1, score=TRUE, mc.cores=1) {

        res <- .xlapply(
            data,
            peakDetection,
            threshold=threshold, width=width, score=score,
            mc.cores=mc.cores
        )

        # Process the result, case with ranges
        if (width > 1) {
            # res is a list of IRanges, conversion is direct
            if (score == FALSE) {
                return(IRangesList(unlist(res)))

            # res is a list of RangedData
            } else {
                # Put the correct space in RangedData by changing range names
                for (name in names(res)) {
                    names(res[[name]]@ranges) <- name
                }

                # Combine RangedData objects into single one
                if (length(res) > 1) {
                    temp <- res[[1]]
                    for (i in 2:length(res)) {
                        temp <- c(temp, res[[i]])
                    }
                } else {
                    temp <- res[[1]]
                }

                return(temp)
            }

        #Case without ranges, anyway return a list
        } else {
            return(res)
        }
    }
)

setMethod(
    "peakDetection",
    signature(data="numeric"),
    function (data, threshold="25%", width=1, score=TRUE, mc.cores=1) {

        if (width < 1) {
            stop("'width' attribute should be greater than 1")
        }

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshdol is given as a string with percentage, convert it
            if(grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        ranges <- IRanges(!is.na(data) & data > threshold)
        if (length(ranges) == 0) {
            return(NULL)
        }

        covers <- lapply(ranges, function(x) data[x])

        # For each range, look for changes of trend and keep the starting
        # position of trend change
        pea <- .xlapply(
            covers,
            function(x)
                if (length(x) == 1) {
                    1
                } else {
                    start(IRanges(x[2:length(x)] <
                                  x[1:(lenth(x) - 1)]))
                },
            mc.cores=mc.cores
        )

        # Some peaks can have only one trend, correct them
        unitrend <- which(sapply(pea, function(x) length(x) == 0))
        pea[unitrend] <- sapply(covers[unitrend],
                                function(x) which(x == max(x)))[1]

        # Add start offset to peaks relative to the start of the range
        starts <- start(ranges)
        res <- unlist(sapply(1:length(starts),
                             function(i) pea[[i]] + starts[[i]]))

        # Extension
        if (width > 1) {
            ext <- floor(width / 2)
            starts <- res - ext
            # Odd/pair correction
            ends <- res + ifelse(width %% 2 == 0, ext - 1, ext)
            res <- IRanges(start=starts, end=ends)
            # Remove out of bounds
            res <- res[start(res) > 1 & end(res) < length(data)]
        }

        if (score) {
            return(peakScoring(peaks=res, data=data, threshold=threshold))
        } else {
            return(res)
        }
    }
)
