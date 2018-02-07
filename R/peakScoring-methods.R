setMethod(
    "peakScoring",
    signature(peaks="list"),
    function (peaks, data, threshold="25%", mc.cores=1)
        # Return the list directly
        .xlapply(
            peaks,
            peakScoring,
            data      = data,
            threshold = threshold,
            mc.cores  = mc.cores
        )
)

setMethod(
    "peakScoring",
    signature(peaks="IRangesList"),
    function (peaks, data, threshold="25%", weight.width=1, weight.height=1,
            dyad.length=38, mc.cores=1) {

        res <- .xlapply(
            names(peaks),
            function(x)
                peakScoring(
                    peaks         = peaks[[x]],
                    data          = data[[x]],
                    threshold     = threshold,
                    dyad.length   = dyad.length,
                    weight.width  = weight.width,
                    weight.height = weight.height
                ),
            mc.cores=mc.cores
        )
        names(res) <- names(peaks)

        # Result should be returned as a single RangedData object
        # Put the correct space in RangedData
        for (name in names(res)) {
            if (nrow(res[[name]])) {
                names(res[[name]]) <- name
            }
        }

        # Combine RangedData objects into single one
        do.call(c, unname(res))
    }
)

setMethod(
    "peakScoring",
    signature(peaks="numeric"),
    function (peaks, data, threshold="25%") {

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshold is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        mean <- mean(data[!is.na(data) & data > threshold], na.rm=TRUE)
        sd <- sd(data[!is.na(data) & data > threshold], na.rm=TRUE)

        res <- pnorm(data[peaks], mean=mean, sd=sd, lower.tail=TRUE)
        return (data.frame(peak=peaks, score=res))
    }
)

setMethod(
    "peakScoring",
    signature(peaks="IRanges"),
    function (peaks, data, threshold="25%", weight.width=1, weight.height=1,
            dyad.length=38) {

        # Calculate the ranges in threshold and get the coverage
        if (!is.numeric(threshold)) {
            # If threshdol is given as a string with percentage, convert it
            if (grep("%", threshold) == 1) {
                threshold <- quantile(
                    data,
                    as.numeric(sub("%", "", threshold)) / 100,
                    na.rm=TRUE
                )
            }
        }

        # Hack for TA, that could have negative values in range calculations
        if (min(data, na.rm=TRUE) < 0) {
            data <- data + abs(min(data, na.rm=TRUE))
        }

        mean <- mean(data[data > threshold], na.rm=TRUE)
        sd <- sd(data[data > threshold], na.rm=TRUE)

        # Calculate dyad range
        dyad.middl <- start(peaks) + floor(width(peaks) / 2)
        dyad.start <- dyad.middl - floor((dyad.length / 2))
        dyad.end <- dyad.start + (dyad.length - 1)
        dyads <- IRanges(start=dyad.start, end=dyad.end)

        sums.range <- lapply(peaks, function(x) mean(data[x], na.rm=TRUE))
        sums.dyad  <- lapply(dyads, function(x) mean(data[x], na.rm=TRUE))

        # Score the heigh of the peak
        scor.heigh <- pnorm(data[dyad.middl], mean=mean, sd=sd, lower.tail=TRUE)

        # Score the width (dispersion) of the peak
        scor.width <- unlist(sums.dyad) / unlist(sums.range)
        scor.width <- scor.width / max(scor.width)

        #Final score
        sum.wei <- weight.width + weight.height
        scor.final <- ((scor.heigh * weight.height) / sum.wei) +
            ((scor.width * weight.width) / sum.wei)

        # Return everything or just merged score
        return (RangedData(
            peaks,
            score   = scor.final,
            score_w = scor.width,
            score_h = scor.heigh
        ))
    }
)
