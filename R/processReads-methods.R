setMethod(
    "processReads",
    signature(data="AlignedRead"),
    function(data, type="single", fragmentLen, trim, ...) {

        # require("ShortRead")
        if (missing(fragmentLen)) {
            if (type == "single") {
                message(" * fragmentLen not provided for strand correction, ",
                        ", infering automatically...")
                fragmentLen <- fragmentLenDetect(data, ...)
                message(paste(" * fragmentLen =", fragmentLen))
            } else {
                fragmentLen <- Inf  # Don't remove anything
            }
        }

        #######################################################################
        if (type == "single") {  # Special case for trim==fragmentLength
            if (!missing(trim)) {
                if (trim == fragmentLen) {
                    start <- position(data)
                    start[strand(data) == "-"] <- position(data) +
                                                  nchar(sread(data))
                    res <- RangedData(ranges=IRanges(start=start,
                                                     width=fragmentLen),
                                      space=as.character(chromosome(data)))
                    return(res)
                }
            }

            # If no trim restriction, use whole read length,
            # else use trim length
            sr_len <- nchar(sread(data))

            if (!missing(trim)) {
                sr_len <- trim
            }

            # Floor, then we will add 1 if needed
            shift  <- floor((fragmentLen - sr_len) / 2)

            strand <- strand(data)

            new_levels <- vector(mode="integer", length=3)
            new_levels[which(levels(strand) == "-")] <- -1
            new_levels[which(levels(strand) == "+")] <- +1
            new_levels[which(levels(strand) == "*")] <- 0
            levels(strand) <- new_levels

            strand <- as.numeric(levels(strand))[strand]

            # Correction for non rounded shifts in negative strand (ie 147/40)
            # If the fragment_length - read_length is odd and negative strand,
            # add 1 extra base
            correct <- rep(0, length(shift))
            correct[(strand == -1) & ((fragmentLen - sr_len) %% 2 == 1)] <- 1

            start <- position(data) + ((shift + correct) * strand)
            res <- RangedData(ranges=IRanges(start=start, width=sr_len),
                              space=as.character(chromosome(data)))

        #######################################################################
        } else if (type == "paired") {
            sr_len <- nchar(sread(data))
            selection <- sr_len < fragmentLen
            space <- as.character(chromosome(data)[selection])

            # Normal start or trimmed version if requested
            start <- position(data)[selection]
            if (!missing(trim)) {
                start <- start + (sr_len[selection] / 2 - floor(trim / 2))
            }

            width <- sr_len[selection]
            if (!missing(trim)) {
                width <- trim
            }
            res <- RangedData(ranges=IRanges(start=start,
                                             width=width),
                              space=space)

        #######################################################################
        } else {
            stop(paste("type must be 'single' for single-ended data or",
                       "'paired' for paired-ended"))
        }

        return(res)
    }
)

setMethod(
    "processReads",
    signature(data="RangedData"),
    function(data, type="single", fragmentLen, trim, ...) {

        if (missing(fragmentLen)) {
            if (type == "single") {
                message(" * fragmentLen not provided for strand correction, ",
                        "infering automatically...")
                fragmentLen <- fragmentLenDetect(data,...)
                message(paste(" * fragmentLen =", fragmentLen))
            } else {
                fragmentLen <- Inf  # Don't remove anything
            }
        }

        #######################################################################
        if (type == "single") {  # Special case for trim==fragmentLength
            if (!missing(trim)) {
                if (trim == fragmentLen) {
                    start <- start(data)
                    negStrand <- data$strand == "-"
                    start[negStrand] <- end(data)[negStrand] - fragmentLen
                    res <- RangedData(ranges=IRanges(start=start,
                                                     width=fragmentLen),
                                      space=space(data))
                    return(res)
                }
            }

            # If no trim restriction, use whole read length,
            # else use trim length
            sr_len <- width(data)
            if (!missing(trim)) {
                sr_len <- trim
            }

            # Floor, then we will add 1 if needed
            shift <- floor((fragmentLen - sr_len) / 2)

            strand <- rep(1, nrow(data))
            strand[data$strand == "-"] <- -1

            # Correction for non rounded shifts in negative strand (ie 147/40)
            # If the fragment_length - read_length is odd and negative strand,
            # add 1 extra base
            correct <- rep(0, length(strand))
            correct[(strand == -1) & ((fragmentLen - sr_len) %% 2 == 1)] <- 1

            start <- start(data) + ((shift + correct) * strand)
            res <- RangedData(ranges=IRanges(start=start,
                                             width=sr_len),
                              space=space(data))

        #######################################################################
        } else if (type=="paired") {
            data <- data[width(data) < fragmentLen, ]
            sr_len <- width(data)

            # Normal start or trimmed version if requested
            start <- start(data)
            if (!missing(trim)) {
                start <- start + ((sr_len / 2) - floor(trim / 2))
            }

            width <- sr_len

            if (!missing(trim)) {
                width <- trim
            }

            res <- RangedData(ranges=IRanges(start=start,
                                             width=width),
                              space=space(data))
        #######################################################################
        } else {
            stop(paste("type must be 'single' for single-ended data or",
                       "'paired' for paired-ended"))
        }

        return(res)
    }
)
