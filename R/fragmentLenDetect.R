setGeneric(
    "fragmentLenDetect",
    function(reads, samples=1000, window=5000, min.shift=1, max.shift=100,
            mc.cores=1, as.shift=FALSE)
        standardGeneric("fragmentLenDetect")
)

setMethod(
    "fragmentLenDetect",
    signature(reads="AlignedRead"), 
    function (reads, samples = 1000, window = 1000, min.shift = 1,
            max.shift = 100, mc.cores = 1, as.shift = FALSE) {

        # Randomly select regions in the available chromosome bounds
        chrSample <- as.character(sample(chromosome(reads), samples))
        chrsLen <- sapply(
            levels(chromosome(reads)),
            function(x)
                max(position(reads[chromosome(reads) == x]))
        )
        position <- round(runif(samples, max=chrsLen[chrSample] - window))
        dd <- data.frame(chrSample, position)

        # For each sampled region, look for the shift with higher correlation
        # between strands
        shiftPos <- function (i) {
            chr <- as.character(dd[i, "chrSample"])
            sta <- as.numeric(dd[i, "position"])
            end <- sta + window
            rea <- reads[
                chromosome(reads) == chr &
                position(reads) > sta &
                position(reads) < end
            ]
            if (length(rea) == 0) {
                return (numeric(0))  # Discard uncovered regions
            }

            cpos <- try(
                as.vector(
                    coverage(rea[as.character(strand(rea)) == "+"])[[1]]
                )[sta:end],
                silent=TRUE
            )
            cneg <- try(
                as.vector(
                    coverage(rea[as.character(strand(rea)) == "-"])[[1]]
                )[sta:end],
                silent=TRUE
            )

            # Avoid unpaired strands (ie, all strands being "+")
            if (class(cpos) == "try-error" | class(cneg) == "try-error") {
                return(numeric(0))
            }

            cpos[is.na(cpos)] <- 0
            cneg[is.na(cneg)] <- 0

            x <- sapply(
                min.shift:max.shift,
                function(s)
                    cor(cpos[1:(window + 1 - s)], cneg[s:window])
            )
            res <- which(x == max(x))[1]

            # We only shifted the negative strand, but in real case both
            # strands will be shifted the half of this amount
            res <- res / 2

            #Discard NAs
            if (is.na(res) | !is.numeric(res)) {
                return(numeric(0))
            }

            return(res + min.shift - 1)
        }

        shift <- round(mean(unlist(.xlapply(1:nrow(dd),
                                            shiftPos,
                                            mc.cores=mc.cores))))

        #Fragment length is the shift * 2 + the length of the read
        fragLen <- shift * 2 + width(reads)[1]

        if (as.shift) {
            return(shift)
        } else {
            return(fragLen)
        }

    }
)

setMethod(
    "fragmentLenDetect",
    signature(reads="RangedData"),
    function (reads, samples=1000, window=1000, min.shift=1, max.shift=100,
            mc.cores=1, as.shift=FALSE) {

        # Calculate the whole coverage here saves cpu and memory later for big
        # genomes. This improves a lot the performance on big genomes if the
        # sampling value is big. For small datasets could be less efficient than
        # obvious approach, but it is worth it
        .doCover <- function(chr) {
            df <- as.data.frame(reads[chr])
            dp <- df[df$strand == "+",]
            dn <- df[df$strand == "-",]
            rp <- IRanges(start=dp$start, width=dp$width)
            rn <- IRanges(start=dn$start, width=dn$width)
            cp <- coverage(rp, method="hash")
            cn <- coverage(rn, method="hash")
            return(list(pos=cp, neg=cn))
        }

        cover <- .xlapply(names(reads), .doCover, mc.cores=mc.cores)
        names(cover) <- names(reads)

        # Randomly select regions in the available chromosome bounds
        count <- sapply(reads, nrow)
        probs <- count / sum(count)
        chrSample <- sample(names(reads), samples, replace=TRUE, prob=probs)
        chrLength <- sapply(
            cover,
            function(x)
                min(length(x$pos), length(x$neg))
        )
        position <- round(runif(samples, max=chrLength[chrSample] - window))
        dd <- data.frame(chrSample, position)


        # For each sampled region, look for the shift with higher correlation
        # between strands
        shiftPos <- function(i) {
            chr <- as.character(dd[i, "chrSample"])
            sta <- as.numeric(dd[i, "position"])
            end <- sta + window

            cpos <- try(as.vector(cover[[chr]][["pos"]][sta:end]), silent=TRUE)
            cneg <- try(as.vector(cover[[chr]][["neg"]][sta:end]), silent=TRUE)

            if (class(cpos) == "try-error" | class(cneg) == "try-error") {
                return(NA)
            }
            if (sum(cpos) == 0 | sum(cneg) == 0) {
                return(NA)
            }

            cpos[is.na(cpos)] <- 0
            cneg[is.na(cneg)] <- 0

            x <- sapply(
                min.shift:max.shift,
                function (s)
                    cor(cpos[1:(window + 1 - s)], cneg[s:window])
            )
            res <- which(x == max(x, na.rm=TRUE))[1]

            # We only shifted the negative strand, but both strands will be
            # shifted the half of this amount
            res <- res / 2

            # Discard NAs
            if (is.na(res) | !is.numeric(res)) {
                return(numeric(0))
            }

            return(res + min.shift - 1)
        }

        shift <- round(mean(unlist(.xlapply(1:nrow(dd),
                                            shiftPos,
                                            mc.cores=mc.cores)),
                            na.rm=TRUE))

        #Fragment length is the shift * 2 + the length of the read
        fragLen <- shift * 2 + width(reads)[1]

        if (as.shift) {
            return(shift)
        } else {
            return(fragLen)
        }
    }
)
