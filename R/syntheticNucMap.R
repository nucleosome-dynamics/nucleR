syntheticNucMap <- function (wp.num = 100, wp.del = 10, wp.var = 20,
                            fuz.num = 50, fuz.var = 50, max.cover = 20,
                            nuc.len = 147, lin.len = 20, rnd.seed = NULL,
                            as.ratio = FALSE, show.plot = FALSE, ...)
{
    # Set random seed if given
    if (!is.null(rnd.seed)) {
        set.seed(rnd.seed)
    }

    # WELL POS NUCLEOSOMES
    # Starting point of putative nucleosomes
    wp.starts <- (nuc.len + lin.len) * seq(0, wp.num - 1) + 1

    # How many times a read is repeated
    wp.nreads <- round(runif(wp.num, min=1, max=max.cover))

    # Delete some reads (set repetition times to 0)
    wp.nreads[round(runif(wp.del, min=0, max=wp.num))] <- 0

    # Set each nucleosome as a repeated single start position
    wp.repstar <- rep(wp.starts, wp.nreads)

    # Add some variance to the starting points
    var <- round(runif(length(wp.repstar), min=-wp.var, max=wp.var))
    wp.varstar <- wp.repstar + var

    # Putative reads
    wp.reads <- IRanges(start=wp.varstar, width=nuc.len)

    # OVERLAPPED (FUZZY) NUCLEOSOMES
    # Starting point of fuzzy nucleosomes (random)
    fuz.starts <- round(runif(
        fuz.num,
        min=1,
        max=(nuc.len + lin.len) * wp.num
    ))

    # How many times a read is repeated
    fuz.nreads <- round(runif(fuz.num, min=1, max=max.cover))

    # Set each nucleosome as a repeated single start position
    fuz.repstar <- rep(fuz.starts, fuz.nreads)

    # Add some variance to the starting points
    var <- round(runif(length(fuz.repstar), min=-fuz.var, max=fuz.var))
    fuz.varstar <- fuz.repstar + var

    # Overlapped reads
    fuz.reads <- IRanges(start=fuz.varstar, width=nuc.len)

    # ALL SYNTHETIC READS
    syn.reads <- c(wp.reads, fuz.reads)

    # RATIO AS HYBRIDIZATION (Tiling Array)
    if (as.ratio) {
        # Just put the same amount of reads as before randomly
        ctr.starts <- round(runif(
            length(syn.reads),
            min=1,
            max=max(start(syn.reads))
        ))

        # This time use a random read length, between 50 and 250 
        ctr.widths <- round(runif(length(syn.reads), min=50, max=250))

        # "Control reads"
        ctr.reads <- IRanges(start=ctr.starts, width=ctr.widths)

        # ratio
        syn.ratio <- suppressWarnings(
            log2(as.vector(coverage(syn.reads))) -
            log2(as.vector(coverage(ctr.reads)))
        )

        syn.ratio[abs(syn.ratio) == Inf] <- NA  # Some lost bases... as reality
        syn.ratio <- Rle(syn.ratio)
    }

    result <- list()

    result[["wp.starts"]] <- wp.starts
    result[["wp.nreads"]] <- wp.nreads
    result[["wp.reads"]] <- wp.reads

    result[["fuz.starts"]] <- fuz.starts
    result[["fuz.nreads"]] <- fuz.nreads
    result[["fuz.reads"]] <- fuz.reads

    result[["syn.reads"]] <- syn.reads

    if(as.ratio) {
        result[["ctr.reads"]] <- ctr.reads
        result[["syn.ratio"]] <- syn.ratio
    }

    if(show.plot) {
        # Y-lim range
        max <- max(coverage(syn.reads), na.rm=TRUE)
        min <- 0
        if (as.ratio) {
            min <- min(syn.ratio, na.rm=TRUE)
        }

        # Plot main (coverage)
        plot(
            as.vector(coverage(syn.reads)), type="h", col="#AADDAA",
            ylim=c(min,max),
            ...
        )

        # Plot ratio, if asked for
        if (as.ratio) {
            lines(as.vector(syn.ratio), type="l", col="darkorange", lwd=2)
        }
        if (as.ratio) {
            abline(h=0, col="darkorange4")
        }

        # Plot nucleosome positions (dyad)
        points(wp.starts + 74, wp.nreads, col="red", pch=19)
        points(fuz.starts + 74, fuz.nreads, col="blue", pch=20)

        # Plot legend
        if (as.ratio) {
            legend(
                "top",
                c("Coverage", "Ratio", "Well-pos", "Fuzzy"),
                fill=c("#AADDAA", "darkorange", "red", "blue"), bty="n",
                horiz=TRUE
            )
        } else {
            legend(
                "top",
                c("Coverage", "Well-pos", "Fuzzy"),
                fill=c("#AADDAA", "red", "blue"), bty="n", horiz=TRUE
            )
        }
    }

    return (result)
}
