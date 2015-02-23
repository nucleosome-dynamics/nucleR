library(multicore)
library(ShortRead)
library(Rsamtools)

int2base <- function(x, b = 2)
{   # Binary conversion
    xi <- as.integer(x)
    if (any(is.na(xi) | ((x - xi) != 0))) {
        print(list(ERROR="x not integer", x=x))
    }
    N <- length(x)
    xMax <- max(x)
    ndigits <- floor(logb(xMax, base=2)) + 1
    Base.b <- array(NA, dim=c(N, ndigits))
    for(i in 1:ndigits) {
        Base.b[, ndigits - i + 1] <- (x %% b)
        x <- (x %/% b)
    }
    if (N == 1) {
        Base.b[1, ]
    } else {
        Base.b
    }
}

bamFlagMatrix <- function(flags, file)
{   # SAM/BAM flag matrix
    message("** ", file, " -> Processing flags")
    bin <- int2base(flags)
    n <- ncol(bin)
    colnames(bin) <- c(rev(names(formals(scanBamFlag))[1:n]))
    bin
}

getPos <- function(flags, s, p)
{   # Given a matrix of flags, a strand ('+' or '-') and a position
    # ('1' or '2'), # return a vector of bools that fullfills the requirements

    # a simple closure to keep things less verbose and repetitive
    getFlag <- function(x) flags[, x]

    is.pair <- getFlag("isPaired") & getFlag("isProperPair")

    if (p == 1) {
        is.position <- getFlag("isFirstMateRead")
    } else if (p == 2) {
        is.position <- getFlag("isSecondMateRead")
    }

    minStrand <- getFlag("isMinusStrand")

    if ((p == 1 && s == "+") ||
        (p == 2 && s == "-")) {
        right.strand <- !minStrand
    } else if ((p == 1 && s == "-") ||
               (p == 2 && s == "+")) {
        right.strand <- minStrand
    }

    is.pair & is.position & right.strand
}

checkConsistency <- function(reads1, reads2)
    # Returns true if the reads positions are consistent
    all(reads1$mpos == reads2$pos &
        reads2$mpos == reads1$pos &
        reads1$rname == reads2$rname)

buildReads <- function(reads1, reads2, strand)
{   # Given a strand, build the reads in GRanges format accordingly
    if (strand == "+") {
        ranges <- IRanges(start=reads1[["pos"]],
                          end=reads2[["pos"]] + reads2[["qwidth"]] - 1)
    } else if (strand == "-") {
        ranges <- IRanges(start=reads2[["pos"]],
                          end=reads1[["pos"]] + reads1[["qwidth"]] - 1)
    }
    names <- as.character(reads1$rname)
    GRanges(seqnames=names, ranges=ranges)
}

processStrand <- function(flags, bam, strand, file)
{   # Given a matrix of flags and bam information, process given flag
    message("** ", file, " -> Processing ", strand, " strand")

    pos1 <- getPos(flags, strand, 1)
    pos2 <- getPos(flags, strand, 2)

    reads1 <- lapply(bam, "[", pos1)
    reads2 <- lapply(bam, "[", pos2)

    if (!checkConsistency(reads1, reads2)) {
        stop(paste("ERROR: Mate selection for", strand, "strand is invalid"))
    }

    buildReads(reads1, reads2, strand)
}

readBamFile <- function(file)
{   # Read a BAM file (only one access to disk, intended for Shared Memory)
    message("** ", file, "-> Reading")
    what <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mrnm",
              "mpos")
    scanBam(file=file, param=ScanBamParam(what=what))[[1]]
}

removeRepeated <- function(bam, flags)
{   # remove repeated reads (marked as isPrimaryRead == 1), there might be none
    x <- "isPrimaryRead"
    if (x %in% colnames(flags)) {
        mmids <- bam[["qname"]] %in% unique(bam$qname[flags[, x] == 1])
        mmpos <- which(bam[["qname"]] %in% mmids)
        list(bam=lapply(bam, "[", -mmpos),
             flags=flags[-mmpos, ])
    } else {
        list(bam=bam, flags=flags)
    }
}

processPairedEnd <- function(file)
{   # Process a paired end Bam file
    bam <- readBamFile(file)
    flags <- bamFlagMatrix(bam$flag, file)
    non.repeated <- removeRepeated(bam, flags)
    with(non.repeated,
         sort(c(processStrand(flags, bam, "+", file),
                processStrand(flags, bam, "-", file))))
}

in.dir <- "/home/rilla/scratch/nucleosome_dynamics/bam_out"
mc.cores <- 5

exps <- dir(in.dir, pattern="bam$", full.names=TRUE)

newreads <- processPairedEnd(exps[1])
