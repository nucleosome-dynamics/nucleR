#' Import reads from a list of BAM files.
#'
#' This function allows to load reads from BAM files from both single and
#' paired-end commming from Next Generation Sequencing nucleosome mapping
#' experiments.
#'
#' @param file List of input BAM files.
#' @param type Describes the type of reads. Values allowed are \code{single}
#' for single-ended reads and \code{paired} for pair-ended.
#' @return List of \code{GRanges} containing the reads of each input BAM file.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}, Ricard Illa
#' \email{ricard.illa@@irbbarcelona.org}
#' @keywords file
#'
#' @examples
#' infile <- system.file(
#'     "extdata", "cellCycleM_chrII_5000-25000.bam", package="nucleR"
#' )
#' reads <- readBAM(infile, type="paired")
#'
#' @export readBAM
#'
readBAM <- function (file, type="paired")
{
    if (type == "single") {
        .loadSingleBam(file)
    } else if (type == "paired") {
        .loadPairedBam(file)
    } else {
        stop("type must be `single` or `paired`")
    }
}

irLs2rd <- function(x)
    # Convert a list of IRanges to a RangedData
    IRanges::RangedData(
        ranges=do.call(c, unname(x)),
        space=rep(names(x), sapply(x, length))
    )

sortReads <- function (reads)
{
    # Sort reads RangedData format. Sort them first by chromosome, then by
    # start and then by end
    sortChrs <- function (rans)
        rans[order(names(rans))]
    sortRans <- function (x) {
        tmp <- x[sort.list(IRanges::end(x))]
        tmp[sort.list(IRanges::start(tmp))]
    }
    irLs2rd(lapply(sortChrs(IRanges::ranges(reads)), sortRans))
}

vectorizedAll <- function(...)
    # Helper function that behaves as a vectorized version of the function
    # `all`
    Reduce(`&`, list(...))

.loadSingleBam <- function (exp)
{
    what <- c("pos", "qwidth", "strand", "rname")
    bam <- Rsamtools::scanBam(exp, param=Rsamtools::ScanBamParam(what=what))[[1]]

    non.na <- Reduce(`&`, lapply(bam, Negate(is.na)))
    filtered.bam <- lapply(bam, `[`, non.na)

    # IRanges
    IRanges::RangedData(
        space  = filtered.bam$rname,
        ranges = IRanges::IRanges(
            start = filtered.bam[["pos"]],
            width = filtered.bam[["qwidth"]]
        ),
        strand = filtered.bam[["strand"]]
    )
}

.processStrand <- function (strand, bam)
{
    message(sprintf("processing strand %s", strand))

    p1 <- ifelse(strand == "+", 99, 163)
    p2 <- ifelse(strand == "+", 147, 83)

    unsorted.reads1 <- bam[bam$flag == p1, ]
    unsorted.reads2 <- bam[bam$flag == p2, ]

    rownames(unsorted.reads1) <- as.vector(unsorted.reads1$qname)
    rownames(unsorted.reads2) <- as.vector(unsorted.reads2$qname)

    # Sort by the name of the reads. Assiming the paired reads will have the
    # same name, this will keep the pairs in the same position
    common <- intersect(rownames(unsorted.reads1), rownames(unsorted.reads2))

    reads1 <- unsorted.reads1[common, ]
    reads2 <- unsorted.reads2[common, ]

    # Consistency check
    test <- all(vectorizedAll(
        reads1$mpos  == reads2$pos,
        reads2$mpos  == reads1$pos,
        reads1$rname == reads2$rname
    ))

    if (!test) {
        stop(sprintf(
            "ERROR: Mate selection for %s strand is invalid",
            strand
        ))
    } else {
        IRanges::RangedData(
            space  = as.character(reads1$rname),
            ranges = IRanges::IRanges(
                start = reads1$pos,
                end   = reads2$pos + reads2$qwidth - 1
            )
        )
    }
}

.loadPairedBam <- function (file)
{
    message(sprintf("reading file %s", file))

    what <- c(
        "qname",
        "flag",
        "rname",
        "strand",
        "pos",
        "qwidth",
        "mrnm",
        "mpos"
    )
    bam <- as.data.frame(Rsamtools::scanBam(
        file=file,
        param=Rsamtools::ScanBamParam(what=what)
    )[[1]])

    message("processing flags")
    bam$flag <- bam$flag %% 256

    # Process both strand and return the reads in sorted order
    sortReads(IRanges::rbind(
        .processStrand("+", bam),
        .processStrand("-", bam)
    ))
}
