#' Import reads from a list of BAM files.
#'
#' This function allows to load reads from BAM files from both single and
#' paired-end commming from Next Generation Sequencing nucleosome mapping
#' experiments.
#'
#' @param file List of input BAM files.
#' @param type Describes the type of reads. Values allowed are `single` for
#'   single-ended reads and `paired` for pair-ended.
#'
#' @return List of `GRanges` containing the reads of each input BAM file.
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
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


#' Vectorized version of `all`
#'
#' Helper function that behaves as a vectorized version of the function `all`
#'
#' @param ... arbitraty amount of `logical` vectors, expected to have the same
#'   length
#' @return `logical` vector
#'
.vectorizedAll <- function(...)
    Reduce(`&`, list(...))


#' Load a single-end BAM
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom Rsamtools scanBam ScanBamParam
#'
.loadSingleBam <- function (exp)
{
    what <- c("pos", "qwidth", "strand", "rname")
    bam <- scanBam(exp, param=ScanBamParam(what=what))[[1]]

    non.na <- Reduce(`&`, lapply(bam, Negate(is.na)))
    filtered.bam <- lapply(bam, `[`, non.na)

    GRanges(
        seqnames = filtered.bam$rname,
        ranges = IRanges(
            start = filtered.bam[["pos"]],
            width = filtered.bam[["qwidth"]]
        ),
        strand = filtered.bam[["strand"]]
    )
}

#' Process a given strand from a BAM file to read
#'
#' @param strand either strand "+" or "-"
#' @param bam `data.frame` representing the BAM file as read by
#'    [Rsamtools::scanBam]
#' @return [GenomicRanges::GRanges] object representing the reads in a given
#'   strand
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
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
    test <- all(.vectorizedAll(
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
        GRanges(
            seqnames = reads1$rname,
            ranges   = IRanges(
                start = reads1$pos,
                end   = reads2$pos + reads2$qwidth - 1
            )
        )
    }
}

#' Load a paired-end-end BAM
#'
#' @importMethodsFrom Rsamtools scanBam ScanBamParam
#' @importMethodsFrom GenomeInfoDb sortSeqlevels
#' @importMethodsFrom BiocGenerics sort
#'
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
    bam <- as.data.frame(scanBam(
        file=file,
        param=ScanBamParam(what=what)
    )[[1]])

    message("processing flags")
    bam$flag <- bam$flag %% 256

    # Process both strand and return the reads in sorted order
    both.strands <- c(.processStrand("+", bam),
                      .processStrand("-", bam))
    sort(sortSeqlevels(both.strands))
}
