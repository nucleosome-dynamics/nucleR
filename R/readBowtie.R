#' Import reads from a vector of Bowtie files
#'
#' This function allows to load reads from Bowtie files from both single and
#' paired-end commming from Next Generation Sequencing nucleosome mapping
#' experiments.
#'
#' @param files List of input Bowtie files.
#' @param type Describes the type of reads. Values allowed are `single` for
#'   single-ended reads and `paired` for pair-ended.
#'
#' @return [GenomicRanges::GRangesList] containing the reads of each input BAM
#'   file.
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
#' @keywords file
#'
#' @importFrom GenomicRanges GRangesList
#'
#' @export readBowtie
#'
readBowtie <- function (files, type="paired")
    GRangesList(lapply(strsplit(unlist(files), "/"), .readBowtieFile, type))

#' @importFrom IRanges IRanges
#' @importMethodsFrom ShortRead position chromosome
#' @importMethodsFrom BiocGenerics width
.readBowtieChr <- function (chr, xplus, xminus)
{
    i <- chromosome(xplus) == chr
    sxplus <- xplus[i]
    sxminus <- xminus[i]
    IRanges(start = position(sxplus),
            end   = position(sxminus) + width(sxminus))
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges IRangesList
#' @importMethodsFrom BiocGenerics width strand
#' @importMethodsFrom ShortRead position chromosome
#' @importMethodsFrom ShortRead readAligned
.readBowtieFile <- function (splt, type)
{

    message("Reading ", splt[length(splt)])
    dirPath <- paste(splt[-length(splt)], collapse="/")
    pattern <- sprintf("^%s$", splt[length(splt)])
    bt <- readAligned(dirPath=dirPath, pattern=pattern, type="Bowtie")

    if (type == "paired") {
        i <- seq(1, length(bt), by=2)
        j <- seq(2, length(bt), by=2)
        chrnames <- levels(chromosome(bt))
        res <- lapply(chrnames, .readBowtieChr, bt[i], bt[j])
        names(res) <- chrnames
        return(GRanges(IRangesList(res)))

    } else if (type == "single") {
        return(GRanges(seqnames = chromosome(bt),
                       ranges   = IRanges(start=position(bt),
                                          width=width(bt)),
                       strand   = strand(bt)))
    }
}
