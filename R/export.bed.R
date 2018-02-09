#' Export ranges in BED format
#'
#' Export ranges in BED format, compatible with UCSC genome browser, IGB, and
#' others.
#'
#' @param ranges Ranges to export, in \code{IRanges}, \code{IRangesList} or
#' \code{RangedData} format
#' @param score Score data if not included in \code{ranges} object. Bed file
#' will put all scores=1000 if scores are not present
#' @param chrom For single \code{IRanges} objects, the chromosome they
#' represent. For other data types, values from \code{names(...)} will be used.
#' @param name Name of the track
#' @param desc Description of the track
#' @param filepath Path and prefix of the file(s) to write. Chromosome number
#' and "bed" extension will be automatically added.
#' @param splitByChrom If multiple chromosomes are given, should they be
#' splitted into one file per chromosome or shall them be saved all together?
#'
#' @return (none)
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @references BED format specification:
#' http://genome.ucsc.edu/FAQ/FAQformat#format1
#' @keywords file
#' @rdname export.bed
#'
#' @examples
#' # Generate random ranges with scores
#' ran <- RangedData(IRanges(start=1:100, end=101:200), score=(1:100) / 100)
#' names(ran) <- "chrX"
#'
#' # Export as bed file
#' export.bed(ran, name="test_track", desc="Just a test track")
#'
#' # If executed, this would create a file named "test_track.chrX.bed" with:
#'
#' #    track name="test_track" description="Just a test track" useScore=0
#' # chrX    1 101 nucl1 0.01
#' # chrX    2 102 nucl2 0.02
#' # chrX    3 103 nucl3 0.03
#' # ...
#'
#' @export
#'
setGeneric(
    "export.bed",
    function(ranges, score=NULL, chrom, name, desc=name, filepath=name,
            splitByChrom=TRUE)
        standardGeneric("export.bed")
)

#' @rdname export.bed
setMethod(
    "export.bed",
    signature(ranges="IRanges"),
    function (ranges, score=NULL, chrom, name, desc=name, filepath=name) {
        dd <- as.data.frame(ranges)
        dd$chrom <- chrom
        dd$score <- 1000
        if (!is.null(score)) {
            dd$score <- floor(score * 1000)
        }
        dd$count <- seq(1, nrow(dd))

        .export.bed(
            df=dd,
            name=name,
            desc=desc,
            filename=paste(filepath, chrom, sep=".")
        )
    }
)

#' @rdname export.bed
setMethod(
    "export.bed",
    signature(ranges="CompressedIRangesList"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
            splitByChrom=TRUE) {

        if (splitByChrom) {
            for(chr in names(ranges)) {
                export.bed(
                    ranges=ranges[[chr]],
                    score=score[[chr]],
                    chrom=chr,
                    name=name,
                    desc=desc,
                    filepath=filepath,
                    splitByChrom=FALSE
                )
            }
        } else {
            dd <- as.data.frame(unlist(ranges))
            dd$chrom <- dd$names
            dd$score <- 1000
            if (!is.null(score)) {
                dd$score <- unlist(score)
            }
            dd$count <- seq(1, nrow(dd))

            .export.bed(df=dd, name=name, desc=desc, filename=filepath)
        }
    }
)

#' @rdname export.bed
setMethod(
    "export.bed",
    signature(ranges="RangedData"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
            splitByChrom=TRUE) {

        if (splitByChrom) {
            for(chr in names(ranges)) {
                export.bed(
                    ranges=ranges[chr],
                    chrom=chr,
                    name=name,
                    desc=desc,
                    filepath=paste(filepath, chr, sep="."),
                    splitByChrom=FALSE
                )
            }
        } else {
            dd <- as.data.frame(ranges)
            dd$chrom <- dd$space
            dd$count <- seq(1, nrow(dd))

            .export.bed(df=dd, name=name, desc=desc, filename=filepath)
        }
    }
)

#' @rdname export.bed
setMethod(
    "export.bed",
    signature(ranges="GRanges"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
            splitByChrom=TRUE) {

        if (splitByChrom) {
            #for(chr in IRanges::levels(seqnames(ranges))) {
            for (chr in levels(GenomeInfoDb::seqnames(ranges))) {
                export.bed(
                    ranges[GenomeInfoDb::seqnames(ranges) == chr],
                    chrom=chr,
                    name=name,
                    desc=desc,
                    filepath=paste(filepath, chr, sep="."),
                    splitByChrom=FALSE
                )
            }
        } else {
            dd <- as.data.frame(ranges)
            dd$chrom <- dd$seqnames
            dd$count <- seq(1, nrow(dd))

            .export.bed(df=dd, name=name, desc=desc, filename=filepath)
        }
    }
)

.export.bed <- function(df, name, desc, filename)
{
    hasScore <- "score" %in% names(df)
    tscore <- ifelse(hasScore, "useScore=1", "useScore=0")

    tdesc <- paste('description="', desc, '"', sep="")
    tname <- paste('name="', name, '"', sep="")

    header <- paste("track", tname, tdesc, tscore, sep="\t")


    .printLine <- function(start, end, chr, score, num)
        cat(paste(
            chr,
            start,
            end,
            paste("nucl", as.numeric(num), sep=""),
            score,
            sep="\t"
        ), sep="\n")

    sink(paste(filename, "bed", sep="."))
    cat(header, sep="\n")
    res <- apply(
        df,
        1,
        function(x)
            .printLine(
                x["start"],
                x["end"],
                x["chrom"],
                x["score"],
                x["count"]
            )
    )
    sink()
}
