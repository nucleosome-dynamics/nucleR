#!/usr/bin/env Rscript

.export.bed <- function(df, name, desc, filename)
{
    hasScore <- "score" %in% names(df)
    tscore <- ifelse(hasScore, "useScore=1", "useScore=0")

    tdesc <- paste('description="', desc, '"', sep="")
    tname <- paste('name="', name, '"', sep="")

    header <- paste("track", tname, tdesc, tscore, sep="\t")


    .printLine <- function(start, end, chr, score, num)
        cat(paste(chr, start, end, paste("nucl", as.numeric(num), sep=""),
                  score, sep="\t"),
            sep="\n")

    sink(paste(filename, "bed", sep="."))
    cat(header, sep="\n")
    res <- apply(
        df,
        1,
        function(x)
            .printLine(x["start"], x["end"], x["chrom"], x["score"], x["count"])
    )
    sink()
}

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

        .export.bed(df=dd, name=name, desc=desc,
                    filename=paste(filepath, chrom, sep="."))
    }
)

setMethod(
    "export.bed",
    signature(ranges="CompressedIRangesList"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
              splitByChrom=TRUE) {

        if (splitByChrom) {
            for(chr in names(ranges)) {
                export.bed(ranges=ranges[[chr]], score=score[[chr]],
                           chrom=chr, name=name, desc=desc, filepath=filepath,
                           splitByChrom=FALSE)
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

setMethod(
    "export.bed",
    signature(ranges="RangedData"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
              splitByChrom=TRUE) {

        if (splitByChrom) {
            for(chr in names(ranges)) {
                export.bed(ranges=ranges[chr], chrom=chr, name=name, desc=desc,
                           filepath=paste(filepath, chr, sep="."),
                           splitByChrom=FALSE)
            }
        } else {
            dd <- as.data.frame(ranges)
            dd$chrom <- dd$space
            dd$count <- seq(1, nrow(dd))

            .export.bed(df=dd, name=name, desc=desc, filename=filepath)
        }
    }
)

setMethod(
    "export.bed",
    signature(ranges="GRanges"),
    function (ranges, score=NULL, name, desc=name, filepath=name,
              splitByChrom=TRUE) {

        if (splitByChrom) {
            #for(chr in IRanges::levels(seqnames(ranges))) {
            for(chr in levels(seqnames(ranges))) {
                export.bed(ranges[seqnames(ranges) == chr],
                           chrom=chr, name=name, desc=desc,
                           filepath=paste(filepath, chr, sep="."),
                           splitByChrom=FALSE)
            }
        } else {
            dd <- as.data.frame(ranges)
            dd$chrom <- dd$seqnames
            dd$count <- seq(1, nrow(dd))

            .export.bed(df=dd, name=name, desc=desc, filename=filepath)
        }
    }
)
