export.wig <- function (data, name, chrom = "", filepath = name)
{
    # Convert any list shaped object to simple list
    # (list, RleList, UncompressedRleList...)
    if(grepl("list", class(data), ignore.case=TRUE)) {
        data <- lapply(data, identity)
    } else {  #Otherwise create a list
        tmp <- data
        data <- list()
        data[[1]] <- tmp
    }

    #Assign chromosome names
    if(chrom == "") {
        if(is.null(names(data))) {
            stop("'chrom' parameter must be provided for unnamed objects")
        }
        chrom <- names(data)
    } else {
        names(data) <- chrom
    }

    for(chr in chrom) {
        file <- paste(filepath, chr, "wig", sep=".")
        values <- as.vector(data[[chr]])

        out <- file(description=file, open="wt", blocking=FALSE)

        tryCatch({
            cat(paste('track type=wiggle_0 name="', name, '"', sep=""),
                sep="\n", file=out, append=FALSE)
            cat(paste('fixedStep chrom=', chr, ' start=1 step=1', sep=""),
                sep="\n", file=out, append=TRUE)
            values <- round(values * 1000) / 1000  #Trim decimals
            cat(values, sep="\n", file=out, append=TRUE)
            close(out)
        }, error = function(e) {
            close(out)  #close the stream
            stop(e)
        })
    }
}
