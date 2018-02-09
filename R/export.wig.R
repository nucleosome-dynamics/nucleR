#' Export values in WIG format
#'
#' Export coverage/intensity values in WIG format, compatible with UCSC genome
#' browser, IGB and others.
#'
#' @param data Coverage/intensity values (numeric vector)
#' @param name Name of the track
#' @param chrom Information about chromosome if not inferrable from \code{data}
#' (only for numeric vectors)
#' @param filepath Filepath where to save the object. Chromosome name and "wig"
#' extension will be automatically added
#' @return (none)
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @references WIG format specification:
#' http://genome.ucsc.edu/FAQ/FAQformat#format6
#' @keywords file
#' @examples
#'
#' # Load data
#' data(nucleosome_htseq)
#' cover <- coverage.rpm(nucleosome_htseq)
#'
#' # Create wig file
#' export.wig(cover, name="example_track")
#'
#' # This would create the file "example_track.chr1.wig" with:
#'
#' # track type=wiggle_0 name="example_track"
#' # fixedStep chrom=chr1 start=1 step=1
#' # 55.55247
#' # 55.55247
#' # 55.55247
#' # 277.7623
#' # 388.8673
#' # ...
#'
#' @export export.wig
#
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
