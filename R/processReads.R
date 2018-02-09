#' Process reads from High-Troughtput Sequencing experiments
#'
#' This method allows the processment of NGS nucleosome reads from different
#' sources and a basic manipulation of them. The tasks includes the correction
#' of strand-specific single-end reads and the trimming of reads to a given
#' length.
#'
#' This function reads a \code{AlignedRead} or a \code{RangedData} object
#' containing the position, length and strand of the sequence reads.
#'
#' It allows the processment of both paired and single ended reads. In the case
#' of single end reads this function corrects the strand-specific mapping by
#' shifting plus strand reads and minus strand reads towards a middle position
#' where both strands are overlaped. This is done by accounting the expected
#' fragment length (\code{fragmentLen}).
#'
#' For paired end reads, mononucleosomal reads could extend more than expected
#' length due to mapping issues or experimental conditions. In this case, the
#' \code{fragmentLen} variable sets the threshold from which reads longer than
#' it should be ignored.
#'
#' If no value is supplied for \code{fragmentLen} it will be calculated
#' automatically (increasing the computing time) using \code{fragmentLenDetect}
#' with default parameters. Performance can be increased by tunning
#' \code{fragmentLenDetect} parameteres in a separated call and passing its
#' result as \code{fragmentLen} parameter.
#'
#' In some cases, could be useful trim the reads to a shorter length to improve
#' the detection of nucleosome dyads, easing its detection and automatic
#' positioning. The parameter \code{trim} allows the selection of how many
#' nucleotides select from each read.
#'
#' A special case for single-ended data is setting the \code{trim} to the same
#' value as \code{fragmentLen}, so the reads will be extended strand-wise
#' towards the 3' direction, creating an artificial map comparable with
#' paired-ended data. The same but opposite can be performed with paired-end
#' data, setting a \code{trim} value equal to the read length from paired
#' ended, so paired-ended data will look like single-ended.
#'
#' @param data Sequence reads objects, probably imported using other packages
#' as \code{ShortRead}. Allowed object types are \code{AlignedRead} and
#' \code{RangedData} with a \code{strand} attribute.
#' @param type Describes the type of reads. Values allowed are \code{single}
#' for single-ended reads and \code{paired} for paired-ended.
#' @param fragmentLen Expected original length of the sequenced fragments. See
#' details.
#' @param trim Length to trim the reads (or extend them if \code{trim} > read
#' length)
#' @param \dots Other parameters passed to \code{fragmentLenDetect} if no fixed
#' \code{fragmentLen} is given.
#'
#' @return \code{RangedData} containing the aligned/trimmed individual reads
#'
#' @note \strong{IMPORTANT}: this information is only used to correct possible
#' strand-specific mapping, this package doesn't link the two ends of paired
#' reads.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link[ShortRead]{AlignedRead}}, \code{\link{RangedData}},
#' \code{\link{fragmentLenDetect}}
#' @keywords manip
#' @rdname processReads
#'
#' @examples
#' # Load data
#' data(nucleosome_htseq)
#' 
#' # Process nucleosome reads, select only those shorter than 200bp
#' pr1 <- processReads(nucleosome_htseq, fragmentLen=200)
#' 
#' # Now process them, but picking only the 40 bases surrounding the dyad
#' pr2 <- processReads(nucleosome_htseq, fragmentLen=200, trim=40)
#' 
#' # Compare the results:
#' par(mfrow=c(2,1), mar=c(3,4,1,1))
#' plot(
#'     as.vector(coverage(pr1)[["chr1"]]), type="l",
#'     ylab="coverage (original)"
#' )
#' plot(
#'     as.vector(coverage(pr2)[["chr1"]]), type="l",
#'     ylab="coverage (trimmed)"
#' )
#'
#' @export
#'
setGeneric(
    "processReads",
    function(data, type="single", fragmentLen, trim, ...)
        standardGeneric("processReads")
)

#' @rdname processReads
setMethod(
    "processReads",
    signature(data="AlignedRead"),
    function (data, type="single", fragmentLen, trim, ...) {

        # require("ShortRead")
        if (missing(fragmentLen)) {
            if (type == "single") {
                message(
                    " * fragmentLen not provided for strand correction, ",
                    ", infering automatically..."
                )
                fragmentLen <- fragmentLenDetect(data, ...)
                message(paste(" * fragmentLen =", fragmentLen))
            } else {
                fragmentLen <- Inf  # Don't remove anything
            }
        }

        #######################################################################
        if (type == "single") {  # Special case for trim==fragmentLength
            if (!missing(trim)) {
                if (trim == fragmentLen) {
                    start <- ShortRead::position(data)
                    start[strand(data) == "-"] <- ShortRead::position(data) +
                        nchar(ShortRead::sread(data))
                    res <- IRanges::RangedData(
                        ranges=IRanges::IRanges(start=start, width=fragmentLen),
                        space=as.character(ShortRead::chromosome(data))
                    )
                    return (res)
                }
            }

            # If no trim restriction, use whole read length,
            # else use trim length
            sr_len <- nchar(ShortRead::sread(data))

            if (!missing(trim)) {
                sr_len <- trim
            }

            # Floor, then we will add 1 if needed
            shift  <- floor((fragmentLen - sr_len) / 2)

            strand <- strand(data)

            new_levels <- vector(mode="integer", length=3)
            new_levels[which(levels(strand) == "-")] <- -1
            new_levels[which(levels(strand) == "+")] <- +1
            new_levels[which(levels(strand) == "*")] <- 0
            levels(strand) <- new_levels

            strand <- as.numeric(levels(strand))[strand]

            # Correction for non rounded shifts in negative strand (ie 147/40)
            # If the fragment_length - read_length is odd and negative strand,
            # add 1 extra base
            correct <- rep(0, length(shift))
            correct[(strand == -1) & ((fragmentLen - sr_len) %% 2 == 1)] <- 1

            start <- ShortRead::position(data) + ((shift + correct) * strand)
            res <- IRanges::RangedData(
                ranges=IRanges::IRanges(start=start, width=sr_len),
                space=as.character(ShortRead::chromosome(data))
            )

        #######################################################################
        } else if (type == "paired") {
            sr_len <- nchar(ShortRead::sread(data))
            selection <- sr_len < fragmentLen
            space <- as.character(ShortRead::chromosome(data)[selection])

            # Normal start or trimmed version if requested
            start <- ShortRead::position(data)[selection]
            if (!missing(trim)) {
                start <- start + (sr_len[selection] / 2 - floor(trim / 2))
            }

            width <- sr_len[selection]
            if (!missing(trim)) {
                width <- trim
            }
            res <- IRanges::RangedData(
                ranges=IRanges::IRanges(start=start, width=width),
                space=space
            )

        #######################################################################
        } else {
            stop(paste(
                "type must be 'single' for single-ended data or 'paired' ",
                "for paired-ended"
            ))
        }

        return (res)
    }
)

#' @rdname processReads
setMethod(
    "processReads",
    signature(data="RangedData"),
    function (data, type="single", fragmentLen, trim, ...) {

        if (missing(fragmentLen)) {
            if (type == "single") {
                message(" * fragmentLen not provided for strand correction, ",
                        "infering automatically...")
                fragmentLen <- fragmentLenDetect(data,...)
                message(paste(" * fragmentLen =", fragmentLen))
            } else {
                fragmentLen <- Inf  # Don't remove anything
            }
        }

        #######################################################################
        if (type == "single") {  # Special case for trim==fragmentLength
            if (!missing(trim)) {
                if (trim == fragmentLen) {
                    start <- start(data)
                    negStrand <- data$strand == "-"
                    start[negStrand] <- IRanges::end(data)[negStrand] - fragmentLen
                    res <- IRanges::RangedData(
                        ranges=IRanges::IRanges(start=start, width=fragmentLen),
                        space=S4Vectors::space(data)
                    )
                    return (res)
                }
            }

            # If no trim restriction, use whole read length,
            # else use trim length
            sr_len <- width(data)
            if (!missing(trim)) {
                sr_len <- trim
            }

            # Floor, then we will add 1 if needed
            shift <- floor((fragmentLen - sr_len) / 2)

            strand <- rep(1, nrow(data))
            strand[data$strand == "-"] <- -1

            # Correction for non rounded shifts in negative strand (ie 147/40)
            # If the fragment_length - read_length is odd and negative strand,
            # add 1 extra base
            correct <- rep(0, length(strand))
            correct[(strand == -1) & ((fragmentLen - sr_len) %% 2 == 1)] <- 1

            start <- start(data) + ((shift + correct) * strand)
            res <- IRanges::RangedData(
                ranges=IRanges::IRanges(start=start, width=sr_len),
                space=S4Vectors::space(data)
            )

        #######################################################################
        } else if (type=="paired") {
            data <- data[width(data) < fragmentLen, ]
            sr_len <- width(data)

            # Normal start or trimmed version if requested
            start <- start(data)
            if (!missing(trim)) {
                start <- start + ((sr_len / 2) - floor(trim / 2))
            }

            width <- sr_len

            if (!missing(trim)) {
                width <- trim
            }

            res <- IRanges::RangedData(
                ranges=IRanges::IRanges(start=start, width=width),
                space=S4Vectors::space(data)
            )
        #######################################################################
        } else {
            stop(paste(
                "type must be 'single' for single-ended data or",
                "'paired' for paired-ended")
            )
        }

        return (res)
    }
)

#' @rdname processReads
setMethod(
    "processReads",
    signature(data="GRanges"),
    function (data, type="single", fragmentLen, trim, ...) {

        if (missing(fragmentLen)) {
            if (type == "single") {
                message(
                    " * fragmentLen not provided for strand correction, ",
                    "infering automatically..."
                )
                fragmentLen <- fragmentLenDetect(data,...)
                message(paste(" * fragmentLen =", fragmentLen))
            } else {
                fragmentLen <- Inf  # Don't remove anything
            }
        }

        #######################################################################
        if (type == "single") {  # Special case for trim==fragmentLength
            if (!missing(trim)) {
                if (trim == fragmentLen) {
                    start <- start(data)
                    negStrand <- data$strand == "-"
                    start[negStrand] <- IRanges::end(data)[negStrand] - fragmentLen
                    res <- IRanges::RangedData(
                        ranges=IRanges::IRanges(start=start, width=fragmentLen),
                        space=S4Vectors::space(data)
                    )
                    return (res)
                }
            }

            # If no trim restriction, use whole read length,
            # else use trim length
            if (missing(trim)) {
                sr_len <- width(data)
            } else {
                sr_len <- trim
            }

            # Floor, then we will add 1 if needed
            shift <- floor((fragmentLen - sr_len) / 2)

            strand <- rep(1, length(data))
            strand[data$strand == "-"] <- -1

            # Correction for non rounded shifts in negative strand (ie 147/40)
            # If the fragment_length - read_length is odd and negative strand,
            # add 1 extra base
            correct <- rep(0, length(strand))
            correct[(strand == -1) & ((fragmentLen - sr_len) %% 2 == 1)] <- 1

            start <- start(data) + ((shift+correct) * strand)

            res <- GenomicRanges::GRanges(
                ranges=IRanges::IRanges(start=start, width=sr_len),
                seqnames=GenomeInfoDb::seqnames(data)
            )

        #######################################################################
        } else if (type == "paired") {
            data <- data[width(data) < fragmentLen, ]
            sr_len <- width(data)

            # Normal start or trimmed version if requested
            start <- start(data)
            if (!missing(trim)) {
                start <- start + ((sr_len / 2) - floor(trim / 2))
            }

            width <- sr_len

            if (!missing(trim)) {
                width <- trim
            }

            res <- GenomicRanges::GRanges(
                ranges=IRanges::IRanges(start=start, width=width),
                seqnames=GenomeInfoDb::seqnames(data)
            )
        #######################################################################
        } else {
            stop(paste(
                "type must be 'single' for single-ended data or",
                "'paired' for paired-ended"
            ))
        }

        return (res)
    }
)
