#' Obtain and clean nucleosome positioning data from tiling array
#'
#' Process and transform the microarray data coming from tiling array
#' nucleosome positioning experiments.
#'
#' The processing of tiling arrays could be complicated as many types exists on
#' the market. This function deals ok with Affymetrix Tiling Arrays in yeast,
#' but hasn't been tested on other species or platforms.
#'
#' The main aim is convert the output of preprocessing steps (supplied by
#' third-parties packages) to a clean genome wide nucleosome occupancy profile.
#'
#' Tiling arrays doesn't use to provide a one-basepair resolution data, so one
#' gets one value per probe in the array, covering X basepairs and shifted
#' (tiled) Y basepairs respect the surrounding ones. So, one gets a piece of
#' information every Y basepairs.
#'
#' This function tries to convert this noisy, low resolution data, to a
#' one-basepair signal, which allows a fast recognition of nucleosomes without
#' using large and artificious statistical machinery as Hidden Markov Models
#' using posterionr noise cleaning process.
#'
#' As example, imagine your array has probes of 20mers and a tiling between
#' probes of 10bp. Starting at position 1 (covering coordinates from 1 to 20),
#' the next probe will be in position 10 (covering the coordinates 10 to 29).
#' This can be represented as two hybridization intensity values on coordinates
#' 1 and 10. This function will try to infer (using a lineal distribution) the
#' values from 2 to 9 using the existing values of probes in coordinate 1 and
#' coordinate 10.
#'
#' The tiling space between adjacent array probes could be not constant, or
#' could be also there are regions not covered in the used microarray. With the
#' function argument \code{inferLen} you can specify wich amout of space (in
#' basepairs) you allow to infer the non-present values.
#'
#' If at some point the range not covered (gap) between two adjacent probes of
#' the array is greater than \code{inferLen} value, then the coordinates
#' between these probes will be setted to NA.
#'
#' @param data \code{ExpressionSet} object wich contains the data of the tiling
#' array.
#' @param exprName Name of the \code{sample} in \code{ExpressionSet} which
#' contains the ratio between nucleosomal and genomic dna (if using
#' \code{Starr}, the \code{description} argument supplied to \code{getRatio}
#' function). If this name is not provided, it is assumed \code{data} has only
#' one column.
#' @param chrPattern Only chromosomes that contain \code{chrPattern} string
#' will be selected from ExpressionSet. Sometimes tilling arrays contain
#' control quality information that is imported as a chromosome. This allows
#' filtering it. If no value is supplied, all chromosomes will be used.
#' @param inferLen Maximum length (in basepairs) for allowing data gaps
#' inference. See \code{details} for further information.
#' @param mc.cores Number of cores available to parallel data processing.
#' @param quiet Avoid printing on going information (TRUE | FALSE)
#' @return RleList with the observed/inferred values for each coordinate.
#' @note This function should be suitable for all \code{data} objects of kind
#' ExpressionSet coding the annotations \code{"chr"} for chromosome and
#' \code{"pos"} for position (acccessible by \code{pData(data@featureData)})
#' and a expression value (accessible by \code{exprs(data)}).
#' @section Warning: This function could not cover all kind of arrays in the
#' market. This package assumes the data is processed and normalized prior this
#' processing, using standard microarray packages existing for R, like
#' \code{Starr}.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link[Biobase]{ExpressionSet}}, \code{\link[Starr]{getRatio}}
#' @keywords manip
#'
#' @examples
#' \dontrun{
#'     # Dataset cannot be provided for size restrictions
#'     # This is the code used to get the hybridization ratio with Starr from
#'     # CEL files
#'     library("Starr")
#'     TA_parsed <- readCelFile(
#'         BPMap, CELfiles, CELnames, CELtype, featureData=TRUE, log.it=TRUE
#'     )
#'     TA_loess <- normalize.Probes(TA_parsed, method="loess")
#'     TA_ratio <- getRatio(
#'         TA_loess, TA_loess$type=="IP", TA_loess$type=="CONTROL", "myRatio"
#'     )
#' 
#'     # From here, we use nucleR:
#' 
#'     # Preprocess the array, using the calculated ratio feature we named
#'     # "myRatio".
#' 
#'     # This will also select only those chromosomes with the pattern
#'     # "Sc:Oct_2003;chr", removing control data present in that tiling
#'     # array.
#' 
#'     # Finally, we allow that loci not covered by a prove being inferred
#'     # from adjacent ones, as far as they are separated by 50bp or less
#'     arr <- processTilingArray(
#'         TA_ratio, "myRatio", chrPattern="Sc:Oct_2003;chr", inferLen=50
#'     )
#' 
#'     # From here we can proceed with the analysis:
#'     arr_fft <- filterFFT(arr)
#'     arr_pea <- peakDetection(arr_fft)
#'     plotPeaks(arr_pea, arr_fft)
#' }
#'
#' @export processTilingArray
#'
processTilingArray <- function (data, exprName, chrPattern, inferLen = 50,
                                mc.cores = 1, quiet = FALSE)
{
    # require("Biobase")

    # Obtain annotation
    if (!quiet) {
        message(" * Parsing annotation")
    }
    df <- Biobase::pData(data@featureData)[, c("chr", "pos")]
    df$chr <- as.character(df$chr)

    # Add expression (intensities) values
    if (missing(exprName)) {
        exprName <- Biobase::sampleNames(data)[1]
    }
    if (!quiet) {
        message(paste(" * Using feature name:", exprName))
    }
    df$value <- Biobase::exprs(data)[, exprName]

    # Select names to keep
    if (!missing(chrPattern)) {
        selection <- unlist(sapply(
            unique(df$chr),
            function(x)
                grep(chrPattern, x, fixed=TRUE, value=TRUE))
        )
        df <- df[df$chr %in% selection, ]
    }

    .fillTheGap <- function(chr_name)
    {   # Function wich fills the gaps by chromosome
        # Select chromosome and sort
        dft <- df[df$chr == chr_name, c("pos", "value")]
        dft <- dft[order(dft$pos), ]

        # Calculate parameters for seq function

        v <- c(dft$pos[2:length(dft$pos)], dft$pos[length(dft$value)])
        dft$len <- v - dft$pos + 1
        dft$from <- dft$value
        dft$to <- c(
            dft$value[2:length(dft$value)],
            dft$value[length(dft$value)]
        )

        # Those values which are out of the gap range set them to a strange
        # value
        # Using some non numeric value has a huge impact on performance,
        # -999e99 is an odd value...
        dft[dft$len > inferLen, c("from", "to")] <- c(-999e99, -999e99)

        # Calculate the values sequence, but discard the first one of each
        # seq to not account it twice
        out <- apply(
            data.matrix(dft),
            1,
            function(x)
                seq(
                    from=x[["from"]],
                    to=x[["to"]],
                    length.out=x[["len"]]
                )[2:x[["len"]]]
        )
        out <- unlist(out)
        names(out) <- NULL
        out <- c(dft$from[1], out)  # Now add the first one hat is never counted
        out[out == - 999e99] <- NA  # Change the "strange" value for NA

        # If the starting point it's not 1, add NA's at the beggining
        if (dft$pos[[1]] > 1) {
            out <- c(rep(NA, dft$pos[[1]] - 1), out)
        }

        return(out)
    }

    # Get the names of all unique chromosomes to iterate on it
    chrs <- unique(df$chr)
    names(chrs) <- chrs

    if(!quiet) {
        message(" * Selected chromosomes:")
    }
    if(!quiet) {
        message(chrs)
    }

    message(" * Inferring missed values")
    res <- .xlapply(
        chrs,
        function(x) {
            if (!quiet) {
                message(paste("   ", x));
            }
            .fillTheGap(x)
        },
        mc.cores=mc.cores
    )

    return(res)
}
