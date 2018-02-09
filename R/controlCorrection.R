#' Correct experimental profiles with control sample
#'
#' This function allows the correction of experimental coverage profiles
#' (usually MNase digested nucleosomal DNAs in this library) with control
#' samples (usually naked DNA sample digested with MNase). This is useful to
#' correct MNase biase.
#'
#' This substracts the enrichment in the control sample respect it's mean from
#' the experimental profile.
#'
#' This is useful for examinating the effect of the MNase digestion in
#' nucleosome experiments using a nucleosomal DNA and a genomic (naked) DNA
#' sample. Notice that genomic DNA samples cannot be strand-corrected using
#' single end data, so only paired end controls are useful for this proupose,
#' despite they can be compared against extended nucleosomal DNA single end
#' reads. Furthermore, both datasets must be converted to reads per milion.
#'
#' This process dificults the nucleosome positioning due the lower sharpness of
#' the peaks, but allows a complementary study of the MNase digestion effect.
#'
#' @param exp,ctr Comparable experimental and control samples (this means same
#' format and equivalent preprocessment)
#' @param mc.cores Number of cores available for parallel list processing
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return Corrected experimental profile
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @keywords manip
#' @rdname controlCorrection
#'
#' @examples
#' map = syntheticNucMap(as.ratio=TRUE)
#' exp = coverage(map$syn.reads)
#' ctr = coverage(map$ctr.reads)
#' corrected = controlCorrection(exp, ctr)
#'
#' @export
#'
setGeneric(
    "controlCorrection",
    function(exp, ctr, ...)
        standardGeneric("controlCorrection")
)

#' @rdname controlCorrection
setMethod(
    "controlCorrection",
    signature(exp="SimpleRleList"),
    function (exp, ctr, mc.cores=1) {

        if (class(exp) != class(ctr)) {
            stop("'exp' and 'ctr' classes must be equal")
        }

        if (length(which(!(names(exp) %in% names(ctr)))) != 0 |
            length(which(!(names(ctr) %in% names(exp)))) != 0) {
            stop("names should be equal in both datasets")
        }

        res <- .xlapply(
            names(exp),
            function(chr)
                controlCorrection(exp[[chr]], ctr[[chr]]),
            mc.cores=mc.cores
        )
        names(res) <- names(exp)
        return(IRanges::RleList(res, compress=FALSE))
    }
)

#' @rdname controlCorrection
setMethod(
    "controlCorrection",
    signature(exp="Rle"),
    function (exp, ctr) {
        if (class(exp) != class(ctr)) {
            stop("'exp' and 'ctr' classes must be equal")
        }
        return(S4Vectors::Rle(controlCorrection(as.vector(exp), as.vector(ctr))))
    }
)

#' @rdname controlCorrection
setMethod(
    "controlCorrection",
    signature(exp="list"),
    function (exp, ctr, mc.cores=1) {

        if (class(exp) != class(ctr)) {
            stop("'exp' and 'ctr' classes must be equal")
        }

        if (length(which(!(names(exp) %in% names(ctr)))) != 0 |
            length(which(!(names(ctr) %in% names(exp)))) != 0) {
            stop("names should be equal in both datasets")
        }

        res <- .xlapply(
            names(exp),
            function(chr)
                controlCorrection(exp[[chr]], ctr[[chr]]),
            mc.cores=mc.cores,
        )
        names(res) <- names(exp)
        return(res)
    }
)

#' @rdname controlCorrection
setMethod(
    "controlCorrection",
    signature(exp="numeric"),
    function (exp, ctr)  {

        if (class(exp) != class(ctr)) {
            stop("'exp' and 'ctr' classes must be equal")
        }

        res <- exp - (ctr - mean(ctr[ctr != 0]))
        res[res  < 0] <- 0
        return(res)
    }
)
