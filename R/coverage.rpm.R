#' Coverage calculation and normalization to reads per million (rpm)
#'
#' Calculates the coverage values from a \code{RangedData} object (or anything
#' with a defined \code{coverage} function associated) and returns the coverage
#' normalized to reads per million, allowing the comparison of experiments with
#' a different absolut number of reads.
#'
#' @param data \code{RangedData} (or compatible) with the reads information
#' @param scale By default, a million (1e6), but you could change this value
#' for abnormal high or low amount of reads.
#' @param \dots Additional arguments to be passed to \code{coverage} function
#'
#' @return \code{RleList} object with the coverage objects
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{processReads}}, \code{\link{coverage}}
#' @keywords manip
#' @rdname coverage.rpm
#'
#' @examples
#' #Load the example dataset and get the coverage
#' data(nucleosome_htseq)
#' cov = coverage.rpm(nucleosome_htseq)
#' print(cov)
#' #Plot it
#' plot(as.vector(cov[["chr1"]]), type="l", ylab="coverage", xlab="position")
#'
#' @export
#'
setGeneric(
    "coverage.rpm",
    function (data, scale=1e6, ...)
        standardGeneric("coverage.rpm")
)

#' @rdname coverage.rpm
setMethod(
    "coverage.rpm",
    signature(data="GRanges"),
    function(data, scale=1e6, ...)
        IRanges::RleList(lapply(
            IRanges::coverage(data, ...),
            function(x) x / length(data) * scale
        ), compress=FALSE)
)

#' @rdname coverage.rpm
setMethod(
    "coverage.rpm",
    signature(data="RangedData"),
    function(data, scale=1e6, ...)
        IRanges::RleList(lapply(
            IRanges::coverage(data, ...),
            function(x) x / nrow(data) * scale
        ), compress=FALSE)
)
