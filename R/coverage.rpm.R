#' Coverage calculation and normalization to reads per million (rpm)
#'
#' Calculates the coverage values from a [GenomicRanges::GRanges] or
#' [IRanges::IRanges] object normalized to reads per million, allowing the
#' comparison of experiments with a different absolut number of reads.
#'
#' @param data [GenomicRanges::GenomicRanges] or [IRanges::IRanges] with the
#'   reads information
#' @param scale By default, a million (1e6), but you could change this value
#'    for abnormal high or low amount of reads.
#' @param \dots Additional arguments to be passed to `coverage` function
#'
#' @return `RleList` object with the coverage objects
#'
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso [processReads()], [IRanges::coverage()]
#' @keywords manip
#' @rdname coverage.rpm
#'
#' @examples
#' # Load the example dataset and get the coverage
#' data(nucleosome_htseq)
#' cov <- coverage.rpm(nucleosome_htseq)
#' print(cov)
#' # Plot it
#' library(ggplot2)
#' cover <- as.vector(cov[["chr1"]])
#' qplot(seq_along(cover), cover, geom="line", ylab="coverage",
#'   xlab="position")
#'
#' @export
#'
setGeneric(
    "coverage.rpm",
    function (data, scale=1e6, ...)
        standardGeneric("coverage.rpm")
)

#' @rdname coverage.rpm
#' @importFrom IRanges RleList
#' @importMethodsFrom IRanges coverage
setMethod(
    "coverage.rpm",
    signature(data="GRanges"),
    function (data, scale=1e6, ...)
        RleList(lapply(
            coverage(data, ...),
            function (x) x / length(data) * scale
        ), compress=FALSE)
)

#' @rdname coverage.rpm
setMethod(
    "coverage.rpm",
    signature(data="CompressedGRangesList"),
    function (data, scale=1e6, ...)
        coverage.rpm(do.call(`c`, unname(data)), scale=scale, ...)
)

#' @rdname coverage.rpm
#' @importMethodsFrom IRanges coverage
setMethod(
    "coverage.rpm",
    signature(data="IRanges"),
    function (data, scale=1e6, ...)
        coverage(data) / length(data) * scale
)

#' @rdname coverage.rpm
#' @importFrom IRanges RleList
#' @importMethodsFrom IRanges coverage
setMethod(
    "coverage.rpm",
    signature(data="RangedData"),
    function (data, scale=1e6, ...)
        RleList(lapply(
            coverage(data, ...),
            function (x) x / nrow(data) * scale
        ), compress=FALSE)
)
