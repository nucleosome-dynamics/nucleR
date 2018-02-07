setMethod(
    "coverage.rpm",
    signature(data="GRanges"),
    function(data, scale=1e6, ...)
        RleList(lapply(
            coverage(data, ...),
            function(x) x / length(data) * scale
        ), compress=FALSE)
)

setMethod(
    "coverage.rpm",
    signature(data="RangedData"),
    function(data, scale=1e6, ...)
        RleList(lapply(
            coverage(data, ...),
            function(x) x / nrow(data) * scale
        ), compress=FALSE)
)
