coverage.rpm <- function(data, scale=1e6, ...)
{
	#Use standard IRange coverage function
	cov = coverage(data, ...)

	#Change rle values to rpm. More efficient
	rle = RleList(lapply(cov, function(x) { x@values=(x@values/nrow(data))*scale; return(x) }))

	return(rle)
}
