coverage.rpm <- function(data)
{
	#Use standard IRange coverage function
	cov = coverage(data)

	#Change rle values to rpm. More efficient
	rle = RleList(lapply(cov, function(x) { x@values=(x@values/nrow(data))*1e6; return(x) }))

	return(rle)
}
