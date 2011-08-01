export.wig <- function(data, name, chrom="", filepath=name)
{
	#Convert any list shaped object to simple list (list, RleList, UncompressedRleList...)
	if(grepl(class(data), "list", ignore.case=TRUE))
	{
		data = lapply(data, identity)
	}
	else{ #Otherwise create a list
		tmp = data
		data = list()
		data[[1]] = tmp
	}
	
	#Assign chromosome names
	if(chrom == "")
	{
		if(is.null(names(data))) stop("'chrom' parameter must be provided for unnamed objects")
		chrom = names(data)
	}else{
		names(data) = chrom
	}
	
	for(chr in chrom)
	{
		file = paste(filepath, chr, "wig", sep=".")
		values = as.vector(data[[chr]])

	  sink(file)
		tryCatch({
			cat(paste('track type=wiggle_0 name="', name, '"', sep=""), sep="\n")
			cat(paste('fixedStep chrom=', chr, ' start=1 step=1', sep=""), sep="\n")
			cat(format(values, nsmall=4), sep="\n")
			sink()
		}, error = function(e) {
			sink()  #close the stream
			stop(e)
		})
	}
}
