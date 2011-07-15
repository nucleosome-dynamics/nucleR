export.wig <- function(data, name, chrom="", filepath=name)
{
	data = lapply(data, identity)
	if(chrom == "") chrom = names(data)

	for(chr in chrom)
	{
		file = paste(filepath, chr, "wig", sep=".")
		values = as.vector(data[[chr]])

	  sink(file)
		tryCatch({
		  cat(paste('track type=wiggle_0 name="', name, '"', sep=""), sep="\n")
  		cat(paste('fixedStep chrom=', chr, ' start=1 step=1', sep=""), sep="\n")
	  	cat(values, sep="\n")
	  	sink()
		}, error = function(e) {
			sink()  #close the stream
			stop(e)
		})
	}
}
