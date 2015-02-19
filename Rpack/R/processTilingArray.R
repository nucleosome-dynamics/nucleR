processTilingArray <- function(data, exprName, chrPattern, inferLen=50, mc.cores=1, quiet=FALSE) {

	require("Biobase")

  #Check if multicore is supported or set to 1
  mc.cores = .check.mc(mc.cores) 

	#Obtain annotation
	if(!quiet) message(" * Parsing annotation")
	df = pData(data@featureData)[,c("chr", "pos")]
	df$chr = as.character(df$chr)
	
	#Add expression (intensities) values
	if(missing(exprName)) exprName = sampleNames(data)[1]
	if(!quiet) message(paste(" * Using feature name:", exprName))
	df$value = exprs(data)[,exprName]

	#Select names to keep
	if(!missing(chrPattern))
	{
		selection = unlist(sapply(unique(df$chr), function(x) grep(chrPattern, x, fixed=TRUE, value=TRUE)))
		df = df[df$chr %in% selection,]
	}

	#Create the function wich fills the gaps by chromosome
	.fillTheGap<-function(chr_name)
	{
		#Select chromosome and sort
		dft = df[df$chr==chr_name, c("pos", "value")]
		dft = dft[order(dft$pos),]
	
		#Calculate parameters for seq function
		dft$len   = (c(dft$pos[2:length(dft$pos)], dft$pos[length(dft$value)]) - dft$pos) + 1
		dft$from  = dft$value
		dft$to    = c(dft$value[2:length(dft$value)], dft$value[length(dft$value)])

		#Those values which are out of the gap range set them to a strange value
		#Using some non numeric value has a huge impact on performance, -999e99 is a odd value...
		dft[dft$len > inferLen, c("from","to")]	= c(-999e99,-999e99)

		#Calculate the values sequence, but discard the first one of each seq to not account it twice
		out = apply(data.matrix(dft), 1, function(x)
								 seq(from=x[["from"]], to=x[["to"]], length.out=x[["len"]])[2:x[["len"]]]) 
		out = unlist(out)
		names(out) = NULL
		out = c(dft$from[1], out) #Now add the first one hat is never counted
		out[out == -999e99] = NA # Change the "strange" value for NA

		#If the starting point it's not 1, add NA's at the beggining
		if(dft$pos[[1]] > 1) out = c(rep(NA, dft$pos[[1]] - 1), out)

		return(out)
	}

	#Get the names of all unique chromosomes to iterate on it
	chrs = unique(df$chr)
	names(chrs) = chrs

	if(!quiet) message(" * Selected chromosomes:")
	if(!quiet) message(chrs)

	if(mc.cores > 1)
	{
		message(" * Inferring missed values (multicore version)")
		res = mclapply(chrs, .fillTheGap, mc.cores=mc.cores)
	}
	else
	{
		message(" * Inferring missed values (single core version): ")
		res = lapply(chrs, function(x){ if(!quiet) message(paste("   ", x)); .fillTheGap(x) })
	}

	return(res)
}

