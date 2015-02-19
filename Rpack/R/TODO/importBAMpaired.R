library(multicore)
library(ShortRead)
library(Rsamtools)

#VARS
INPUT_DIR="01_bowtie"
OUTPUT_DIR="02_reads"
CORES=8   #Each core will process 1 BAM file, beware of the disk access

system(paste("mkdir", OUTPUT_DIR))

############################################################################################

#Binary conversion 
int2base <- function(x, b=2){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0))) print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1 
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b
}

#SAM/BAM flag matrix
bamFlagMatrix <- function(flags){
  bin = int2base(flags)
  n = ncol(bin)
  colnames(bin) = c(rev(names(formals(scanBamFlag))[1:n]))
  return(bin)
}

#SAM/BAM flag summary
bamFlagSummary <- function(flags){
  tab = table(flags)
  fla = as.integer(names(tab))
  bin = int2base(fla)
  rownames(bin) = fla
  bin = data.frame(bin)
  bin$count = tab
  n = ncol(bin)
  names = c(rev(names(formals(scanBamFlag))[1:n-1]), "count")
  names(bin) = names
  return(bin[,c(n:1,n)])
}

#This is the main function to process one given file, allows paralelism
process <- function(file)
{
	message("** ", file, " -> Reading")
	#Read BAM file (only one access to disk, intended for Shared Memory)
	what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mrnm", "mpos")
	bam  = scanBam(file=file, param=ScanBamParam(what=what))[[1]] #Only read one file
	
  message("** ", file, " -> Processing flags")
	#We will process the flags in R (an alternative is multiple scanBam calls...)
	flags = bamFlagMatrix(bam$flag)

  message("** ", file, " -> Remove multiple matches")
	#Remove multiple mappings
	mmids = unique(bam$qname[flags[,"isPrimaryRead"]==1]) #isPrimaryRead=1 means the read is repeated...
	mmids = which(bam$qname %in% mmids)
	bam2   = lapply(bam, "[", -mmids)      #Remove repeated IDs from bam object and flags
	flags = flags[-mmids,] 

	#######################################################################################
	#Process + strand
  message("** ", file, " -> Processing + strand")
	pos_1 = flags[,"isPaired"] & flags[,"isProperPair"] & 
					flags[,"isFirstMateRead"] & !flags[,"isMinusStrand"]

	pos_2 = flags[,"isPaired"] & flags[,"isProperPair"] & 
          flags[,"isSecondMateRead"] & flags[,"isMinusStrand"]

	reads1 = lapply(bam2, "[", pos_1) 
	reads2 = lapply(bam2, "[", pos_2)

	#Consistency check
	test = reads1$mpos == reads2$pos & reads2$mpos == reads1$pos & reads1$rname == reads2$rname
	if(!all(test)) stop("ERROR: Mate selection for + strand is invalid")
	
	#Create reads for the positive strand
	pos_ranges = IRanges(start=reads1$pos, end=reads2$pos+reads2$qwidth-1)
	pos_names  = as.character(reads1$rname)

  #######################################################################################
  #Process - strand
  message("** ", file, " -> Processing - strand")
  pos_1 = flags[,"isPaired"] & flags[,"isProperPair"] &
          flags[,"isFirstMateRead"] & flags[,"isMinusStrand"]

  pos_2 = flags[,"isPaired"] & flags[,"isProperPair"] &
          flags[,"isSecondMateRead"] & !flags[,"isMinusStrand"]

  reads1 = lapply(bam2, "[", pos_1)
  reads2 = lapply(bam2, "[", pos_2)

  #Consistency check
  test = reads1$mpos == reads2$pos & reads2$mpos == reads1$pos & reads1$rname == reads2$rname
  if(!all(test)) stop("ERROR: Mate selection for - strand is invalid")
  
  #Create reads for the positive strand
	neg_ranges = IRanges(start=reads2$pos, end=reads1$pos+reads1$qwidth-1)
	neg_names  = as.character(reads1$rname)


	############ Join
  message("** ", file, " -> Joining and sorting")
	ranges = c(pos_ranges, neg_ranges)
	ord = order(ranges)
	ranges = ranges[ord]
	spaces = c(pos_names, neg_names)
	spaces = spaces[ord]

	reads = RangedData(space=spaces, ranges=ranges)

	############ Save
  message("** ", file, " -> Saving")
	base = sub(".bam", "", basename(file))
	save(reads, file=paste(OUTPUT_DIR, "/", base, ".RData", sep=""))
}	

#######################################################################
#Main call

#Get filenames and remove extension
exps = dir(INPUT_DIR, pattern="bam$", full.names=TRUE) #loads all files ending with *.bam
res = mclapply(exps, process, mc.cores=CORES)
