\name{controlCorrection}
\alias{controlCorrection}
\alias{controlCorrection,SimpleRleList-method}
\alias{controlCorrection,Rle-method}
\alias{controlCorrection,list-method}
\alias{controlCorrection,numeric-method}
\title{
	Correct experimental profiles with control sample
}
\description{
	This function allows the correction of experimental coverage profiles (usually MNase digested nucleosomal DNAs in this library) with control samples (usually naked DNA sample digested with MNase). This is useful to correct MNase biase.
}
\usage{
\S4method{controlCorrection}{SimpleRleList}(exp, ctr, mc.cores=1)
\S4method{controlCorrection}{Rle}(exp, ctr)
\S4method{controlCorrection}{list}(exp, ctr, mc.cores=1)
\S4method{controlCorrection}{numeric}(exp, ctr)
}
\arguments{
  \item{exp, ctr}{
	Comparable experimental and control samples (this means same format and equivalent preprocessment)
}
	\item{mc.cores}{
	Number of cores available for parallel list processing
}
}
\details{
	This substracts the enrichment in the control sample respect it's mean from the experimental profile.

	This is useful for examinating the effect of the MNase digestion in nucleosome experiments using a nucleosomal DNA and a genomic (naked) DNA sample. Notice that genomic DNA samples cannot be strand-corrected using single end data, so only paired end controls are useful for this proupose, despite they can be compared against extended nucleosomal DNA single end reads. Furthermore, both datasets must be converted to reads per milion.
	
	This process dificults the nucleosome positioning due the lower sharpness of the peaks, but allows a complementary study of the MNase digestion effect.
}
\value{
	Corrected experimental profile
}
\author{
	Oscar Flores \email{oflores@mmb.pcb.ub.es}
}
\keyword{ manip }
\examples{

		#Toy example
		map = syntheticNucMap(as.ratio=TRUE)

		exp = coverage(map$syn.reads)	
		ctr = coverage(map$ctr.reads)
		
		corrected = controlCorrection(exp, ctr)
}
