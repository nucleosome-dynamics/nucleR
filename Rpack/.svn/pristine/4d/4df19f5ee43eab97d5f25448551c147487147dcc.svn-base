\name{pcKeepCompDetect}
\alias{pcKeepCompDetect}
\title{
	Auto detection of a fitted \code{pcKeepComp} param for filterFFT function
}
\description{
	This function tries to obtain the minimum number of components needed in a FFT filter to
	achieve or get as close as possible to a given correlation value. Usually you don't need
	to call directly this function, is used in \code{filterFFT} by default.
}
\usage{
pcKeepCompDetect(data, pc.min=0.01, pc.max=0.1, max.iter=20, verbose=FALSE, 
			cor.target=0.98, cor.tol=1e-3, smpl.num=25, smpl.min.size=2^10, smpl.max.size=2^14)
}
\arguments{
  \item{data}{
	Numeric vector to be filtered
}
  \item{pc.min, pc.max}{
	Range of allowed values for pcKeepComp (minimum and maximum), in the range 0:1.
}
	\item{max.iter}{
	Maximum number of iterations
}
	\item{verbose}{
	Extra information (debug)
}
	\item{cor.target}{
	Target correlation between the filtered and the original profiles. A value around 0.99 is recommeded for Next Generation Sequencing data and around 0.7 for Tiling Arrays.
}
	\item{cor.tol}{
	Tolerance allowed between the obtained correlation an the target one.
}
  \item{smpl.num}{
	If \code{data} is a large vector, some samples from the vector will be used instead the whole dataset. This parameters tells the number of samples to pick.
}
	\item{smpl.min.size, smpl.max.size}{
	Minimum and maximum size of the samples. This is used for selection and sub-selection of ranges with meaningful values (i,e, different from 0 and NA). Power of 2 values are recommended, despite non-mandatory.
}
	\item{\dots}{
	Parameters to be pass to \code{autoPcKeepComp}
}
}
\details{
	This function predicts a suitable \code{pcKeepComp} value for \code{filterFFT} function. This is the recommended amount of components (in percentage) to keep in the \code{filterFFT} function to obtain a correlation of (or near of) \code{cor.target}.
	
	The search starts from two given values \code{pc.min, pc.max} and uses linial interpolation to quickly reach a value that gives a corelation between the filtered and the original near \code{cor.target} within the specified tolerance \code{cor.tol}.

	To allow a quick detection without an exhaustive search, this function uses a subset of the data by randomly sampling those regions with meaningful coverage values (i,e, different from 0 or NA) larger than \code{smpl.min.size}. If it's not possible to obtain \code{smpl.max.size} from this region (this could be due to flanking 0's, for example) at least \code{smpl.min.size} will be used to check correlation. Mean correlation between all sampled regions is used to test the performance of the pcKeepComp parameter.

	If the number of meaningful bases in \code{data} is less than \code{smpl.min.size * (smpl.num/2)} all the \code{data} vector will be used instead of using sampling.
}
\value{
	Fitted \code{pcKeepComp} value
}
\author{
  Oscar Flores \email{oflores@mmb.pcb.ub.es}, David Rosell \email{david.rosell@irbbarcelona.org}
}
\keyword{ attribute }
\examples{

	#Load dataset
	data(nucleosome_htseq)
	data = as.vector(coverage.rpm(nucleosome_htseq)[[1]])

	#Get recommended pcKeepComp value
	pckeepcomp = pcKeepCompDetect(data, cor.target=0.99)
	print(pckeepcomp)

	#call filterFFT
	f1 = filterFFT(data, pcKeepComp=pckeepcomp)

	#Also this can be called directly
	f2 = filterFFT(data, pcKeepComp="auto", cor.target=0.99)

	#Plot
	plot(data[1:2000], col="black", type="l", lwd=2)
	lines(f1[1:2000], col="red", lwd=2)
	lines(f2[1:2000], col="blue", lwd=2, lty=2)
	legend("bottom", c("original", "two calls", "one call"), col=c("black", "red", "blue"), lty=c(1,1,2), horiz=TRUE, bty="n")
}
