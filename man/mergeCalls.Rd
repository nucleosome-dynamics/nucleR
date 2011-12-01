\name{mergeCalls}
\alias{mergeCalls}
\title{
Automatic merging of overlapped nucleosome calls
}
\description{
	This function joints close nucleosome calls into one larger, fuzzy nucleosome.
}
\usage{
mergeCalls(calls, dyad.dist=70, discard.low=0.2, mc.cores=1)
}
\arguments{
  \item{calls}{
	\code{RangedData} with scored and ranged nucleosome calls from \code{peakScoring} or \code{peakDetection(..., score=TRUE)}.
}
	\item{dyad.dist}{
	Maximum distance dyad-to-dyad from which the peaks will be joined
}
  \item{discard.low}{
	Discard low covered calls (i.e. calls with \code{score_h < discard.low} will be discarded)
}
  \item{mc.cores}{
	Number of cores available to parallel data processing. 
}
}
\details{
	This functions looks for close (normally overlapped) calls and join those with their dyads separated less than \code{dyad.dist}. More than two reads can be joined in one single call.
	
	Joining is performed in chain, so if nucleosome call A is close to B and B is close to C, the final call will comprise the range A-B-C. The resulting scores (mixed, width, height) of the final joined call will be the average value of the individual scores.

	The parameter \code{discard.low} allows to ignore the small peaks that could be merged with larger ones, originating large calls. In the case that all of the overlapped reads in a given position have \code{score_h} less than \code{discard.low}, all of them will be selected instead of deleting that call. 
}
\value{
	\code{RangedData} with merged calls and the additional data column \code{nmerge}, with the count of how many original ranges are merged in the resulting range.
}
\author{
	Oscar Flores \email{oflores@mmb.pcb.ub.es}
}

\seealso{
	\code{\link{peakScoring}}
}
\keyword{manip}
\examples{
	#Generate a synthetic coverage map (assuming reads of 40bp and fragments of 130)
	map = syntheticNucMap(wp.num=20, fuz.num=20,  nuc.len=40, lin.len=130, rnd.seed=1)
	cover = filterFFT(coverage(map$syn.reads))

	#Find peaks over FFT filtered coverage
	calls = peakDetection(filterFFT(cover, pcKeepComp=0.02), width=130, score=TRUE)

	#Merge overlapped calls
	merged_calls = mergeCalls(calls)

	plotPeaks(merged_calls, cover)
}
