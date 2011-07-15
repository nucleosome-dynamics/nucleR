\name{peakScoring}
\alias{peakScoring}
\alias{peakScoring,list-method}
\alias{peakScoring,IRangesList-method}
\alias{peakScoring,numeric-method}
\alias{peakScoring,IRanges-method}
\title{
	Peak scoring function
}
\description{
	Scores peaks detected with function \code{peakDetection} according the height and the sharpness of the peak. This function can be called automatically from \code{peakDetection} if \code{score=TRUE}.
}
\usage{
	\S4method{peakScoring}{numeric}(peaks, data, threshold=0.25)	
	\S4method{peakScoring}{list}(peaks, data, threshold=0.25, mc.cores=1)
	\S4method{peakScoring}{IRanges}(peaks, data, threshold=0.25, dyad.length=38)
	\S4method{peakScoring}{IRangesList}(peaks, data, threshold=0.25, dyad.length=38, mc.cores=1)
}
\arguments{
  \item{peaks}{
	The identified peaks resulting from \code{peakDetection}. Could be a \code{numeric} vector with the position of the peaks, or a \code{IRanges} object with the extended range of the peak.
	For both types, list support is implemented as a \code{numeric} list or a \code{IRangesList}
}
	\item{data}{
	Data of nucleosome coverage or intensites.
}
	\item{threshold}{
	The non-default \code{threshold} previously used in \code{peakDetection} function, if applicable.
}
	\item{dyad.length}{
	How many bases account in the nucleosome dyad for sharpness description. If working with NGS data, works best with the reads width value for single-ended data or the \code{trim} value given to the \code{processReads} function.
}
  \item{mc.cores}{
	If input is a \code{list} or \code{IRangeList}, and multiple cores support is available, the maximum number of cores for parallel processing
}
}
\details{
	This function scores each previously identified peak according its height and sharpness. 

	The height score tells how large is a peak, higher means more coverage or intensity, so better positioned nucleosome. This score is obtained by checking the observed peak value in a Normal distribution with the mean and sd of \code{data}. This value is between 0 and 1.

	The width score is a mesure of how sharp is a peak. With a NGS coverage in mind, a perfect phased (well-positioned) nucleosome is this that starts and ends exactly in the same place many times. The shape of this ideal peak will be a rectangular shape of the lenght of the read. A wider top of a peak could indicate fuzzyness. The parameter \code{dyad.length} tells how long should be the "flat" region of an ideal peak. The optimum value for this parameter is the lenght of the read in single-ended data or the \code{trim} value of the function \code{processReads}. For Tiling Array, the default value should be fine.
	This score is obtained calculating the ratio between the mean of the nucleosome scope (the one provided by range in the elements of \code{peaks}) and the \code{dyad.length} central bases. This value is normalized between 0 and 1.

	For punctual, single points peaks (provided by \code{numeric} vector or list as \code{peaks} attribute) the score returned is the height score. For range \code{peaks} the product between heigh and width score is returned.
}
\value{
	In the case of \code{numeric} input, the value returned is a \code{data.frame} containing a 'peak' and a 'score' column. If the input is a \code{list}, the result will be a \code{list} of \code{data.frame}.

	If input is a \code{IRanges} or \code{IRangesList}, the result will be a RangedData object with one or multiple spaces respectively and a data column with the score.
}
\author{
	Oscar Flores \email{oflores@mmb.cpb.ub.es}
}
\seealso{
	\code{\link{peakDetection}}, \code{\link{processReads}}, 
}
\keyword{ manip }
\examples{

	#Generate a synthetic map
	map = syntheticNucMap(nuc.len=40, lin.len=130) #Trimmed length nucleosome map
	
	#Get the information of dyads and the coverage
	peaks = c(map$wp.starts, map$fz.starts)
	cover = filterFFT(coverage(map$syn.reads))

	#Calculate the scores
	scores = peakScoring(peaks, cover)
	plotPeaks(scores$peak, cover, scores=scores$score, start=5000, end=10000)
}
