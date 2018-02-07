\name{plotPeaks}
\alias{plotPeaks}
\alias{plotPeaks,numeric-method}
\alias{plotPeaks,data.frame-method}
\alias{plotPeaks,RangedData-method}
\alias{plotPeaks,IRanges-method}
\title{
Nucleosome calling plot function
}
\description{
Helper function for a quick and convenient overview of nucleosome calling data.

This function is intended to plot data previously processed with \code{nucleR}
pipeline. It shows a coverage/intensity profile toghether with the identified
peaks. If available, score of each peak is also shown.
}
\usage{
plotPeaks(peaks, data, ...)
\S4method{plotPeaks}{IRanges}(
    peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
    dyn.pos=TRUE, xlab="position", type="l", col.points="red", thr.lty=1,
    thr.lwd="1", thr.col="darkred", rect.thick=2, rect.lwd=1,
    rect.border="black", scor.col=col.points, scor.font=2, scor.adj=c(0.5,0),
    scor.cex=0.75, scor.digits=2, indiv.scores=TRUE, ...
)
\S4method{plotPeaks}{numeric}(
    peaks, data, threshold=0, scores=NULL, start=1, end=length(data),
    xlab="position", type="l", col.points="red", thr.lty=1,  thr.lwd="1",
    thr.col="darkred", scor.col=col.points, scor.font=2, scor.adj=c(0.5,0),
    scor.cex=0.75, scor.digits=2,...
)
}
\arguments{
\item{peaks}{
\code{numeric}, \code{data.frame}, \code{IRanges} or \code{RangedData} object
containing the detected peaks information. See help of
\code{\link{peakDetection}} or \code{\link{peakScoring}} for more details.
}
\item{data}{
Coverage or Tiling Array intensities
}
\item{threshold}{
Threshold applied in \code{peakDetection}
}
\item{scores}{
If \code{peaks} is a \code{data.frame} or a \code{RangedData} it's obtained
from 'score' column, otherwise, \code{scores} can be given here as a
\code{numeric} vector.
}
\item{start, end}{
Start and end points defining a subset in the range of \code{data}. This is a
convenient way to plot only a small region of data, without dealing with
subsetting of range or score objects.
}
\item{dyn.pos}{
If peaks are ranges, should they be positioned dynamicaly on top of the peaks
or staticaly at \code{threshold} baseline. Spacing of overlapping ranges is
automatically applied if \code{FALSE}.
}
\item{xlab, type, col.points}{
Default values to be passed to \code{plot} and \code{points}
}
\item{thr.lty, thr.lwd, thr.col}{
Default values to be passed to \code{abline} for threshold representation
}
\item{rect.thick, rect.lwd, rect.border}{
Default values for \code{rect} representation or ranges. \code{rect.thick}
indicates the thickness (in percentage relative to y-axis range) of the
rectangles.
}
\item{scor.col, scor.font, scor.adj, scor.cex, scor.digits}{
Default values for \code{text} representation for score numbers, if available.
}
\item{indiv.scores}{
Show or hide individual scores for width and height in brakets besides the
mixed score.
}
\item{\dots}{
Other parameters passed to \code{\link{plot}} function
}
}
\value{
(none)
}
\author{
Oscar Flores \email{oflores@mmb.pcb.ub.es}
}
\seealso{
\code{\link{peakDetection}}, \code{\link{peakScoring}}, \code{\link{plot}},
}
\keyword{hplot}
\examples{

    # Generate a random peaks profile
    reads = syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
    cover = coverage(reads)

    # Filter them
    cover_fft = filterFFT(cover)

    # Detect peaks
    peaks = peakDetection(cover_fft, threshold="40\%", score=TRUE, width=140)

    # Plot peaks and coverage profile (show only a window)
    plotPeaks(peaks, cover_fft, threshold="40\%", start=1000, end=6000)
}

