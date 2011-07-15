\name{coverage.rpm}
\alias{coverage.rpm}
\title{
	Coverage calculation and normalization to reads per million (rpm)
}
\description{
	Calculates the coverage values from a \code{RangedData} object (or anything with a defined \code{coverage} function associated) and returns the coverage normalized to reads per million, allowing the comparison of experiments with a different absolut number of reads.
}
\usage{
coverage.rpm(data)
}
\arguments{
  \item{data}{
	\code{RangedData} (or compatible) with the reads information
}
}
\value{
	\code{RleList} object with the coverage objects
}
\author{
	Oscar Flores \email{oflores@mmb.pcb.ub.es}
}
\seealso{
	\code{\link{processReads}}, \code{\link{coverage}}
}
\keyword{manip}
\examples{
	
	#Load the example dataset and get the coverage
	data(nucleosome_htseq)
	cov = coverage.rpm(nucleosome_htseq)

	print(cov)

	#Plot it
	plot(as.vector(cov[["chr1"]]), type="l", ylab="coverage", xlab="position")
}
