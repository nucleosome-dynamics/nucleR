\name{nucleosome_tiling}
\alias{nucleosome_tiling}
\docType{data}
\title{
	Example intensities from Tiling Microarray nucleosome positioning experiment
}
\description{
	Some bases from S.cerevisiae tiling microarray where mononucleosomes were sequenced and hybridizated
with histone-free naked DNA. The intensity is the normalized ratio between the intensities from nucleosomic and naked DNA.

	Due to the difficulty of providing a raw file, this file has been preprocessed. See details.
}
\usage{data(nucleosome_tiling)}
\format{
	\code{numeric} vector with the intensities.
}
\details{
	The raw .CEL files from Affymetrix S.Cerevisiae Tilling 1.0R Array (3 nucleosomal + 3 naked DNA) has been merged using package \code{Starr} and the resulting \code{ExpressionSet} object has been passed to \code{processTilingArray} function from this package as follows:

	\code{processTilingArray(data, exprName, chrPAttern="Sc:Oct_2003;chr1", closeGaps=50)}

	The first 8000bp of the chr1 have been saved as this example dataset.
}
\source{
	Publication pending
}
\keyword{datasets}
