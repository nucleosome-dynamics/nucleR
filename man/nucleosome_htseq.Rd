\name{nucleosome_htseq}
\alias{nucleosome_htseq}
\docType{data}
\title{
	Example reads from high-troughtput sequencing nucleosome positioning experiment
}
\description{
	Few reads from paired-ended MNase-seq experiment in S.cerevisiae where mononucleosomes were sequenced
}
\usage{data(nucleosome_htseq)}
\format{
	\code{RangedData} with the range of the reads and a data column with the strand information.
}
\details{
	This data is obtained from MNase digested nucleosomal DNA and sequenced with Illumina platform. Paired-ended reads where mapped to SacCer1 genome using Bowtie, and imported to R using the package \code{ShortRead} and paired ends where merged into a single range.

	Reads were sorted by chromosome and starting position and only a few reads from the starting positions of chromosome 1 are presented.
}
\source{
	Publication pending
}
\keyword{datasets}
