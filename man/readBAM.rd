\name{readBAM}
\alias{readBAM}
\title{
    Import reads from a list of BAM files.
}
\description{
    This function allows to load reads from BAM files from both single and
    paired-end commming from Next Generation Sequencing nucleosome mapping
    experiments.
}
\usage{
  readBAM(file, type = "paired")
}
\arguments{
    \item{file}{
        List of input BAM files.
    }
    \item{type}{
        Describes the type of reads. Values allowed are \code{single} for
        single-ended reads and \code{paired} for pair-ended.
    }
}
\value{
    List of \code{GRanges} containing the reads of each input BAM file.
}
\author{
    Oscar Flores \email{oflores@mmb.pcb.ub.es},
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{
    file
}
