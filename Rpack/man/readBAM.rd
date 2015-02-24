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
  readBAM(file, type = "single", mc.cores = 1)
}
\arguments{
    \item{file}{
        List of input BAM files.
    }
    \item{type}{
        Describes the type of reads. Values allowed are \code{single} for
        single-ended reads and \code{paired} for pair-ended.
    }
    \item{mc.cores}{
        If \code{multicore} support, the number of cores available.
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
