#' Example reads from high-troughtput sequencing nucleosome positioning
#' experiment
#'
#' Few reads from paired-ended MNase-seq experiment in S.cerevisiae where
#' mononucleosomes were sequenced
#'
#' This data is obtained from MNase digested nucleosomal DNA and sequenced with
#' Illumina platform. Paired-ended reads where mapped to SacCer1 genome using
#' Bowtie, and imported to R using the package `ShortRead` and paired ends
#' where merged into a single range.
#'
#' Reads were sorted by chromosome and starting position and only a few reads
#' from the starting positions of chromosome 1 are presented.
#'
#' @name nucleosome_htseq
#' @docType data
#' @format `RangedData` with the range of the reads and a data column with the
#'   strand information.
#' @source Publication pending
#' @keywords datasets
NULL

#' Example intensities from Tiling Microarray nucleosome positioning experiment
#'
#' Some bases from S.cerevisiae tiling microarray where mononucleosomes were
#' sequenced and hybridizated with histone-free naked DNA. The intensity is the
#' normalized ratio between the intensities from nucleosomic and naked DNA.
#'
#' Due to the difficulty of providing a raw file, this file has been
#' preprocessed. See details.
#'
#' The raw .CEL files from Affymetrix S.Cerevisiae Tilling 1.0R Array (3
#' nucleosomal + 3 naked DNA) has been merged using package `Starr` and the
#' resulting `ExpressionSet` object has been passed to `processTilingArray`
#' function from this package as follows:
#'
#' `processTilingArray(data, exprName, chrPAttern="Sc:Oct_2003;chr1",
#' closeGaps=50)`
#'
#' The first 8000bp of the chr1 have been saved as this example dataset.
#'
#' @name nucleosome_tiling
#' @docType data
#' @format `numeric` vector with the intensities.
#' @source Publication pending
#' @keywords datasets
NULL
