

#' Correct experimental profiles with control sample
#' 
#' This function allows the correction of experimental coverage profiles
#' (usually MNase digested nucleosomal DNAs in this library) with control
#' samples (usually naked DNA sample digested with MNase). This is useful to
#' correct MNase biase.
#' 
#' This substracts the enrichment in the control sample respect it's mean from
#' the experimental profile.
#' 
#' This is useful for examinating the effect of the MNase digestion in
#' nucleosome experiments using a nucleosomal DNA and a genomic (naked) DNA
#' sample. Notice that genomic DNA samples cannot be strand-corrected using
#' single end data, so only paired end controls are useful for this proupose,
#' despite they can be compared against extended nucleosomal DNA single end
#' reads. Furthermore, both datasets must be converted to reads per milion.
#' 
#' This process dificults the nucleosome positioning due the lower sharpness of
#' the peaks, but allows a complementary study of the MNase digestion effect.
#' 
#' @aliases controlCorrection controlCorrection,SimpleRleList-method
#' controlCorrection,Rle-method controlCorrection,list-method
#' controlCorrection,numeric-method
#' @param exp,ctr Comparable experimental and control samples (this means same
#' format and equivalent preprocessment)
#' @param mc.cores Number of cores available for parallel list processing
#' @return Corrected experimental profile
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @keywords manip
#' @examples
#' 
#'     #Toy example
#'     map = syntheticNucMap(as.ratio=TRUE)
#' 
#'     exp = coverage(map$syn.reads)
#'     ctr = coverage(map$ctr.reads)
#' 
#'     corrected = controlCorrection(exp, ctr)
#' 
NULL





#' Coverage calculation and normalization to reads per million (rpm)
#' 
#' Calculates the coverage values from a \code{RangedData} object (or anything
#' with a defined \code{coverage} function associated) and returns the coverage
#' normalized to reads per million, allowing the comparison of experiments with
#' a different absolut number of reads.
#' 
#' 
#' @aliases coverage.rpm coverage.rpm,RangedData-method
#' coverage.rpm,GRanges-method
#' @param data \code{RangedData} (or compatible) with the reads information
#' @param scale By default, a million (1e6), but you could change this value
#' for abnormal high or low amount of reads.
#' @param \dots Additional arguments to be passed to \code{coverage} function
#' @return \code{RleList} object with the coverage objects
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{processReads}}, \code{\link{coverage}}
#' @keywords manip
#' @examples
#' 
#'     #Load the example dataset and get the coverage
#'     data(nucleosome_htseq)
#'     cov = coverage.rpm(nucleosome_htseq)
#' 
#'     print(cov)
#' 
#'     #Plot it
#'     plot(as.vector(cov[["chr1"]]), type="l", ylab="coverage", xlab="position")
#' 
NULL





#' Export ranges in BED format
#' 
#' Export ranges in BED format, compatible with UCSC genome browser, IGB, and
#' others.
#' 
#' 
#' @aliases export.bed .export.bed export.bed,IRanges-method
#' export.bed,CompressedIRangesList-method export.bed,RangedData-method
#' export.bed,GRanges-method
#' @param ranges Ranges to export, in \code{IRanges}, \code{IRangesList} or
#' \code{RangedData} format
#' @param score Score data if not included in \code{ranges} object. Bed file
#' will put all scores=1000 if scores are not present
#' @param chrom For single \code{IRanges} objects, the chromosome they
#' represent. For other data types, values from \code{names(...)} will be used.
#' @param name Name of the track
#' @param desc Description of the track
#' @param filepath Path and prefix of the file(s) to write. Chromosome number
#' and "bed" extension will be automatically added.
#' @param splitByChrom If multiple chromosomes are given, should they be
#' splitted into one file per chromosome or shall them be saved all together?
#' @return (none)
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @references BED format specification:
#' http://genome.ucsc.edu/FAQ/FAQformat#format1
#' @keywords file
#' @examples
#' 
#'     # Generate random ranges with scores
#'     ran <- RangedData(IRanges(start=1:100, end=101:200), score=(1:100) / 100)
#'     names(ran) <- "chrX"
#' 
#'     # Export as bed file
#'     export.bed(ran, name="test_track", desc="Just a test track")
#' 
#'     # If executed, this would create a file named "test_track.chrX.bed" with:
#' 
#'     #    track name="test_track" description="Just a test track" useScore=0
#'     # chrX    1 101 nucl1 0.01
#'     # chrX    2 102 nucl2 0.02
#'     # chrX    3 103 nucl3 0.03
#'     # ...
#' 
NULL





#' Clean noise and smoothing for genomic data using Fourier-analysis
#' 
#' Remove noise from genomic data smoothing and cleaning the observed signal.
#' This function doesn't alter the shape or the values of the signal as much as
#' the traditional method of sliding window average does, providing a great
#' correlation within the original and filtered data (>0.99).
#' 
#' Fourier-analysis principal components selection is widely used in signal
#' processing theory for an unbiased cleaning of a signal over the time.
#' 
#' Other procedures, as the traditional sliding window average, can change too
#' much the shape of the results in function of the size of the window, and
#' moreover they don't only smooth the noise without removing it.
#' 
#' With a Fourier Transform of the original signal, the input signal is
#' descomposed in diferent wavelets and described as a combination of them.
#' Long frequencies can be explained as a function of two ore more periodical
#' shorter frequecies. This is the reason why long, unperiodic sequences are
#' usually identified as noise, and therefore is desireable to remove them from
#' the signal we have to process.
#' 
#' This procedure here is applied to genomic data, providing a novel method to
#' obtain perfectly clean values wich allow an efficient detection of the peaks
#' which can be used for a direct nucleosome position recognition.
#' 
#' This function select a certain number of components in the original power
#' spectrum (the result of the Fast Fourier Transform which can be seen with
#' \code{showPowerSpec=TRUE}) and sets the rest of them to 0 (component
#' knock-out).
#' 
#' The amout of components to keep (given as a percentage of the input lenght)
#' can be set by the \code{pcKeepComp}. This will select the first components
#' of the signal, knock-outing the rest. If this value is close to 1, more
#' components will be selected and then more noise will be allowed in the
#' output. For an effective filtering which removes the noise keeping almost
#' all relevant peaks, a value between 0.01 and 0.05 is usually sufficient.
#' Lower values can cause merging of adjacent minor peaks.
#' 
#' This library also allows the automatic detection of a fitted value for
#' \code{pcKeepComp}. By default, if uses the \code{pcKeepCompDetect} function,
#' which looks which is the minimum percentage of components than can reproduce
#' the original signal with a corelation between the filtered and the original
#' one of 0.99. See the help page of \code{pcKeepCompDetect} for further
#' details and reference of available parameters.
#' 
#' One of the most powerful features of \code{nucleR} is the efficient
#' implementation of the FFT to genomic data. This is achived trought few
#' tweaks that allow an optimum performance of the Fourier Transform. This
#' includes a by-range filtering, an automatic detection of uncovered regions,
#' windowed execution of the filter and padding of the data till nearest power
#' of 2 (this ensures an optimum case for FFT due the high factorization of
#' components). Internal testing showed up that in specific datasets, these
#' optimizations lead to a dramatic improvement of many orders of magnitude
#' (from 3 days to few seconds) while keeping the correlation between the
#' native \code{fft} call and our \code{filterFFT} higher than 0.99. So, the
#' use of these optimizations is highly recomended.
#' 
#' If for some reason you want to apply the function without any kind of
#' optimizations you can specify the parameter \code{useOptim=FALSE} to bypass
#' them and get the pure knockout inverse from native FFT call. All other
#' parameters can be still applyied in this case.
#' 
#' @aliases filterFFT filterFFT,SimpleRleList-method filterFFT,Rle-method
#' filterFFT,list-method filterFFT,numeric-method
#' @param data Coverage or intensities values representing the results of the
#' NGS of TA experiment. This attribute could be a individual vector
#' representing a chromosome (\code{Rle} or \code{numeric} object) or a list of
#' them.
#' @param pcKeepComp Number of components to select, in percentage respect
#' total length of the sample. Allowed values are numeric (in range 0:1) for
#' manual setting or "auto" for automatic detection. See details.
#' @param showPowerSpec Plot the Power Spectrum of the Fast Fourier Transform
#' to visually identify the selected components (see details).
#' @param useOptim This function implements tweaks to a standard fft call to
#' improve (dramatically) the performance in large genomic data. These
#' optimizations can be bypassed by setting this parameter to \code{FALSE}.
#' @param mc.cores If multiple cores are available, maximum number of them to
#' use for parallel processing of \code{data} elements (only useful if
#' \code{data} is a list of elements)
#' @param \dots Other parameters to be passed to \code{pcKeepCompDetect}
#' function
#' @return Numeric vector with cleaned/smoothed values
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}, David Rosell
#' \email{david.rossell@@irbbarcelona.org}
#' @references Smith, Steven W. (1999), The Scientist and Engineer's Guide to
#' Digital Signal Processing (Second ed.), San Diego, Calif.: California
#' Technical Publishing, ISBN 0-9660176-3-3 (availabe online:
#' http://www.dspguide.com/pdfbook.htm)
#' @keywords manip
#' @examples
#' 
#'     #Load example data, raw hybridization values for Tiling Array
#'     raw_data = get(data(nucleosome_tiling))
#' 
#'     #Filter data
#'     fft_data = filterFFT(raw_data, pcKeepComp=0.01)
#' 
#'     #See both profiles
#'     par(mfrow=c(2,1), mar=c(3, 4, 1, 1))
#'     plot(raw_data, type="l", xlab="position", ylab="Raw intensities")
#'     plot(fft_data, type="l", xlab="position", ylab="Filtered intensities")
#' 
#'     #The power spectrum shows a visual representation of the components
#'     fft_data = filterFFT(raw_data, pcKeepComp=0.01, showPowerSpec=TRUE)
#' 
NULL





#' Fragments length detection from single-end sequencing samples
#' 
#' When using single-ended sequencing, the resulting partial sequences map only
#' in one strand, causing a bias in the coverage profile if not corrected. The
#' only way to correct this is knowing the average size of the real fragments.
#' \code{nucleR} uses this information when preprocessing single-ended
#' sequences. You can provide this information by your own (usually a 147bp
#' length is a good aproximation) or you can use this method to automatically
#' guess the size of the inserts.
#' 
#' This function shifts one strand downstream one base by one from
#' \code{min.shift} to \code{max.shift}. In every step, the correlation on a
#' random position of length \code{window} is checked between both strands. The
#' maximum correlation is returned and averaged for \code{samples} repetitions.
#' 
#' The final returned length is the best shift detected plus the width of the
#' reads. You can increase the performance of this function by reducing the
#' \code{samples} value and/or narrowing the shift range. The \code{window}
#' size has almost no impact on the performance, despite a to small value can
#' give biased results.
#' 
#' @aliases fragmentLenDetect fragmentLenDetect,AlignedRead-method
#' fragmentLenDetect,RangedData-method
#' @param reads Raw single-end reads (\code{AlignedRead} or \code{RangedData}
#' format)
#' @param samples Number of samples to perform the analysis (more = slower but
#' more accurate)
#' @param window Analysis window. Usually there's no need to touch this
#' parameter.
#' @param min.shift,max.shift Minimum and maximum shift to apply on the strands
#' to detect the optimal fragment size. If the range is too big, the
#' performance decreases.
#' @param as.shift If TRUE, returns the shift needed to align the middle of the
#' reads in opposite strand. If FALSE, returns the mean inferred fragment
#' length.
#' @param mc.cores If multicore support, maximum number of cores allowed to
#' use.
#' @return Inferred mean lenght of the inserts by default, or shift needed to
#' align strands if \code{as.shift=TRUE}
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @keywords attribute
#' @examples
#' 
#'     # Create a sinthetic dataset, simulating single-end reads, for positive and
#'     # negative strands
#'     # Positive strand reads
#'     pos = syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#'     # Negative strand (shifted 147bp)
#'     neg = IRanges(end=start(pos)+147, width=40)
#'     sim = RangedData(
#'         c(pos, neg),
#'         strand=c(rep("+", length(pos)), rep("-", length(neg)))
#'     )
#' 
#'     #Detect fragment lenght (we know by construction it is really 147)
#'     fragmentLenDetect(sim, samples=50)
#' 
#'     #The function restrict the sampling to speed up the example
#' 
NULL





#' Example reads from high-troughtput sequencing nucleosome positioning
#' experiment
#' 
#' Few reads from paired-ended MNase-seq experiment in S.cerevisiae where
#' mononucleosomes were sequenced
#' 
#' This data is obtained from MNase digested nucleosomal DNA and sequenced with
#' Illumina platform. Paired-ended reads where mapped to SacCer1 genome using
#' Bowtie, and imported to R using the package \code{ShortRead} and paired ends
#' where merged into a single range.
#' 
#' Reads were sorted by chromosome and starting position and only a few reads
#' from the starting positions of chromosome 1 are presented.
#' 
#' @name nucleosome_htseq
#' @docType data
#' @format \code{RangedData} with the range of the reads and a data column with
#' the strand information.
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
#' nucleosomal + 3 naked DNA) has been merged using package \code{Starr} and
#' the resulting \code{ExpressionSet} object has been passed to
#' \code{processTilingArray} function from this package as follows:
#' 
#' \code{processTilingArray(data, exprName, chrPAttern="Sc:Oct_2003;chr1",
#' closeGaps=50)}
#' 
#' The first 8000bp of the chr1 have been saved as this example dataset.
#' 
#' @name nucleosome_tiling
#' @docType data
#' @format \code{numeric} vector with the intensities.
#' @source Publication pending
#' @keywords datasets
NULL





#' Nucleosome positioning package for R
#' 
#' Nucleosome positioning from Tiling Arrays and High-Troughput Sequencing
#' Experiments
#' 
#' \tabular{ll}{ Package: \tab nucleR\cr Type: \tab Package\cr License: \tab
#' LGPL (>= 3)\cr LazyLoad: \tab yes\cr } This package provides a convenient
#' pipeline to process and analize nucleosome positioning experiments from
#' High-Troughtput Sequencing or Tiling Arrays.
#' 
#' Despite it's use is intended to nucleosome experiments, it can be also
#' useful for general ChIP experiments, such as ChIP-on-ChIP or ChIP-Seq.
#' 
#' See following example for a brief introduction to the available functions
#' 
#' @name nucleR-package
#' @aliases nucleR-package nucleR
#' @docType package
#' @author Oscar Flores Ricard Illa
#' 
#' Maintainer: Ricard Illa <ricard.illa@@irbbarcelona.org>
#' @keywords package
#' @examples
#' 
#'     #Load example dataset:
#'     # some NGS paired-end reads, mapped with Bowtie and processed with R
#'     # it is a RangedData object with the start/end coordinates for each read.
#'     reads = get(data(nucleosome_htseq))
#' 
#'     #Process the paired end reads, but discard those with length > 200
#'     preads_orig = processReads(reads, type="paired", fragmentLen=200)
#' 
#'     #Process the reads, but now trim each read to 40bp around the dyad
#'     preads_trim = processReads(reads, type="paired", fragmentLen=200, trim=40)
#' 
#'     #Calculate the coverage, directly in reads per million (r.p.m)
#'     cover_orig = coverage.rpm(preads_orig)
#'     cover_trim = coverage.rpm(preads_trim)
#' 
#'     #Compare both coverages, the dyad is much more clear in trimmed version
#'     t1 = as.vector(cover_orig[[1]])[1:2000]
#'     t2 = as.vector(cover_trim[[1]])[1:2000]
#'     t1 = (t1-min(t1))/max(t1-min(t1)) #Normalization
#'     t2 = (t2-min(t2))/max(t2-min(t2)) #Normalization
#'     plot(
#'         t1, type="l", lwd="2", col="blue", main="Original vs Trimmed coverage"
#'     )
#'     lines(t2, lwd="2", col="red")
#'     legend(
#'         "bottomright", c("Original coverage", "Trimmed coverage"), lwd=2,
#'         col=c("blue","red"), bty="n"
#'     )
#' 
#'     #Let's try to call nucleosomes from the trimmed version
#'     #First of all, let's remove some noise with FFT
#'     #Power spectrum will be plotted, look how with a 2%
#'     #of the components we capture almost all the signal
#'     cover_clean = filterFFT(cover_trim, pcKeepComp=0.02, showPowerSpec=TRUE)
#' 
#'     #How clean is now?
#'     plot(
#'         as.vector(cover_trim[[1]])[1:4000], t="l", lwd=2, col="red",
#'         main="Noisy vs Filtered coverage"
#'     )
#'     lines(cover_clean[[1]][1:4000], lwd=2, col="darkgreen")
#'     legend(
#'         "bottomright", c("Input coverage", "Filtered coverage"), lwd=2,
#'         col=c("red","darkgreen"), bty="n"
#'     )
#' 
#'     #And how similar? Let's see the correlation
#'     cor(cover_clean[[1]], as.vector(cover_trim[[1]]))
#' 
#'     #Now it's time to call for peaks, first just as points
#'     #See that the score is only a measure of the height of the peak
#'     peaks = peakDetection(cover_clean, threshold="25%", score=TRUE)
#'     plotPeaks(peaks[[1]], cover_clean[[1]], threshold="25%")
#' 
#'     #Do the same as previously, but now we will create the nucleosome calls:
#'     peaks = peakDetection(cover_clean, width=147, threshold="25%", score=TRUE)
#'     plotPeaks(peaks, cover_clean[[1]], threshold="25%")
#' 
#'     #This is all. From here, you can filter, merge or work with the nucleosome
#'     #calls using standard IRanges functions and R/Bioconductor manipulation
#' 
NULL





#' Detect peaks (local maximum) from values series
#' 
#' This function allows a efficient recognition of the local maximums (peaks)
#' in a given numeric vector.
#' 
#' It's recommended to smooth the input with \code{filterFFT} prior the
#' detection.
#' 
#' 
#' @aliases peakDetection peakDetection,list-method
#' peakDetection,numeric-method
#' @param data Input numeric values, or a list of them
#' @param threshold Threshold value from which the peaks will be selected. Can
#' be given as a percentage string (i.e., \code{"25%"} will use the value in
#' the 1st quantile of \code{data}) or as an absolute coverage numeric value
#' (i.e., \code{20} will not look for peaks in regions without less than 20
#' reads (or reads per milion)).
#' @param width If a positive integer > 1 is given, the peaks are returned as a
#' range of the given width centered in the local maximum. Useful for
#' nucleosome calling from a coverage peak in the dyad.
#' @param score If TRUE, the results will be scored using \code{peakScoring}
#' function.
#' @param min.cov Minimum coverage that a peak needs in order to be considered
#' as a nucleosome call.
#' @param mc.cores The number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. Parallelization requires at least two
#' cores.
#' @return The type of the return depends on the input parameters:
#' 
#' \code{numeric} (or a list of them) if \code{width==1 & score==FALSE}
#' containing the position of the peaks.
#' 
#' \code{data.frame} (or list of them) if \code{width==1 & score==TRUE}
#' containing a 'peak' column with the position of the peak plus a 'score'
#' column with its score.
#' 
#' \code{IRanges} (or \code{IRangesList}) if \code{width>1 & score==FALSE}
#' containing the ranges of the peaks.
#' 
#' \code{RangedData} if \code{width>1 & score==TRUE} containing the ranges of
#' the peaks and the assigned score.
#' @note If \code{width} > 1, those ranges outside the range
#' \code{1:length(data)} will be skipped.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{filterFFT}}, \code{\link{peakScoring}}
#' @keywords manip
#' @examples
#' 
#'     # Generate a random peaks profile
#'     reads = syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#'     cover = coverage(reads)
#' 
#'     # Filter them
#'     cover_fft = filterFFT(cover)
#' 
#'     # Detect and plot peaks (up a bit the threshold for accounting synthetic
#'     # data)
#'     peaks = peakDetection(cover_fft, threshold="40%", score=TRUE)
#'     plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#' 
#'     #Now use ranges version, which accounts for fuzziness when scoring
#'     peaks = peakDetection(cover_fft, threshold="40%", score=TRUE, width=147)
#'     plotPeaks(peaks, cover_fft, threshold="40%", start=10000, end=15000)
#' 
NULL





#' Peak scoring function
#' 
#' Scores peaks detected with function \code{peakDetection} according the
#' height and the sharpness (width) of the peak. This function can be called
#' automatically from \code{peakDetection} if \code{score=TRUE}.
#' 
#' This function scores each previously identified peak according its height
#' and sharpness.
#' 
#' The height score (\code{score_h}) tells how large is a peak, higher means
#' more coverage or intensity, so better positioned nucleosome. This score is
#' obtained by checking the observed peak value in a Normal distribution with
#' the mean and sd of \code{data}. This value is between 0 and 1.
#' 
#' The width score (\code{score_w}) is a mesure of how sharp is a peak. With a
#' NGS coverage in mind, a perfect phased (well-positioned) nucleosome is this
#' that starts and ends exactly in the same place many times. The shape of this
#' ideal peak will be a rectangular shape of the lenght of the read. A wider
#' top of a peak could indicate fuzzyness. The parameter \code{dyad.length}
#' tells how long should be the "flat" region of an ideal peak. The optimum
#' value for this parameter is the lenght of the read in single-ended data or
#' the \code{trim} value of the function \code{processReads}. For Tiling Array,
#' the default value should be fine.
#' 
#' This score is obtained calculating the ratio between the mean of the
#' nucleosome scope (the one provided by range in the elements of \code{peaks})
#' and the \code{dyad.length} central bases. This value is normalized between 0
#' and 1.
#' 
#' For punctual, single points peaks (provided by \code{numeric} vector or list
#' as \code{peaks} attribute) the score returned is the height score.
#' 
#' For range \code{peaks} the weighted sum of the heigth and width scores is
#' used. This is: \code{((score_h * weigth.height) / sum.wei) + ((score_w *
#' weigth.width) / sum.wei)}. Note that you can query for only one score by
#' weting its weight to 1 and the other to 0.
#' 
#' @aliases peakScoring peakScoring,list-method peakScoring,IRangesList-method
#' peakScoring,numeric-method peakScoring,IRanges-method
#' @param peaks The identified peaks resulting from \code{peakDetection}. Could
#' be a \code{numeric} vector with the position of the peaks, or a
#' \code{IRanges} object with the extended range of the peak. For both types,
#' list support is implemented as a \code{numeric} list or a \code{IRangesList}
#' @param data Data of nucleosome coverage or intensites.
#' @param threshold The non-default \code{threshold} previously used in
#' \code{peakDetection} function, if applicable. Can be given as a percentage
#' string (i.e., \code{"25%"} will use the value in the 1st quantile of
#' \code{data}) or as an absolute coverage numeric value (i.e., \code{20} will
#' not look for peaks in regions without less than 20 reads (or reads per
#' milion)).
#' @param dyad.length How many bases account in the nucleosome dyad for
#' sharpness description. If working with NGS data, works best with the reads
#' width value for single-ended data or the \code{trim} value given to the
#' \code{processReads} function.
#' @param weight.height,weight.width If the score is a range, the height and
#' the widht score (coverage and fuzzynes) can be defined with different
#' weigths with these parameters. See details.
#' @param mc.cores If input is a \code{list} or \code{IRangeList}, and multiple
#' cores support is available, the maximum number of cores for parallel
#' processing.
#' @return In the case of \code{numeric} input, the value returned is a
#' \code{data.frame} containing a 'peak' and a 'score' column. If the input is
#' a \code{list}, the result will be a \code{list} of \code{data.frame}.
#' 
#' If input is a \code{IRanges} or \code{IRangesList}, the result will be a
#' RangedData object with one or multiple spaces respectively and a 3 data
#' column with the mixed, width and heigh score.
#' @author Oscar Flores \email{oflores@@mmb.cpb.ub.es}
#' @seealso \code{\link{peakDetection}}, \code{\link{processReads}},
#' @keywords manip
#' @examples
#' 
#'     # Generate a synthetic map
#' 
#'     # Trimmed length nucleosome map
#'     map = syntheticNucMap(nuc.len=40, lin.len=130)
#' 
#'     # Get the information of dyads and the coverage
#'     peaks = c(map$wp.starts, map$fz.starts)
#'     cover = filterFFT(coverage(map$syn.reads))
#' 
#'     # Calculate the scores
#'     scores = peakScoring(peaks, cover)
#'     plotPeaks(scores$peak, cover, scores=scores$score, start=5000, end=10000)
#' 
NULL





#' Nucleosome calling plot function
#' 
#' Helper function for a quick and convenient overview of nucleosome calling
#' data.
#' 
#' This function is intended to plot data previously processed with
#' \code{nucleR} pipeline. It shows a coverage/intensity profile toghether with
#' the identified peaks. If available, score of each peak is also shown.
#' 
#' 
#' @aliases plotPeaks plotPeaks,numeric-method plotPeaks,data.frame-method
#' plotPeaks,RangedData-method plotPeaks,IRanges-method
#' @param peaks \code{numeric}, \code{data.frame}, \code{IRanges} or
#' \code{RangedData} object containing the detected peaks information. See help
#' of \code{\link{peakDetection}} or \code{\link{peakScoring}} for more
#' details.
#' @param data Coverage or Tiling Array intensities
#' @param threshold Threshold applied in \code{peakDetection}
#' @param scores If \code{peaks} is a \code{data.frame} or a \code{RangedData}
#' it's obtained from 'score' column, otherwise, \code{scores} can be given
#' here as a \code{numeric} vector.
#' @param start,end Start and end points defining a subset in the range of
#' \code{data}. This is a convenient way to plot only a small region of data,
#' without dealing with subsetting of range or score objects.
#' @param dyn.pos If peaks are ranges, should they be positioned dynamicaly on
#' top of the peaks or staticaly at \code{threshold} baseline. Spacing of
#' overlapping ranges is automatically applied if \code{FALSE}.
#' @param xlab,type,col.points Default values to be passed to \code{plot} and
#' \code{points}
#' @param thr.lty,thr.lwd,thr.col Default values to be passed to \code{abline}
#' for threshold representation
#' @param rect.thick,rect.lwd,rect.border Default values for \code{rect}
#' representation or ranges. \code{rect.thick} indicates the thickness (in
#' percentage relative to y-axis range) of the rectangles.
#' @param scor.col,scor.font,scor.adj,scor.cex,scor.digits Default values for
#' \code{text} representation for score numbers, if available.
#' @param indiv.scores Show or hide individual scores for width and height in
#' brakets besides the mixed score.
#' @param \dots Other parameters passed to \code{\link{plot}} function
#' @return (none)
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link{peakDetection}}, \code{\link{peakScoring}},
#' \code{\link{plot}},
#' @keywords hplot
#' @examples
#' 
#' 
#'     # Generate a random peaks profile
#'     reads = syntheticNucMap(nuc.len=40, lin.len=130)$syn.reads
#'     cover = coverage(reads)
#' 
#'     # Filter them
#'     cover_fft = filterFFT(cover)
#' 
#'     # Detect peaks
#'     peaks = peakDetection(cover_fft, threshold="40%", score=TRUE, width=140)
#' 
#'     # Plot peaks and coverage profile (show only a window)
#'     plotPeaks(peaks, cover_fft, threshold="40%", start=1000, end=6000)
#' 
NULL





#' Process reads from High-Troughtput Sequencing experiments
#' 
#' This method allows the processment of NGS nucleosome reads from different
#' sources and a basic manipulation of them. The tasks includes the correction
#' of strand-specific single-end reads and the trimming of reads to a given
#' length.
#' 
#' This function reads a \code{AlignedRead} or a \code{RangedData} object
#' containing the position, length and strand of the sequence reads.
#' 
#' It allows the processment of both paired and single ended reads. In the case
#' of single end reads this function corrects the strand-specific mapping by
#' shifting plus strand reads and minus strand reads towards a middle position
#' where both strands are overlaped. This is done by accounting the expected
#' fragment length (\code{fragmentLen}).
#' 
#' For paired end reads, mononucleosomal reads could extend more than expected
#' length due to mapping issues or experimental conditions. In this case, the
#' \code{fragmentLen} variable sets the threshold from which reads longer than
#' it should be ignored.
#' 
#' If no value is supplied for \code{fragmentLen} it will be calculated
#' automatically (increasing the computing time) using \code{fragmentLenDetect}
#' with default parameters. Performance can be increased by tunning
#' \code{fragmentLenDetect} parameteres in a separated call and passing its
#' result as \code{fragmentLen} parameter.
#' 
#' In some cases, could be useful trim the reads to a shorter length to improve
#' the detection of nucleosome dyads, easing its detection and automatic
#' positioning. The parameter \code{trim} allows the selection of how many
#' nucleotides select from each read.
#' 
#' A special case for single-ended data is setting the \code{trim} to the same
#' value as \code{fragmentLen}, so the reads will be extended strand-wise
#' towards the 3' direction, creating an artificial map comparable with
#' paired-ended data. The same but opposite can be performed with paired-end
#' data, setting a \code{trim} value equal to the read length from paired
#' ended, so paired-ended data will look like single-ended.
#' 
#' @aliases processReads processReads,AlignedRead-method
#' processReads,RangedData-method processReads,GRanges-method
#' @param data Sequence reads objects, probably imported using other packages
#' as \code{ShortRead}. Allowed object types are \code{AlignedRead} and
#' \code{RangedData} with a \code{strand} attribute.
#' @param type Describes the type of reads. Values allowed are \code{single}
#' for single-ended reads and \code{paired} for paired-ended.
#' @param fragmentLen Expected original length of the sequenced fragments. See
#' details.
#' @param trim Length to trim the reads (or extend them if \code{trim} > read
#' length)
#' @param \dots Other parameters passed to \code{fragmentLenDetect} if no fixed
#' \code{fragmentLen} is given.
#' @return \code{RangedData} containing the aligned/trimmed individual reads
#' @note \strong{IMPORTANT}: this information is only used to correct possible
#' strand-specific mapping, this package doesn't link the two ends of paired
#' reads.
#' @author Oscar Flores \email{oflores@@mmb.pcb.ub.es}
#' @seealso \code{\link[ShortRead]{AlignedRead}}, \code{\link{RangedData}},
#' \code{\link{fragmentLenDetect}}
#' @keywords manip
#' @examples
#' 
#'     # Load data
#'     data(nucleosome_htseq)
#' 
#'     # Process nucleosome reads, select only those shorter than 200bp
#'     pr1 = processReads(nucleosome_htseq, fragmentLen=200)
#' 
#'     # Now process them, but picking only the 40 bases surrounding the dyad
#'     pr2 = processReads(nucleosome_htseq, fragmentLen=200, trim=40)
#' 
#'     # Compare the results:
#'     par(mfrow=c(2,1), mar=c(3,4,1,1))
#'     plot(
#'         as.vector(coverage(pr1)[["chr1"]]), type="l",
#'         ylab="coverage (original)"
#'     )
#'     plot(
#'         as.vector(coverage(pr2)[["chr1"]]), type="l",
#'         ylab="coverage (trimmed)"
#'     )
#' 
NULL



