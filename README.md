# nucleR

[NucleR](http://bioconductor.org/packages/nucleR/) is an R/Bioconductor package for working with next generation sequencing and tilling arrays. It uses a novel aproach in this field which comprises a deep profile cleaning using Fourier Transform and peak scoring for a quick and flexible nucleosome calling.

The aim of this package is not providing an all-in-one data analysis pipeline but complement those existing specialized libraries for low-level data importation and pre-processment into R/Bioconductor framework.

[NucleR](http://bioconductor.org/packages/nucleR/) works with data from high-troughput technologies MNase-seq and ChIP-seq, and Tiling Microarrays (ChIP-on-Chip).

This is a brief summary of the main functions:

* Data import: `readBAM`, `processReads`, `processTilingArray`
* Data transformation: `coverage.rpm`, `filterFFT`, `controlCorrection`
* Nucleosome calling: `peakDetection`, `peakScoring`
* Visualization: `plotPeaks`
* Data generation: `syntheticNucMap`

This software was published in Bioinformatics Journal: Flores, O., and Orozco, M. (2011). nucleR: a package for non-parametric nucleosome positioning. Bioinformatics 27, 2149–2150.


## Installation
---------------
&nbsp;

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("nucleR", version = "3.8")


## Usage 
---------------


1- Load the package in R
&nbsp;

    library(nucleR)

2- Load the example data provided with the package, containing position of MNase-seq reads mapped to S. cerevisiae genome
&nbsp;

    data(nucleosome_htseq)
    class(nucleosome_htseq)
    nucleosome_htseq

3- Filter reads and remove noise:  discard the reads longer than 200bp (threshold given to only keep mononucleosomes), remove noise due to MNase efficiency by trimming reads to use only its central part (50bp around the dyad) 

&nbsp;
  
   reads_trim <- processReads(nucleosome_htseq, type="paired", fragmentLen=200, trim=50)

4- Obtain the normalized coverage (the count of how many reads are mapped to each position, divided by the total number of reads and multiplied by one milion)  

&nbsp;

    cover_trim <- coverage.rpm(reads_trim)

5- Smooth the coverage signal using the Fast Fourier Transformation

&nbsp;

    fft_ta <- filterFFT(cover_trim, pcKeepComp=0.01, showPowerSpec=TRUE)

6- Detect peaks in the smoothed coverage which correspond to nucleosome dyads and score them according to their fuzziness level

&nbsp;

   peaks <- peakDetection(fft_ta, threshold="25%", score=TRUE, width=147)


For more details about the functions and additional analyses refer to nucleR [manual](https://bioconductor.org/packages/release/bioc/manuals/nucleR/man/nucleR.pdf) and [vignette](https://bioconductor.org/packages/release/bioc/vignettes/nucleR/inst/doc/nucleR.pdf).



