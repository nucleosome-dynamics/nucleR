### R code from vignette source 'nucleR.Rnw'

###################################################
### code chunk number 1: loadTA
###################################################
require(IRanges)
library(nucleR)
data(nucleosome_tiling)
head(nucleosome_tiling, n=25)


###################################################
### code chunk number 2: import
###################################################
data(nucleosome_htseq)
class(nucleosome_htseq)
nucleosome_htseq


###################################################
### code chunk number 3: processReads
###################################################
# Process the paired end reads, but discard those with length > 200
reads_pair <- processReads(nucleosome_htseq, type="paired", fragmentLen=200)

# Process the reads, but now trim each read to 40bp around the dyad
reads_trim <- processReads(nucleosome_htseq, type="paired",
                           fragmentLen=200, trim=40)


###################################################
### code chunk number 4: coverage
###################################################
# Calculate the coverage, directly in reads per million (r.p.m)
cover_pair <- coverage.rpm(reads_pair)
cover_trim <- coverage.rpm(reads_trim)


###################################################
### code chunk number 5: figcover
###################################################
# Compare both coverages
t1 <- as.vector(cover_pair[[1]])[1:2000]
t2 <- as.vector(cover_trim[[1]])[1:2000]
t1 <- (t1 - min(t1)) / max(t1 - min(t1)) # Normalization
t2 <- (t2 - min(t2)) / max(t2 - min(t2)) # Normalization
par(mar=c(5, 4, 1, 1), cex=0.7)
plot(t1, type="l", lwd="1", col="blue", ylab="norm. coverage", xlab="position")
lines(t2, lwd="1", col="orange")


###################################################
### code chunk number 6: figcover
###################################################
# Compare both coverages
t1 <- as.vector(cover_pair[[1]])[1:2000]
t2 <- as.vector(cover_trim[[1]])[1:2000]
t1 <- (t1 - min(t1)) / max(t1 - min(t1)) # Normalization
t2 <- (t2 - min(t2)) / max(t2 - min(t2)) # Normalization
par(mar=c(5, 4, 1, 1), cex=0.7)
plot(t1, type="l", lwd="1", col="blue", ylab="norm. coverage", xlab="position")
lines(t2, lwd="1", col="orange")


###################################################
### code chunk number 7: mnase
###################################################
# Toy example
map <- syntheticNucMap(as.ratio=TRUE, wp.num=50, fuz.num=25)

exp <- coverage(map$syn.reads)
ctr <- coverage(map$ctr.reads)

corrected <- controlCorrection(exp, ctr)


###################################################
### code chunk number 8: figmnase
###################################################
# Toy example
par(mar=c(5,4,1,1), cex=0.7)
plot(exp, col="darkgrey", type="l", ylab="coverage", xlab="position")
lines(corrected, col="darkblue", lty=2)
legend("top", c("normal", "corrected"), hor=TRUE, bty="n",
       fill=c("darkgrey", "blue"))


###################################################
### code chunk number 9: fignoise
###################################################
par(mar=c(5,4,1,1), mfcol=c(1,4), cex=0.4)
plot(nucleosome_tiling[1:2000], type="l", main="original", ylim=c(-3,2),
     xlab="position", ylab="intensity")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,20))/20), type="l",
     main="sliding w. 20bp", ylim=c(-3,2), xlab="position", ylab="")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,50))/50), type="l",
     main="sliding w. 50bp", ylim=c(-3,2), xlab="position", ylab="")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,100))/100), type="l",
     main="sliding w. 100bp", ylim=c(-3,2), xlab="position", ylab="")


###################################################
### code chunk number 10: fignoise
###################################################
par(mar=c(5,4,1,1), mfcol=c(1,4), cex=0.4)
plot(nucleosome_tiling[1:2000], type="l", main="original", ylim=c(-3,2),
     xlab="position", ylab="intensity")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,20))/20), type="l",
     main="sliding w. 20bp", ylim=c(-3,2), xlab="position", ylab="")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,50))/50), type="l",
     main="sliding w. 50bp", ylim=c(-3,2), xlab="position", ylab="")
plot(filter(nucleosome_tiling[1:2000], c(rep(1,100))/100), type="l",
     main="sliding w. 100bp", ylim=c(-3,2), xlab="position", ylab="")


###################################################
### code chunk number 11: fft
###################################################
fft_ta = filterFFT(nucleosome_tiling, pcKeepComp=0.01, showPowerSpec=TRUE)


###################################################
### code chunk number 12: figfft
###################################################
par(mar=c(6,4,1,1), cex=0.7)
fft_ta <- filterFFT(nucleosome_tiling, pcKeepComp=0.01, showPowerSpec=TRUE)


###################################################
### code chunk number 13: figfft
###################################################
par(mar=c(6,4,1,1), cex=0.7)
fft_ta <- filterFFT(nucleosome_tiling, pcKeepComp=0.01, showPowerSpec=TRUE)


###################################################
### code chunk number 14: figfft3
###################################################
par(mar=c(5, 4, 1, 1), mfcol=c(2, 1), cex=0.7)
tiling_raw <- nucleosome_tiling[1:3000]
tiling_fft <- filterFFT(tiling_raw, pcKeepComp=0.01)
htseq_raw <- as.vector(cover_trim[[1]])[1:3000]
htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.02)
plot(tiling_raw, type="l", col="darkblue", lwd=3, ylab="intensity", xlab="")
lines(tiling_fft, type="l", col="cyan")
plot(htseq_raw, type="l", col="darkred", lwd=3, ylab="coverage",
     xlab="position")
lines(htseq_fft, type="l", col="pink")


###################################################
### code chunk number 15: figfft3
###################################################
par(mar=c(5, 4, 1, 1), mfcol=c(2, 1), cex=0.7)
tiling_raw <- nucleosome_tiling[1:3000]
tiling_fft <- filterFFT(tiling_raw, pcKeepComp=0.01)
htseq_raw <- as.vector(cover_trim[[1]])[1:3000]
htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.02)
plot(tiling_raw, type="l", col="darkblue", lwd=3, ylab="intensity", xlab="")
lines(tiling_fft, type="l", col="cyan")
plot(htseq_raw, type="l", col="darkred", lwd=3, ylab="coverage",
     xlab="position")
lines(htseq_fft, type="l", col="pink")


###################################################
### code chunk number 16: corfft
###################################################
tiling_raw <- nucleosome_tiling
tiling_fft <- filterFFT(tiling_raw, pcKeepComp=0.01)
htseq_raw <- as.vector(cover_trim[[1]])
htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.02)

cor(tiling_raw, tiling_fft, use="complete.obs")
cor(htseq_raw, htseq_fft, use="complete.obs")


###################################################
### code chunk number 17: peaks
###################################################
peaks <- peakDetection(htseq_fft, threshold="25%", score=FALSE)
peaks


###################################################
### code chunk number 18: plotpeak
###################################################
plotPeaks(peaks, htseq_fft, threshold="25%")


###################################################
### code chunk number 19: figpeak
###################################################
par(cex=0.7, mar=c(5,4,1,1))
plotPeaks(peaks, htseq_fft, threshold="25%", yaxt="n", ylab="coverage")


###################################################
### code chunk number 20: peaks2
###################################################
peaks <- peakDetection(htseq_fft, threshold="25%", score=TRUE)
head(peaks)


###################################################
### code chunk number 21: figpeak2
###################################################
par(cex=0.7, mar=c(5,4,1,1))
plotPeaks(peaks, htseq_fft, threshold="25%", yaxt="n", xlim=c(1,4000),
          ylab="coverage")


###################################################
### code chunk number 22: peaks3
###################################################
peaks <- peakDetection(htseq_fft, threshold="25%", score=TRUE, width=140)
peaks


###################################################
### code chunk number 23: figpeak3
###################################################
par(cex=0.7, mar=c(5,4,1,1))
plotPeaks(peaks, htseq_fft, threshold="25%", yaxt="n", xlim=c(1,4000),
          rect.lwd=0.5, rect.thick=2.5, ylab="coverage")


###################################################
### code chunk number 24: ranges
###################################################
nuc_calls <- ranges(peaks[peaks$score > 0.1,])[[1]]
red_calls <- reduce(nuc_calls)
red_class <- RangedData(red_calls, isFuzzy=width(red_calls) > 140)
red_class


###################################################
### code chunk number 25: figranges
###################################################
par(cex=0.7, mar=c(5,4,1,1))
plotPeaks(red_calls, htseq_fft, threshold="25%", yaxt="n", xlim=c(1,4000), rect.lwd=0.5, rect.thick=3, ylab="coverage")


###################################################
### code chunk number 26: syn
###################################################
syntheticNucMap(wp.num=100, wp.del=10, wp.var=30, fuz.num=20, fuz.var=50,
                max.cover=20, nuc.len=147, lin.len=20, rnd.seed=1,
                as.ratio=TRUE, show.plot=TRUE)


###################################################
### code chunk number 27: figsyn
###################################################
par(mar=c(5,4,1,1), cex=0.7)
syn <- syntheticNucMap(wp.num=100, wp.del=10, wp.var=30, fuz.num=20,
                       fuz.var=50, max.cover=20, nuc.len=147, lin.len=20,
                       rnd.seed=1, as.ratio=TRUE, show.plot=TRUE, ylab="",
                       xlab="position")


###################################################
### code chunk number 28: figsyn
###################################################
par(mar=c(5,4,1,1), cex=0.7)
syn <- syntheticNucMap(wp.num=100, wp.del=10, wp.var=30, fuz.num=20,
                       fuz.var=50, max.cover=20, nuc.len=147, lin.len=20,
                       rnd.seed=1, as.ratio=TRUE, show.plot=TRUE, ylab="",
                       xlab="position")


