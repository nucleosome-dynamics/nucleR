context("calculating coverage in reads per million")

test_that("test coverage.rpm on GRanges", {
    data("nucleosome_htseq")
    cov.sum <- sum(coverage.rpm(nucleosome_htseq))
    expect_equal(cov.sum, c(chr1=175075607))
})

test_that("test coverage.rpm on RangedData", {
    data("nucleosome_htseq")
    rd <- RangedData(ranges=ranges(nucleosome_htseq),
                     space=seqnames(nucleosome_htseq))
    cov.sum <- sum(coverage.rpm(rd))
    expect_equal(cov.sum, c(chr1=175075607))
})

test_that("test coverage.rpm on IRanges", {
    data("nucleosome_htseq")
    ir <- ranges(nucleosome_htseq)
    cov.sum <- sum(coverage.rpm(ir))
    expect_equal(cov.sum, 175075607)
})
