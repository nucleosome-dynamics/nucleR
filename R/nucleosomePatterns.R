#' Nucleosome patterns arround transcription start sites
#'
#' The transcription start site for each gene can be classified according to
#' the properties of the nucleosomes surrounding it (minus 1 and plus 1).
#' Those nucleosomes can be either _well-positioned_ (`W`), _fuzzy_ (`W`) or
#' _missing_. And, when present, they can be at a distance less than a given
#' threshold (by default `215` bp.), which creates an _close conformation_ or
#' at a distance bigger than the threshold, which makes an _open conformation_.
#'
#' @param genes a `data.frame`, a [GenomicRanges::GRanges] or a
#'   [GenomicFeatures::TxDb] object representing the genes of the used organism
#' @param calls nucleosome calls performed with `nucleR`
#' @param window window in bp. around the start of the gene where the plus1 and
#'   minus1 are searched. If no nucleosome is found in this window, it's
#'   reported as `missing`.
#' @param p1.max.downstream the nucleosome plus1 can be upstream of the
#'   transcription start site, but only this amount of base pairs.
#' @param open.thresh open/closed conformation threshold. Distances smaller
#'   than this will be classified as _close conformation_, while distances
#'   bigger than this will be classified as _open conformation_.
#' @param position whether our we are computing transcription start sites
#'   (`tss`) or transcription termination sites (`tts`).
#' @param col.id,col.pos,col.strand,col.chrom names of the columns in the genes
#'   `data.frame` containing the information corresponding to the gene id, the
#'   position to compute (usually the transcription start site), the strand and
#'   the chromosome for each gene.
#'
#' @return a `data.frame` with the following columns:
#'   * seqnames Name of the chromosome of the gene
#'   * id Name ID of the gene
#'   * p1.pos position of the plus 1 nucleosome
#'   * m1.pos position of the minus 1 nucleosome
#'   * dist distance in base pairs between `p1.pos` and `m1.pos`
#'   * descr descriptor of the conformation (W-open-W, W-close-F, +1_missing,
#'     -1_missing, etc).
#'
#' @author Ricard Illa \email{ricard.illa@irbbarcelona.org}
#' @keywords mainp
#' @rdname nucleosomePatterns
#'
#' @export nucleosomePatterns
#'
setGeneric(
    "nucleosomePatterns",
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id=NULL, col.chrom=NULL, col.strand=NULL,
              col.pos=NULL)
        standardGeneric("nucleosomePatterns")
)

#' @rdname nucleosomePatterns
setMethod(
    "nucleosomePatterns",
    signature(genes="data.frame", calls="data.frame"),
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id="name", col.chrom="chrom",
              col.strand="strand", col.pos="start")
        .nucleosomePatternsDF(genes = genes,
                              calls = calls,

                              window            = window,
                              p1.max.downstream = p1.max.downstream,
                              open.thresh       = open.thresh,
                              position          = position,

                              col.id     = col.id,
                              col.chrom  = col.chrom,
                              col.strand = col.strand,
                              col.pos    = col.pos)
)

#' @rdname nucleosomePatterns
#' @importFrom dplyr bind_rows
setMethod(
    "nucleosomePatterns",
    signature(genes="data.frame", calls="list"),
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id="name", col.chrom="chrom",
              col.strand="strand", col.pos="start")
        .nucleosomePatternsDF(genes = genes,
                              calls = bind_rows(calls),

                              window            = window,
                              p1.max.downstream = p1.max.downstream,
                              open.thresh       = open.thresh,
                              position          = position,

                              col.id     = col.id,
                              col.chrom  = col.chrom,
                              col.strand = col.strand,
                              col.pos    = col.pos)
)

#' @rdname nucleosomePatterns
#' @importFrom dplyr rename mutate
#' @importFrom magrittr %>%
setMethod(
    "nucleosomePatterns",
    signature(genes="data.frame", calls="GRanges"),
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id="name", col.chrom="chrom",
              col.strand="strand", col.pos="start")
    {
        calls %>%
            data.frame %>%
            rename(chromosome=seqnames) %>%
            mutate(peak=.mid(.)) -> df

        nucleosomePatterns(genes = genes,
                           calls = df,

                           window            = window,
                           p1.max.downstream = p1.max.downstream,
                           open.thresh       = open.thresh,
                           position          = position,

                           col.id     = col.id,
                           col.chrom  = col.chrom,
                           col.strand = col.strand,
                           col.pos    = col.pos)
    }
)
globalVariables(c(".", "seqnames"))

#' @rdname nucleosomePatterns
setMethod(
    "nucleosomePatterns",
    signature(genes="GRanges"),
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id="name", col.chrom="seqnames",
              col.strand="strand", col.pos=NULL)
    {
        df <- data.frame(genes)
        if (is.null(col.pos)) {
            df$pos <- ifelse(df$strand == "+", df$start,
                      ifelse(df$strand == "-", df$end,
                                               NA))
            col.pos <- "pos"
        }
        # call the data.frame method
        nucleosomePatterns(genes = df,
                           calls = calls,

                           window            = window,
                           p1.max.downstream = p1.max.downstream,
                           open.thresh       = open.thresh,
                           position          = position,

                           col.id     = col.id,
                           col.chrom  = col.chrom,
                           col.strand = col.strand,
                           col.pos    = col.pos)
    }
)

#' @rdname nucleosomePatterns
#' @importClassesFrom GenomicFeatures TxDb
#' @importMethodsFrom GenomicFeatures transcripts
setMethod(
    "nucleosomePatterns",
    signature(genes="TxDb"),
    function (genes, calls, window=3000, p1.max.downstream=20, open.thresh=215,
              position="tss", col.id="tx_name", col.chrom="seqnames",
              col.strand="strand", col.pos=NULL)
        # Call the GRanges method
        nucleosomePatterns(genes = transcripts(genes),
                           calls = calls,

                           window            = window,
                           p1.max.downstream = p1.max.downstream,
                           open.thresh       = open.thresh,
                           position          = position,

                           col.id     = col.id,
                           col.chrom  = col.chrom,
                           col.strand = col.strand,
                           col.pos    = col.pos)
)

#' @importFrom dplyr group_by group_by_ arrange arrange_ filter mutate
#' @importFrom magrittr %>%
#' @importFrom purrr pmap
#' @importFrom tidyr nest unnest
.nucleosomePatternsDF <- function (genes, calls, window=300,
                                   p1.max.downstream=20, open.thresh=215,
                                   position="tss", col.id="name",
                                   col.pos="start", col.strand="strand",
                                   col.chrom="chrom")
{
    genes %>%
        group_by_(col.chrom) %>%
        nest %>%
        arrange_(col.chrom) -> genes.t

    calls %>%
        data.frame %>%
        group_by(chromosome) %>%
        nest %>%
        arrange(chromosome) -> calls.t

    chrs <- intersect(genes.t[[col.chrom]], calls.t$chromosome)
    if (length(chrs) == 0) {
        warning("no chromosome name in common between the genes and the calls")
    }
    genes.t %>% filter(.[[col.chrom]] %in% chrs) -> genes.t
    calls.t %>% filter(.$chromosome %in% chrs) -> calls.t

    genes.t %>%
        mutate(data=pmap(list(.$data, calls.t$data, .[[col.chrom]]),
                         .nucleosomePatternsChr,
                         window        = window,
                         open.thresh   = open.thresh,
                         position      = position,
                         col.id        = col.id,
                         col.pos       = col.pos,
                         col.strand    = col.strand)) %>%
        unnest %>%
        data.frame
}
globalVariables(c(".", "chromosome"))

#' @importFrom dplyr rowwise do
#' @importFrom magrittr %>%
.nucleosomePatternsChr <- function (genes, calls, chrom, window=300,
                                    p1.max.downstream=20, open.thresh=215,
                                    position="tss", col.id="name",
                                    col.pos="start", col.strand="strand")
{
    message("computing chromosome ", chrom)
    genes %>%
        rowwise %>%
        do(.getPattern(calls, .[[col.id]], .[[col.pos]], .[[col.strand]],
                       window            = window,
                       p1.max.downstream = p1.max.downstream,
                       open.thresh       = open.thresh,
                       position          = position)) %>%
        data.frame
}
globalVariables(".")

.getPattern <- function (calls, id, pos, strand="+", window=300,
                         p1.max.downstream=20, open.thresh=215, position="tss")
{   # The start and end will be the closest nucleosome found nearby the TSS
    # those positions will also be the p1 and m1 positions if they are within
    # a -/+ window
    p1 <- .detectNuc(calls, pos, p1.max.downstream, strand, "p1", position)

    no.p1s <- nrow(p1) == 0  # no nucleosome found upstream of the TSS
    if (no.p1s) {
        p1.pos <- pos
    } else {
        p1.pos <- p1$peak
    }

    if (!no.p1s && .checkPos(p1.pos, pos, window)) {
        # there's a nucleosome upstream and within the window
        p1.nuc <- p1$peak
        p1.class <- p1$class

        # look for m1 relative to p1
        m1 <- .detectNuc(calls, p1.nuc, 0, strand, "m1", position)

        no.m1s <- nrow(m1) == 0
        if (no.m1s) {
            m1.pos <- pos
        } else {
            m1.pos <- m1$peak
        }

        if (!no.m1s && .checkPos(m1.pos, pos, window)) {
            # and there's also one downstream
            m1.class <- m1$class
            m1.nuc <- m1$peak
        } else {
            m1.nuc <- NA
            m1.class <- NULL
        }
    } else {
        p1.nuc <- NA
        m1.nuc <- NA
        p1.class <- NULL
        m1.class <- NULL

        # no p1... so look for m1 relative to the TSS
        m1.pos <- .mid(.detectNuc(calls, pos, 0, strand, "m1", position))
        if (length(m1.pos) == 0) {
            m1.pos <- pos
        }
    }

    dist <- abs(p1.nuc - m1.nuc)
    descr <- .getDescr(dist, m1.class, p1.class, open.thresh)

    data.frame(id     = id,
               p1.pos = p1.nuc,
               m1.pos = m1.nuc,
               dist   = dist,
               descr  = descr,
               stringsAsFactors=FALSE)
}

.detectNuc <- function (calls, pos, margin, strand, nucpos, position="tss")
{
    a <- strand == "+" && nucpos == "p1"
    b <- strand == "-" && nucpos == "m1"
    c <- strand == "-" && nucpos == "p1"
    d <- strand == "+" && nucpos == "m1"

    flipper <- ifelse(position == "tss", `(`, `!`)

    after <- flipper(a || b)
    before <- flipper(c || d)

    if (after) {
        shift <- `-`
        comp <- `>`
        closest <- which.min
    } else if (before) {
        shift <- `+`
        comp <- `<`
        closest <- which.max
    }

    subxs <- calls[comp(calls$peak, shift(pos, margin)), ]
    return(subxs[closest(subxs$peak), ])
}

.gimmeDist <- function(x, open.thresh)
    ifelse(is.na(x),        "-",
    ifelse(x > open.thresh, "open",
    ifelse(x < 120,         "overlap",
                            "close")))

.getDescr <- function (dist, m1.class, p1.class, open.thresh)
{
    if (!is.null(m1.class) & !is.null(p1.class)) {
        dist.class <- .gimmeDist(dist, open.thresh)
        return(paste(m1.class, dist.class, p1.class, sep="-"))
    } else if (is.null(p1.class)) {
        return("+1_missing")
    } else if (is.null(m1.class)) {
        return("-1_missing")
    } else {
        return(NA)
    }
}

.checkPos <- function (nuc.pos, pos, window)
    nuc.pos > (pos-window) && nuc.pos < (pos+window)
