
test_that("cnaCalling() returns a correct output", {

    data("muscadet_obj", package = "muscadet")

    muscadet_cna <- cnaCalling(muscadet_obj,
                               omics.coverage = "ATAC",
                               depthmin.a.clusters = 3, # set low thresholds for example data
                               depthmin.c.clusters = 5,
                               depthmin.a.allcells = 3,
                               depthmin.c.allcells = 5,
                               depthmin.c.nor = 0,
                               min.nhet = 0,
                               quiet = TRUE)

    expect_s4_class(muscadet_cna, "muscadet")
    expect_length(muscadet_cna$cnacalling, 21)
})


test_that("preProcSample2() returns a correct ouput", {

    data("muscadet_obj", package = "muscadet")

    counts <- muscadet_obj$cnacalling$combined.counts
    counts <- counts[complete.cases(counts),]
    counts_clus <- counts[which(counts$cluster == 1),]

    result <- preProcSample2(counts_clus)
    expect_length(result, 8)
    expect_equal(ncol(result[["pmat"]]), 10)
})


test_that("getSegConsensus() returns a correct output", {

    segs <- data.frame(
      chrom = c("chr1", "chr1", "chr2", "chr2"),
      start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
      end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
      cluster = c("1", "2", "1", "2")
    )

    consensus_segs <- getSegConsensus(segs,
                                      ncells = c("1" = 50, "2" = 30),
                                      dist.breakpoints = 1e6)
    expect_equal(dim(consensus_segs), c(2,3))
})

test_that("annotateSegments() returns a correct ouput", {

    segs <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2"),
        start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
        end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
        cluster = c("1", "2", "1", "2"),
        cf.em = c(1, 1, 1, 1),
        tcn.em = c(1, 2, 3, 3),
        lcn.em = c(0, 1, 1, 1)
    )

    consensus_segs <- getSegConsensus(
        x = segs,
        ncells = c("1" = 50, "2" = 30)
    )

    table <- annotateSegments(
        x = segs,
        consensus_segs = consensus_segs,
        ncells = c("1" = 50, "2" = 30)
    )

    expect_equal(dim(table), c(4, 19))
})
