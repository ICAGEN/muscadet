

test_that("assignClusters() returns an updated muscadet object", {

    data("muscadet_obj", package = "muscadet")
    # remove cna calling result
    muscadet_obj@cnacalling <- list()

    # Select clustering result for partition 0.6
    obj1 <- assignClusters(muscadet_obj, partition = 0.6)

    expect_true(length(obj1@cnacalling) != 0)
    expect_identical(length(table(obj1@cnacalling$clusters)), as.integer(2)) # 2 clusters

    # Assign custom clusters
    cell_names <- Reduce(union, SeuratObject::Cells(muscadet_obj))
    n1 <- sample(1:length(cell_names), 1)
    n2 <- length(cell_names) - n1
    custom_clusters <- c(rep.int(1, n1), rep.int(2, n2))
    names(custom_clusters) <- cell_names

    obj2 <- assignClusters(muscadet_obj, clusters = custom_clusters, redo_imputation = FALSE)

    expect_true(length(obj2@cnacalling) != 0)
    expect_length(table(obj2@cnacalling$clusters), 2) # 2 clusters

    # Assign clusters with remapping
    # example to remap from 5 clusters to 4 by merging clusters 1 and 2
    clusters <- muscadet_obj$clustering$clusters[["1"]] # res = 1
    mapping <- c("1" = 1, "2" = 1, "3" = 2, "4" = 3, "5" = 4)

    obj3 <- assignClusters(muscadet_obj, clusters = clusters, mapping = mapping)

    expect_true(length(obj3@cnacalling) != 0)
    expect_length(table(obj3@cnacalling$clusters), 4) # 4 clusters

    obj4 <- assignClusters(muscadet_obj, partition = 1, mapping = mapping)

    expect_true(length(obj4@cnacalling) != 0)
    expect_length(table(obj4@cnacalling$clusters), 4) # 4 clusters
})


test_that("mergeCounts() returns an updated muscadet object with data frames in the cnacalling slot", {

    # Load example muscadet objects
    data(muscadet_obj)
    data(muscadet_obj_ref)

    # remove cna calling result
    muscadet_obj@cnacalling <- list()
    muscadet_obj <- assignClusters(muscadet_obj, partition = 0.6)

    # Merge counts from all omics from both sample and reference
    obj <- mergeCounts(muscadet_obj, muscadet_obj_ref)

    expect_true(length(obj@cnacalling) != 1)
    expect_identical(class(obj@cnacalling$allelic.counts), "data.frame")
    expect_identical(class(obj@cnacalling$coverage.counts), "data.frame")
    expect_identical(class(obj@cnacalling$combined.counts), "data.frame")

})


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
