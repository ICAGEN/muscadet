
# heatmapMuscadet() ------------------------------------------------------------

test_that("heatmapMuscadet() saves PNG with show_missing = FALSE", {
    skip_on_cran()

    # Temp file
    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    # partition
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        partition = 0.3,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # custom clusters
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        clusters = exdata_muscadet$clustering$clusters[["0.3"]],
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # assigned clusters (partition = NULL, clusters = NULL)
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # averages
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        averages = TRUE,
                        filename = out_file,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))
})

test_that("heatmapMuscadet() saves PNG with show_missing = TRUE", {
    skip_on_cran()

    # Temp file
    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    # partition
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        partition = 0.3,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # custom clusters
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        clusters = exdata_muscadet$clustering$clusters[["0.3"]],
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # assigned clusters
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # averages
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        averages = TRUE,
                        filename = out_file,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))
})

test_that("heatmapMuscadet() saves PDF", {
    skip_on_cran()

    out_file <- tempfile(fileext = ".pdf")
    on.exit(unlink(out_file))

    expect_no_error(
        heatmapMuscadet(
            exdata_muscadet,
            filename = out_file,
            partition = 0.3,
            quiet = TRUE
        )
    )
    expect_true(file.exists(out_file))
    expect_gt(file.size(out_file), 0)
})

test_that("heatmapMuscadet() saves SVG", {
    skip_on_cran()
    skip_on_ci()

    out_file <- tempfile(fileext = ".svg")
    on.exit(unlink(out_file))

    expect_no_error(
        heatmapMuscadet(
            exdata_muscadet,
            filename = out_file,
            partition = 0.3,
            quiet = TRUE
        )
    )
    expect_true(file.exists(out_file))
})

test_that("heatmapMuscadet() returns plot/width/height when no filename", {
    skip_on_cran()

    ht <- heatmapMuscadet(exdata_muscadet, partition = 0.3, quiet = TRUE)

    expect_named(ht, c("plot", "width", "height"))
    expect_s3_class(ht$plot, "gTree")
    expect_true(is.numeric(ht$width))
    expect_true(is.numeric(ht$height))
})

test_that("heatmapMuscadet() errors if partition not found", {
    skip_on_cran()

    data("exdata_muscadet", package = "muscadet")

    expect_error(heatmapMuscadet(exdata_muscadet, partition = "a"))
})

test_that("heatmapMuscadet() dim_scale affects output dimensions", {
    skip_on_cran()

    ht_full <- heatmapMuscadet(exdata_muscadet, partition = 0.3,
                               dim_scale = 1,   quiet = TRUE)
    ht_half <- heatmapMuscadet(exdata_muscadet, partition = 0.3,
                               dim_scale = 0.5, quiet = TRUE)

    expect_lt(as.numeric(ht_half$width),  as.numeric(ht_full$width))
    expect_lt(as.numeric(ht_half$height), as.numeric(ht_full$height))
})

test_that("heatmapMuscadet() runs with quiet = FALSE", {
    skip_on_cran()

    expect_message(heatmapMuscadet(exdata_muscadet, partition = 0.3, quiet = FALSE))
})

test_that("heatmapMuscadet() runs with quiet = TRUE", {
    skip_on_cran()

    expect_no_message(heatmapMuscadet(exdata_muscadet, partition = 0.3, quiet = TRUE))
})

test_that("heatmapMuscadet() accepts white_scale as list", {
    skip_on_cran()

    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    n_omics <- length(exdata_muscadet@omics)
    ws_list <- rep(list(c(0.3, 0.7)), n_omics)

    expect_no_error(
        heatmapMuscadet(
            exdata_muscadet,
            filename = out_file,
            partition = 0.3,
            white_scale = ws_list,
            quiet = TRUE
        )
    )
})

test_that("heatmapMuscadet() accepts custom colors", {
    skip_on_cran()

    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    expect_no_error(
        heatmapMuscadet(
            exdata_muscadet,
            filename = out_file,
            partition = 0.3,
            colors = c("tomato", "steelblue"),
            quiet = TRUE
        )
    )
})

test_that("heatmapMuscadet() accepts add_bulk_lrr = FALSE", {
    skip_on_cran()

    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    expect_no_error(
        heatmapMuscadet(
            exdata_muscadet,
            filename = out_file,
            partition = 0.3,
            add_bulk_lrr = FALSE,
            quiet = TRUE
        )
    )
})

test_that("heatmapMuscadet() errors with invalid muscadet input", {
    skip_on_cran()

    expect_error(heatmapMuscadet(data.frame(a = 1:5)))
})

test_that("heatmapMuscadet() errors with invalid inputs", {
    skip_on_cran()

    # partition not found
    expect_error(heatmapMuscadet(exdata_muscadet, partition = "a"))

    # invalid filename extension
    out_file <- tempfile(fileext = ".txt")
    on.exit(unlink(out_file))
    expect_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = tempfile(fileext = ".txt"),
                        partition = 0.3)
    )

    # invalid show_missing
    expect_error(
        heatmapMuscadet(exdata_muscadet,
                        partition = 0.3,
                        show_missing = "yes")
    )

    # invalid white_scale
    expect_error(
        heatmapMuscadet(exdata_muscadet,
                        partition = 0.3,
                        white_scale = c(1.5, 0.3))
    )
    expect_error(
        heatmapMuscadet(exdata_muscadet,
                        partition = 0.3,
                        white_scale = "invalid")
    )

    # invalid scale
    expect_error(
        heatmapMuscadet(exdata_muscadet, partition = 0.3, dim_scale = -1)
    )
    expect_error(
        heatmapMuscadet(exdata_muscadet, partition = 0.3, dim_scale = 0)
    )

    # invalid raster_quality
    expect_error(
        heatmapMuscadet(exdata_muscadet, partition = 0.3, raster_quality = 0)
    )

    # row_annots not a list of HeatmapAnnotation
    expect_error(
        heatmapMuscadet(exdata_muscadet,
                        partition = 0.3,
                        row_annots = list("not_an_annotation"))
    )
})

# plotSil() --------------------------------------------------------------------

test_that("plotSil() returns a ggplot object", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    expect_no_error(p <- plotSil(exdata_muscadet, partition = 0.3, title = "Silhouette Test"))
    expect_s3_class(p, "ggplot")
})

test_that("plotSil() works with annotations = FALSE", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    expect_no_error(p <- plotSil(
        exdata_muscadet,
        partition = 0.3,
        annotations = FALSE
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotSil() accepts custom colors", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    expect_no_error(p <- plotSil(
        exdata_muscadet,
        partition = 0.3,
        colors = c("red", "blue", "green")
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotSil() errors for wrong inputs", {
    expect_error(plotSil(list()))
    expect_error(plotSil(exdata_muscadet))
    expect_error(plotSil(exdata_muscadet, partition = "z"))
})

test_that("plotSil() errors when silhouette data is missing", {
    dummy <- exdata_muscadet
    dummy@clustering$silhouette <- list()
    expect_error(plotSil(dummy))
    expect_error(plotSil(dummy, partition = 0.3))
})


# plotIndexes() ----------------------------------------------------------------

test_that("plotIndexes() returns ggplot for silhouette", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    expect_no_error(p <- plotIndexes(exdata_muscadet, index = "silhouette", title = "Test"))
    expect_s3_class(p, "ggplot")
})

test_that("plotIndexes() returns ggplot for multiple indexes", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("fpc")
    skip_if_not_installed("clusterSim")

    expect_no_error(p <- plotIndexes(exdata_muscadet, index = c("silhouette", "dunn2", "c")))
    expect_s3_class(p, "ggplot")

    expect_no_error(p <- plotIndexes(exdata_muscadet, index = c("pearsongamma")))
    expect_s3_class(p, "ggplot")
})

test_that("plotIndexes() accepts custom colors", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    expect_no_error(p <- plotIndexes(
        exdata_muscadet,
        index = "silhouette",
        colors = c("red")
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotIndexes() errors with invalid muscadet object", {
    expect_error(plotIndexes(data.frame(a = 1)))
})

test_that("plotIndexes() errors with muscadet object missing clustering", {
    dummy <- exdata_muscadet
    dummy@clustering <- list()
    expect_error(plotIndexes(dummy))
})


# plotProfile() ----------------------------------------------------------------

test_that("plotProfile() runs for allcells", {
    skip_on_cran()

    tmpfile <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmpfile))

    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(exdata_muscadet, data = "allcells", title = "Test allcells"))
    dev.off()
    expect_gt(file.info(tmpfile)$size, 1000)
})

test_that("plotProfile() runs for a specific cluster", {
    skip_on_cran()

    tmpfile <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmpfile))

    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(exdata_muscadet, data = 1, title = "Cluster 1"))
    dev.off()
    expect_gt(file.info(tmpfile)$size, 1000)
})

test_that("plotProfile() runs with allelic.type = 'lor'", {
    skip_on_cran()

    tmpfile <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmpfile))

    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(
        exdata_muscadet,
        data = "allcells",
        allelic.type = "lor"
    ))
    dev.off()
})

test_that("plotProfile() runs with var.colors = 'none'", {
    skip_on_cran()

    tmpfile <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmpfile))

    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(exdata_muscadet, data = "allcells", var.colors = "none"))
    dev.off()
})

test_that("plotProfile() errors with invalid inputs", {
    expect_error(plotProfile(list()))
    expect_error(plotProfile(exdata_muscadet))
    expect_error(plotProfile(exdata_muscadet, data = 99))
    expect_error(plotProfile(
        exdata_muscadet,
        data = "allcells",
        allelic.type = "invalid"
    ))
    expect_error(plotProfile(exdata_muscadet, data = "allcells", point.cex = "a"))
})


# plotCNA() --------------------------------------------------------------------

test_that("plotCNA() returns ggplot with default args", {
    skip_on_cran()
    expect_no_error(p <- plotCNA(exdata_muscadet, title = "CNA Test"))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() labels argument - keyword modes", {
    skip_on_cran()

    expect_no_error(p <- plotCNA(exdata_muscadet, labels = "auto"))
    expect_s3_class(p, "ggplot")

    expect_no_error(p <- plotCNA(exdata_muscadet, labels = "clusters"))
    expect_s3_class(p, "ggplot")

    expect_no_error(p <- plotCNA(exdata_muscadet, labels = "cells"))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() labels argument - NULL", {
    skip_on_cran()

    expect_no_error(p <- plotCNA(exdata_muscadet, labels = NULL))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() labels argument - custom vector", {
    skip_on_cran()

    n_clusters <- length(unique(exdata_muscadet@cnacalling$ncells))
    custom_labels <- paste0("clone_", seq_len(n_clusters))

    expect_no_error(p <- plotCNA(exdata_muscadet, labels = custom_labels))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() labels argument errors", {
    skip_on_cran()

    # wrong length vector
    expect_error(plotCNA(exdata_muscadet, labels = c("a", "b", "c", "d", "e",
                                                     "f", "g", "h", "i", "j")))
    # invalid type
    expect_error(plotCNA(exdata_muscadet, labels = TRUE))
})

test_that("plotCNA() cf.gradient = TRUE", {
    skip_on_cran()

    expect_no_error(p <- plotCNA(exdata_muscadet, cf.gradient = TRUE))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() cna.colors - named vector with missing state", {
    skip_on_cran()

    # Missing cnloh — should be filled with default
    expect_no_error(p <- plotCNA(
        exdata_muscadet,
        cna.colors = c("gain" = "red", "loss" = "blue")
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() cna.colors - unnamed vector", {
    skip_on_cran()

    expect_no_error(p <- plotCNA(exdata_muscadet, cna.colors = c("red", "blue", "green")))
    expect_s3_class(p, "ggplot")

    # Unnamed with only 2 colors — padded with default
    expect_no_error(p <- plotCNA(exdata_muscadet, cna.colors = c("red", "blue")))
    expect_s3_class(p, "ggplot")
})

test_that("plotCNA() errors with invalid object", {
    expect_error(plotCNA(list()))
    expect_error(plotCNA(data.frame(a = 1:5)))
})

# plotUMAP() -------------------------------------------------------------------

test_that("plotUMAP() returns ggplot with partition", {
    skip_on_cran()

    expect_no_error(p <- plotUMAP(exdata_muscadet, partition = 0.3))
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP() returns ggplot with assigned clusters", {
    skip_on_cran()

    expect_no_error(p <- plotUMAP(exdata_muscadet))
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP() returns ggplot with custom clusters", {
    skip_on_cran()

    clus <- exdata_muscadet$clustering$clusters[["0.3"]]
    expect_no_error(p <- plotUMAP(exdata_muscadet, clusters = clus))
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP() returns ggplot with add_clusters_labels = TRUE", {
    skip_on_cran()
    skip_if_not_installed("ggrepel")

    expect_no_error(p <- plotUMAP(
        exdata_muscadet,
        partition = 0.3,
        add_clusters_labels = TRUE
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP() accepts custom colors", {
    skip_on_cran()

    expect_no_error(p <- plotUMAP(
        exdata_muscadet,
        partition = 0.3,
        colors = c("red", "blue", "green")
    ))
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP errors with wrong input class", {
    expect_error(plotUMAP(list()))
})

test_that("plotUMAP() errors with invalid inputs", {
    expect_error(plotUMAP(list()))

    # missing umap data
    dummy <- exdata_muscadet
    dummy@clustering$umap <- NULL
    expect_error(plotUMAP(dummy))

    # partition not found
    expect_error(plotUMAP(exdata_muscadet, partition = "z"))

    # no assigned clusters
    dummy2 <- exdata_muscadet
    dummy2@cnacalling$clusters <- NULL
    expect_error(plotUMAP(dummy2))
})


# add_labels() -----------------------------------------------------------------

test_that("add_labels() with repel = TRUE adds layers", {
    skip_on_cran()
    skip_if_not_installed("ggrepel")

    p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
        geom_point()

    labeled <- add_labels(p, labels = cyl, repel = TRUE)
    expect_s3_class(labeled, "ggplot")
    expect_gt(length(labeled$layers), length(p$layers))
})

test_that("add_labels() with repel = FALSE adds layers", {
    skip_on_cran()

    p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
        geom_point()

    labeled <- add_labels(p, labels = cyl, repel = FALSE)
    expect_s3_class(labeled, "ggplot")
    expect_gt(length(labeled$layers), length(p$layers))
})

test_that("add_labels() with label.box = TRUE", {
    skip_on_cran()

    p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
        geom_point()

    labeled <- add_labels(p, labels = cyl, repel = FALSE, label.box = TRUE)
    expect_s3_class(labeled, "ggplot")
})

test_that("add_labels() with repel = TRUE and label.box = TRUE", {
    skip_on_cran()
    skip_if_not_installed("ggrepel")

    p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
        geom_point()

    labeled <- add_labels(p, labels = cyl, repel = TRUE, label.box = TRUE)
    expect_s3_class(labeled, "ggplot")
})

# heatmapStep() ----------------------------------------------------------------

# Shared setup for heatmapStep tests
.build_heatmapstep_obj <- function() {
    data("exdata_muscadet",     package = "muscadet")
    data("exdata_muscadet_ref", package = "muscadet")
    data("exdata_peaks",        package = "muscadet")

    computeLogRatioATAC(
        matTumor  = matCounts(exdata_muscadet)$ATAC,
        matRef    = matCounts(exdata_muscadet_ref)$ATAC,
        peaksCoord = exdata_peaks,
        genome    = slot(exdata_muscadet, "genome"),
        minReads  = 1,
        minPeaks  = 1,
        all_steps = TRUE,
        quiet     = TRUE
    )
}

test_that("heatmapStep() saves PNG", {
    skip_on_cran()
    skip_if_not_installed("patchwork")

    obj <- .build_heatmapstep_obj()
    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    expect_no_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = out_file,
        title = "Test PNG"
    ))
    expect_true(file.exists(out_file))
    expect_gt(file.size(out_file), 0)
})

test_that("heatmapStep() saves PDF", {
    skip_on_cran()
    skip_if_not_installed("patchwork")

    obj <- .build_heatmapstep_obj()
    out_file <- tempfile(fileext = ".pdf")
    on.exit(unlink(out_file))

    expect_no_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = out_file
    ))
    expect_true(file.exists(out_file))
    expect_gt(file.size(out_file), 0)
})

test_that("heatmapStep() saves SVG", {
    skip_on_cran()
    skip_on_ci()
    skip_if_not_installed("patchwork")

    obj <- .build_heatmapstep_obj()
    out_file <- tempfile(fileext = ".svg")
    on.exit(unlink(out_file))

    expect_no_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = out_file
    ))
    expect_true(file.exists(out_file))
})

test_that("heatmapStep() works with col_breaks instead of col_quantiles", {
    skip_on_cran()
    skip_if_not_installed("patchwork")

    obj <- .build_heatmapstep_obj()
    out_file <- tempfile(fileext = ".png")
    on.exit(unlink(out_file))

    expect_no_error(
        heatmapStep(
            obj = obj,
            step = "step08",
            filename      = out_file,
            col_quantiles = NULL,
            col_breaks    = c(-2, -0.5, 0.5, 2)
        )
    )
    expect_true(file.exists(out_file))
})

test_that("heatmapStep() errors with invalid inputs", {
    skip_on_cran()
    skip_if_not_installed("patchwork")

    obj <- .build_heatmapstep_obj()

    # not a list
    expect_error(heatmapStep(
        obj = "not_a_list",
        step = "step08",
        filename = tempfile(fileext = ".png")
    ))

    # step not found
    expect_error(heatmapStep(
        obj = obj,
        step = "stepXX",
        filename = tempfile(fileext = ".png")
    ))

    # invalid filename extension
    expect_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = tempfile(fileext = ".txt")
    ))

    # neither col_quantiles nor col_breaks
    expect_error(
        heatmapStep(
            obj = obj,
            step = "step08",
            filename = tempfile(fileext = ".png"),
            col_quantiles = NULL,
            col_breaks = NULL
        )
    )

    # invalid colors length
    expect_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = tempfile(fileext = ".png"),
        colors = c("red", "blue")
    ))

    # invalid col_quantiles values
    expect_error(heatmapStep(
        obj = obj,
        step = "step08",
        filename = tempfile(fileext = ".png"),
        col_quantiles = c(0.1, 0.4, 0.6, 1.5)
    ))
})




