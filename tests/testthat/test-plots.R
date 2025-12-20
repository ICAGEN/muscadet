
# heatmapMuscadet() ------------------------------------------------------------
test_that("heatmapMuscadet() generates a heatmap with show_missing=FALSE", {
    skip_on_cran()
    skip_if_not_installed("ComplexHeatmap")

    data("exdata_muscadet", package = "muscadet")

    # Create a temp file for saving the heatmap
    out_file <- tempfile(fileext = ".png")


    # show_missing = FALSE
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        partition = 0.3,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        clusters = exdata_muscadet$clustering$clusters[["0.3"]],
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        averages = TRUE,
                        filename = out_file,
                        show_missing = FALSE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # Clean up
    unlink(out_file)
})

test_that("heatmapMuscadet() generates a heatmap with show_missing=TRUE", {
    skip_on_cran()
    skip_if_not_installed("ComplexHeatmap")

    data("exdata_muscadet", package = "muscadet")

    # Create a temp file for saving the heatmap
    out_file <- tempfile(fileext = ".png")


    # show_missing = TRUE
    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        partition = 0.3,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        clusters = exdata_muscadet$clustering$clusters[["0.3"]],
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        filename = out_file,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    expect_no_error(
        heatmapMuscadet(exdata_muscadet,
                        averages = TRUE,
                        filename = out_file,
                        show_missing = TRUE,
                        quiet = TRUE)
    )
    expect_true(file.exists(out_file))

    # Clean up
    unlink(out_file)
})

test_that("heatmapMuscadet() errors if partition not found", {
    skip_on_cran()

    data("exdata_muscadet", package = "muscadet")

    expect_error(heatmapMuscadet(exdata_muscadet, partition = "a"))
})

test_that("heatmapMuscadet() errors with invalid muscadet input", {
    skip_on_cran()

    expect_error(heatmapMuscadet(data.frame(a = 1:5)))
})

# plotSil() --------------------------------------------------------------------
test_that("plotSil() returns a ggplot object for valid input", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    data("exdata_muscadet", package = "muscadet")

    expect_no_error(p <- plotSil(exdata_muscadet, partition = 0.3, title = "Silhouette Test"))

    expect_s3_class(p, "ggplot")
})

test_that("plotSil() errors for wrong inputs", {
    data("exdata_muscadet", package = "muscadet")

    expect_error(plotSil(list()))
    expect_error(plotSil(exdata_muscadet))
})

test_that("plotSil() errors when silhouette data is missing", {

    data("exdata_muscadet", package = "muscadet")
    dummy <- exdata_muscadet
    dummy@clustering$silhouette <- list()

    expect_error(plotSil(dummy))
})



# plotIndexes() ----------------------------------------------------------------
test_that("plotIndexes() returns a ggplot object - silhouette", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")

    data("exdata_muscadet", package = "muscadet")

    expect_no_error(p <- plotIndexes(
        exdata_muscadet,
        index = "silhouette",
        title = "Test Plot"
    ))

    expect_s3_class(p, "ggplot")
})

test_that("plotIndexes() returns a ggplot object - other indexes", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("fpc")
    skip_if_not_installed("clusterSim")

    data("exdata_muscadet", package = "muscadet")

        expect_no_error(p <- plotIndexes(
        exdata_muscadet,
        index = c("silhouette", "dunn2", "c"),
        title = "Test Plot"
    ))

    expect_s3_class(p, "ggplot")

    expect_no_error(p <- plotIndexes(
        exdata_muscadet,
        index = c("pearsongamma"),
        title = "Test Plot"
    ))

    expect_s3_class(p, "ggplot")
})


test_that("plotIndexes() errors with invalid muscadet object", {
    expect_error(plotIndexes(data.frame(a = 1)))
})

test_that("plotIndexes() errors with muscadet object missing clustering", {
    data("exdata_muscadet", package = "muscadet")
    dummy <- exdata_muscadet
    dummy@clustering <- list()
    expect_error(plotIndexes(dummy))
})


# plotProfile ------------------------------------------------------------------
test_that("plotProfile() runs without error for all cells", {
    skip_on_cran()

    data("exdata_muscadet", package = "muscadet")

    tmpfile <- tempfile(fileext = ".pdf")

    # all cells
    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(exdata_muscadet, data = "allcells", title = "Test plot"))
    dev.off()
    # Check file is not empty
    expect_gt(file.info(tmpfile)$size, 1000)

    # cluster 1
    pdf(tmpfile, width = 6, height = 3)
    expect_no_error(plotProfile(exdata_muscadet, data = 1, title = "Test plot"))
    dev.off()
    # Check file is not empty
    expect_gt(file.info(tmpfile)$size, 1000)

    unlink(tmpfile)
})

test_that("plotProfile() errors with wrong inputs", {
    data("exdata_muscadet", package = "muscadet")

    expect_error(plotProfile(list()))
    expect_error(plotProfile(exdata_muscadet))
    expect_error(plotProfile(exdata_muscadet, data = 17))
})


# plotCNA() --------------------------------------------------------------------

test_that("plotCNA() returns a ggplot object and runs without error", {
    skip_on_cran()

    data("exdata_muscadet", package = "muscadet")

    expect_no_error(p <- plotCNA(exdata_muscadet, title = "Copy Number Alterations in Test Data"))

    # Check that it returns a ggplot object
    expect_s3_class(p, "ggplot")

})



# plotUMAP() -------------------------------------------------------------------
test_that("plotUMAP returns a ggplot object", {
    skip_on_cran()

    data("exdata_muscadet", package = "muscadet")

    expect_true(!is.null(exdata_muscadet@clustering$umap))
    expect_true(length(exdata_muscadet@clustering$clusters) > 0)

    p <- plotUMAP(exdata_muscadet, partition = 0.3)
    expect_s3_class(p, "ggplot")
})

test_that("plotUMAP errors with wrong class", {
    expect_error(plotUMAP(list()))
})

test_that("plotUMAP detects missing clustering", {
    dummy <- exdata_muscadet
    dummy@clustering$umap <- NULL
    expect_error(plotUMAP(dummy))
})


# add_labels() -----------------------------------------------------------------
test_that("add_labels() adds text layers to ggplot", {
    skip_on_cran()

    # Base plot
    p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
        geom_point()

    # Apply add_labels
    labeled_plot <- add_labels(p, labels = cyl)

    # Check: class still ggplot
    expect_s3_class(labeled_plot, "ggplot")

    # Check: additional layer(s) were added
    expect_gt(length(labeled_plot$layers), length(p$layers))
})


# heatmapStep() ----------------------------------------------------------------

test_that("heatmapStep() generates a PNG heatmap", {
    skip_on_cran()
    skip_if_not_installed("patchwork")

    data("exdata_muscadet", package = "muscadet")
    data("exdata_muscadet_ref", package = "muscadet")
    data("exdata_peaks", package = "muscadet")

    # Example data
    expect_no_error(
        obj_atac_all <- computeLogRatioATAC(
            matTumor = matCounts(exdata_muscadet)$ATAC,
            matRef = matCounts(exdata_muscadet_ref)$ATAC,
            peaksCoord = exdata_peaks,
            genome = slot(exdata_muscadet, "genome"),
            minReads = 1,
            minPeaks = 1,
            all_steps = TRUE,
            quiet = TRUE
        )
    )

    expect_true("step08" %in% names(obj_atac_all))

    # Temp file
    heatmap_file <-tempfile(fileext = ".png")

    # Plot heatmap
    expect_no_error(
        heatmapStep(
            obj = obj_atac_all,
            step = "step08",
            filename = heatmap_file,
            title = "Test - Step08"
        )
    )

    # Test file exists and is not empty
    expect_true(file.exists(heatmap_file))
    expect_gt(file.size(heatmap_file), 0)

    # Clean up
    unlink(heatmap_file)
})




