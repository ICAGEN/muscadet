#' Heatmap plot for `muscadet` object
#'
#' This function generates a heatmap to visualize log R ratio (LRR) data
#' contained in [`muscadet`][muscadet-class] objects. One heatmap is generated
#' per omic, rows are cells and columns are chromosomes, for
#' [`muscadet`][muscadet-class] object containing multiple omics, the heatmaps
#' are plotted horizontally aligned. The cells can be clustered for a specific
#' clustering partition following the clustering step of
#' [`muscadet`][muscadet-class] object, or custom cluster assignments can be
#' used. Additionally, LRR values from bulk sequencing data can be plotted as an
#' annotation under the heatmaps.
#'
#' @param x A [`muscadet`][muscadet-class] object containing LRR data for all
#'   omics (with [computeLogRatio()]) and clustering data (with
#'   [clusterMuscadet()]) (`muscadet`).
#'
#' @param filename (Optional) Character string specifying the file path to save
#'   the heatmap image in the PNG (if it ends by .png), PDF (if it ends by
#'   .pdf) or SVG (if it ends by .svg) format (`character` string).
#'
#' @param partition (Optional) Value specifying the clustering partition to
#'   plot (`numeric` or `character`). It should be either the resolution or the
#'   number of cluster (k) used for clustering depending on the clustering
#'   method (`res_range` or `k_range` with [clusterMuscadet()]).
#'   If both `partition` and `clusters` arguments are `NULL` (default), the
#'   assigned clusters for CNA calling (`x@cnacalling$cluster`) are used if
#'   available in `x` (see [assignClusters()]).
#'
#' @param clusters (Optional) A custom named vector of cluster assignments
#'   (`integer` named vector). Names must corresponds to cell names within the
#'   muscadet object `x`. If it contains less cells than the muscadet object
#'   `x`, the missing cells are filtered out and not displayed in the heatmap.
#'   If `show_missing = FALSE` only the provided cells with data in all omics
#'   will be displayed. If both `partition` and `clusters` arguments are `NULL`
#'   (default), the assigned clusters for CNA calling (`x@cnacalling$cluster`)
#'   are used if available in `x` (see [assignClusters()]).
#'
#' @param add_bulk_lrr Logical. If `TRUE` (default), adds bulk log R ratio (LRR) data as
#'   annotation if available in the muscadet object.
#'
#' @param show_missing Logical. If `TRUE` (default), missing cells (i.e., cells
#'   with missing data in at least one omic) are displayed in the heatmaps.
#'
#' @param averages Logical. If `TRUE`, plots the average log R ratio per
#'   cluster. Default is `FALSE`.
#'
#' @param title Character string for the title of the plot (`character`
#'   string). Default is an empty character string.
#'
#' @param row_annots Optional. A list of
#'   [`HeatmapAnnotation`][ComplexHeatmap::HeatmapAnnotation-class] objects from
#'   the [`ComplexHeatmap`][ComplexHeatmap::ComplexHeatmap-package] package,
#'   specifying row annotations to add on left part of the heatmap. Each element
#'   in the list must be of class
#'   [`HeatmapAnnotation`][ComplexHeatmap::HeatmapAnnotation-class], must be a row
#'   annotation (using [`rowAnnotation()`][ComplexHeatmap::rowAnnotation()] or
#'   [`HeatmapAnnotation()`][ComplexHeatmap::HeatmapAnnotation()] with
#'   `which = 'row'`), and must have a unique name (`name` argument in
#'   [`rowAnnotation()`][ComplexHeatmap::rowAnnotation()] or
#'   [`HeatmapAnnotation()`][ComplexHeatmap::HeatmapAnnotation()]). If
#'   `averages =FALSE`, annotations must concern cells, while if `averages = TRUE`
#'   they must concern clusters. Default is `NULL`, no row annotations is added.
#'
#' @param white_scale Numeric vector of length 2 or a list of numeric vectors
#'   (`numeric` vector or `list`).
#' - If a numeric vector of length 2, the same white color boundaries are
#'   applied to all omics in the muscadet object. E.g. `c(0.3, 0.7)` (default).
#' - If a list (named or not with omics name), it must have the same length as
#'   the number of omics in the muscadet object, where each vector element
#'   applies the white color boundaries for a specific omic. E.g. `list(c(0.3,
#'   0.7), c(0.4, 0.6))` uses 0.3 and 0.7 quantiles of LRR ref data for the 1st
#'   omic heatmap, and 0.4 and 0.6 quantiles for the second.
#'
#'   Values of the vectors must be between 0 and 1, specifying the quantiles of
#'   the LRR reference data that define the boundaries for the white color in
#'   each heatmap. LRR values falling within this range are considered close to
#'   the majority of the LRR reference data, indicating no significant gain or
#'   loss of coverage, and are represented as white on the heatmap. Default is
#'   `c(0.3, 0.7)`.
#'
#' @param colors Vector of colors for the cluster annotation (`character`
#'   vector). If `NULL` (default), it uses predefined colors.
#'
#' @param png_res Resolution in ppi for [grDevices::png()] if `filename` ends
#'   with the .png extension (`numeric`). Default is `300`.
#'
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return A list containing:
#' - `plot`: A [`gTree`][grid::gTree] object created with [grid::grid.grab()]
#' ([`gTree`][grid::gTree]).
#' - `width`: Width of the heatmap plot in mm ([`unit`][grid::unit()]).
#' - `height`: Height of the heatmap plot in mm ([`unit`][grid::unit()]).
#'
#' If the `filename` argument is provided, the heatmap is directly saved as a
#' PNG image at the provided path.
#'
#' @include objects.R
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom methods slot
#' @importFrom stats median
#' @importFrom grDevices pdf png svg palette dev.off
#' @importFrom grid gpar unit grid.rect grid.text grid.grab
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#'
#' # --- Method "seurat" ---
#'
#' print(exdata_muscadet$clustering$params$method)
#'
#' # Perform clustering if not already done
#' # exdata_muscadet <- clusterMuscadet(
#' #     x = exdata_muscadet,
#' #     method = "seurat",
#' #     res_range = c(0.1, 0.3),
#' #     dims_list = list(1:8, 1:8),
#' #     knn_seurat = 10, # adapted for low number of cells in example data
#' #     knn_range_seurat = 30 # adapted for low number of cells in example data
#' # )
#'
#' # Plot a single partition
#' heatmapMuscadet(exdata_muscadet,
#'                 filename = file.path("heatmap_res0.3.png"),
#'                 partition = 0.3,
#'                 show_missing = FALSE) # only displaying cells without missing data
#'
#' ht <- heatmapMuscadet(exdata_muscadet, partition = 0.3)
#' pdf(
#'     file = file.path("heatmap_res0.3.pdf"),
#'     width = ht$width * 0.0393701, # convert to inches
#'     height = ht$height * 0.0393701, # convert to inches
#' )
#' grid.draw(ht$plot)
#' dev.off()
#'
#'
#' # Loop over partitions
#' for (p in names(exdata_muscadet$clustering$clusters)) {
#'     filename <- paste0("heatmap_res", p, ".png")
#'     title <- paste(
#'         "Example |",
#'         paste0("method=", exdata_muscadet$clustering$params[["method"]]), "|",
#'         paste0("omics=", paste0(exdata_muscadet$clustering$params[["omics"]], collapse = ",")), "|",
#'         paste0("dims=", "1:10,1:10"), "|",
#'         paste0("res=", p)
#'     )
#'     heatmapMuscadet(exdata_muscadet, filename, partition = p, title = title)
#' }
#'
#' # --- Plot Averages per Clusters ---
#'
#' heatmapMuscadet(exdata_muscadet,
#'                 filename = file.path("heatmap_res0.3_averages.png"),
#'                 partition = 0.3,
#'                 averages = TRUE,
#'                 add_bulk_lrr = FALSE)
#'
#' # --- Add Row Annotation ---
#'
#' library("ComplexHeatmap")
#' library("grid")
#'
#' # Define example cell annotation
#' muscadet_cells <- Reduce(union, SeuratObject::Cells(exdata_muscadet))
#' cells_origin <- setNames(c(
#'     rep("sample1", ceiling(length(muscadet_cells) / 2)),
#'     rep("sample2", floor(length(muscadet_cells) / 2))
#'     ),
#'     muscadet_cells
#' )
#' cells_origin <- cells_origin[sort(names(cells_origin))]
#' # IMPORTANT: annotation names (cells) must be sorted to match heatmap
#' # matrices (column names of log ratio matrix are sorted in heatmapMuscadet())
#'
#' # Create row annotation
#' ha <- rowAnnotation(
#'     annot = anno_simple(
#'         cells_origin[sort(names(cells_origin))],
#'         col = c("sample1" = "cadetblue3", "sample2" = "orchid3")),
#'     name = "origin", # unique name
#'     annotation_label = "origin", # label displayed on heatmap
#'     annotation_name_gp = gpar(fontsize = 10) # change font size
#' )
#'
#' # Plot heatmap with supplementary row annotation
#' heatmapMuscadet(exdata_muscadet,
#'                 filename = file.path("heatmap_res0.3_annot.png"),
#'                 partition = 0.3,
#'                 row_annots = list(ha))
#'
#' # --- Add Row Annotation for averages ---
#'
#' library("ComplexHeatmap")
#' library("grid")
#'
#' # Define example cluster annotation
#' clus <- setNames(c("annot1", "annot2"), c(1, 2)) # 2 clusters for partition 0.6
#' clus <- clus[order(names(clus))]
#' # IMPORTANT: annotation names (clusters) must be sorted to match heatmap
#' # matrices (column names of log ratio averages matrix are sorted in heatmapMuscadet())
#'
#'
#' # Create row annotation
#' ha2 <- rowAnnotation(
#'     annot = anno_simple(clus, col = c("annot1" = "tomato", "annot2" = "gold2")),
#'     name = "annot", # unique name
#'     annotation_label = "annot", # label displayed on heatmap
#'     annotation_name_gp = gpar(fontsize = 10) # change font size
#' )
#' heatmapMuscadet(exdata_muscadet,
#'                 averages = TRUE,
#'                 filename = file.path("heatmap_res0.3_annot_averages.png"),
#'                 partition = 0.3,
#'                 row_annots = list(ha2))
#'
#'
#' # --- Method "hclust" ---
#'
#' # Perform clustering if not already done
#' exdata_muscadet2 <- clusterMuscadet(
#'     x = exdata_muscadet,
#'     method = "hclust",
#'     k_range = 2:4,
#'     dist_method = "euclidean",
#'     hclust_method = "ward.D"
#' )
#'
#' print(exdata_muscadet2$clustering$params$method)
#'
#' # Plot a single partition
#' heatmapMuscadet(exdata_muscadet2,
#'                 filename = file.path("heatmap_k2.png"),
#'                 partition = 2,
#'                 show_missing = FALSE)
#'
#' # Loop over partitions
#' for (p in names(exdata_muscadet2$clustering$clusters)) {
#'
#'     filename <- paste0("heatmap_k", p, ".png")
#'     title <- paste(
#'         "Example |",
#'         paste0("method=", exdata_muscadet2$clustering$params[["method"]]), "|",
#'         exdata_muscadet2$clustering$params[["dist_method"]],
#'         exdata_muscadet2$clustering$params[["hclust_method"]], "|",
#'         paste0("weights=",
#'                paste0(exdata_muscadet2$clustering$params[["weights"]], collapse = ",")),
#'         "|",
#'         paste0("k=", p)
#'     )
#'
#'     heatmapMuscadet(exdata_muscadet2, filename, partition = p, title = title)
#' }
#' }
#'
heatmapMuscadet <- function(x,
                            filename = NULL,
                            partition = NULL,
                            clusters = NULL,
                            add_bulk_lrr = TRUE,
                            show_missing = TRUE,
                            averages = FALSE,
                            title = "",
                            row_annots = NULL,
                            white_scale = c(0.3, 0.7),
                            colors = NULL,
                            png_res = 300,
                            quiet = FALSE) {

    ## x ------
    # Validate input: x must be a muscadet object
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Validate log ratio present in the muscadet object
    stopifnot("Input object `x` must contain log ratio data (use computeLogRatio())." =
                  !is.null(matLogRatio(x)))

    ## filename ------
    # Check if output directory exists
    if (!is.null(filename)) {
        stopifnot("`filename`: the directory doesn't exist" = file.exists(dirname(filename)))
        stopifnot("The `filename` argument must end with either .png, .pdf or svg." = grepl(".(png|pdf|svg)$", filename))
    }

    ## Validate clusters & partition ------
    # Get common and all cells
    common_cells <- sort(Reduce(intersect, lapply(muscadet::matLogRatio(x), rownames)))
    all_cells <- sort(Reduce(union, lapply(muscadet::matLogRatio(x), rownames)))

    # Filter cells based on provided `clusters` argument
    if (!is.null(clusters)) {

        # Check `clusters` cells
        stopifnot("The `clusters` argument must have names, corresponding to cell names within the muscadet object `x`." = !is.null(names(clusters)))

        stopifnot(
            "Names of `clusters` don't match cell names within the muscadet object `x`." = length(setdiff(names(clusters), all_cells)) == 0
        )

        # cells not in `clusters`
        cells_filtered <- setdiff(all_cells, names(clusters))

        if (length(cells_filtered) > 0) {
            warning(
                paste(
                    "The `clusters` argument does not contain cluster assignments for all cells.",
                    length(cells_filtered),
                    "cells in the muscadet object `x` are filtered out."
                )
            )
            # filter cells based on `clusters` cells
            common_cells <- common_cells[common_cells %in% names(clusters)]
            all_cells <- all_cells[all_cells %in% names(clusters)]
        }

    } else if (!is.null(partition)) {
        stopifnot(
            "The muscadet object `x` must contain the clustering results for the provided `partition`." =
                as.character(partition) %in% names(x@clustering$clusters)
        )

    } else if (is.null(partition) && is.null(clusters)) {

        stopifnot(
            "The muscadet object `x` must contain the assigned clusters (x@cnacalling$clusters) (use assignClusters())
            if `partition` and `clusters` arguments are NULL." =
                !is.null(x@cnacalling$clusters)
        )
    }

    ## Validate show_missing ------
    stopifnot("`show_missing` must be a single TRUE/FALSE value." =
                  is.logical(show_missing) && length(show_missing) == 1 && !is.na(show_missing))
    # Set to no missing cells if only one omic
    if (length(x@omics) == 1) {
        show_missing <- FALSE
    }

    ## Validate add_bulk_lrr ------
    # Default addition of bulk data depending on its presence in muscadet object
    if (add_bulk_lrr | !is.logical(show_missing)) {
        add_bulk_lrr <- !is.null(x@bulk.data$logratio)
    }

    ## Validate row_annots ------
    if (!is.null(row_annots)) {

        # Check that it is a list of HeatmapAnnotation objects
        if (!is.list(row_annots) || !all(sapply(row_annots, function(x) inherits(x, "HeatmapAnnotation")))) {
            stop("`row_annots` must be a list of HeatmapAnnotation objects.")
        }

        # Ensure all annotations are row annotations
        if (!all(sapply(row_annots, function(x) x@anno_list[["annot"]]@fun@which == "row"))) {
            stop("All elements in `row_annots` must be row annotations (use rowAnnotation() or HeatmapAnnotation(..., which = 'row').")
        }

        # Ensure annotation names are unique
        annot_names <- sapply(row_annots, function(x) x@name)
        if (length(annot_names) != length(unique(annot_names))) {
            stop("All annotation names in `row_annots` must be unique (`name` argument in rowAnnotation() or HeatmapAnnotation()).")
        }
    }

    ## Validate white_scale ------
    if (is.numeric(white_scale) && length(white_scale) == 2) {
        stopifnot("`white_scale` must be a numeric vector of length 2." = length(white_scale) == 2)
        stopifnot(
            "`white_scale` values must be between 0 and 1." = all(white_scale >= 0 & white_scale <= 1)
        )
        # Single pair of values applied to all omics
        white_scale <- round(sort(white_scale), 2)
        white_scale <- rep(list(white_scale), length(x@omics))  # Ensure list matches omics count
        names(white_scale) <- names(x@omics)  # Assign omic names for consistency

    } else if (is.list(white_scale)) {
        stopifnot(
            "`white_scale` must have the same length as the number of omics in muscadet object `x`." =
                length(white_scale) == length(x@omics)
        )

        stopifnot("All elements of `white_scale` must be numeric vectors of length 2." =
                      all(
                          vapply(white_scale, function(x)
                              is.numeric(x) && length(x) == 2, logical(1))
                      ))

        white_scale <- lapply(white_scale, function(x) {
            x <- round(x, 2)
            stopifnot("The two elements of `white_scale` vectors must not be equal." = x[1] != x[2])
            if (x[1] > x[2])
                x <- sort(x)
            return(x)
        })

        # If white_scale is unnamed, assign names based on omics order
        if (is.null(names(white_scale))) {
            names(white_scale) <- names(x@omics)
        } else {
            stopifnot(
                "Named `white_scale` list must have names matching the omics in muscadet object `x`." =
                    all(names(white_scale) %in% names(x@omics))
            )
        }
    } else {
        stop(
            "`white_scale` must be either a numeric vector of length 2 or a list of such vectors."
        )
    }

    # Set default color palette for clusters if not provided
    if (is.null(colors)) {
        colors <- rep(
            c("#FABC2A", "#7FC97F", "#EE6C4D", "#39ADBD", "#BEAED4", "#FEE672", "#F76F8E",
              "#487BEA", "#B67BE6", "#F38D68", "#7FD8BE", "#F2AFC0"), 3)
    }

    # Set palette
    palette(colors)

    # Averages ------
    if (averages) {

        # retrieve clusters if `partition` is defined and not `clusters`
        if (is.null(clusters) & !is.null(partition)) {
            clusters <- x@clustering$clusters[[as.character(partition)]]
        } else if (is.null(clusters) & is.null(partition)) {
            clusters <- x@cnacalling$clusters
        }

        # Order clusters in the same way as the matrix
        clusters <- setNames(as.vector(clusters), names(clusters))
        clusters_order <- unique(clusters) # save original order

        # Modify matrix of log ratio into a matrix of average log ratio per cluster
        for (omic in names(x@omics)) {

            # Retrieve log ratio matrix
            mat_lrr <- matLogRatio(x@omics[[omic]])
            mat_lrr <- mat_lrr[rownames(mat_lrr) %in% names(clusters), , drop = FALSE]

            clusters_mat <- clusters[rownames(mat_lrr)]

            # Compute matrix of averages per cluster
            mat_lrr_av <- t(sapply(unique(clusters_mat), function(cl) {
                colMeans(mat_lrr[clusters_mat == cl, , drop = FALSE], na.rm = TRUE)
            }))

            # Replace cell names by cluster names
            rownames(mat_lrr_av) <- unique(clusters_mat)
            mat_lrr_av <- mat_lrr_av[rownames(mat_lrr_av), , drop = FALSE]

            # Replace log ratio matrix per cell by average matrix per cluster
            x@omics[[omic]]@coverage$logratio$mat <- mat_lrr_av
        }

        # Define new clusters in correct order
        clusters <- setNames(clusters_order, clusters_order)

        # Define new "cells" becoming clusters
        common_cells <- Reduce(intersect, lapply(muscadet::matLogRatio(x), rownames))
        all_cells <- Reduce(union, lapply(muscadet::matLogRatio(x), rownames))

        common_cells <- common_cells[order(match(common_cells, clusters))]
        all_cells <- all_cells[order(match(all_cells, clusters))]
    }

    if (!quiet) {
        # Print information messages
        message("---- Heatmap Parameters ----")
        message("Omics in the muscadet object: ",
                paste(sapply(slot(x, "omics"), function(n)
                    slot(n, "label.omic")), collapse = ", "))


        if (!averages) {
            omic_dims <- sapply(slot(x, "omics"), function(m) {
                paste0(
                    m@label.omic,
                    ": ",
                    nrow(muscadet::matLogRatio(m)),
                    " cells x ",
                    ncol(muscadet::matLogRatio(m)),
                    " features"
                )
            })
            message("Omics log R ratio data dimensions:\n  ",
                    paste(omic_dims, collapse = "\n  "))

            if (!is.null(partition))
                message("Clustering partition: ", partition)
            if (!is.null(clusters))
                message("Custom clusters provided: ",
                        length(unique(clusters)),
                        " clusters.")
            if (is.null(partition) && is.null(clusters))
                message("Assigned clusters: ",
                        length(unique(x@cnacalling$clusters)),
                        " clusters.")

            message(
                "Number of cells: ",
                length(all_cells),
                " total (",
                length(common_cells),
                " common across all omics)."
            )

            message("Show missing cells: ", show_missing)

        } else {
            message("Average log ratio per cluster: ", averages)

            omic_dims <- sapply(slot(x, "omics"), function(m) {
                paste0(
                    m@label.omic,
                    ": ",
                    nrow(muscadet::matLogRatio(m)),
                    " clusters x ",
                    ncol(muscadet::matLogRatio(m)),
                    " features"
                )
            })
            message("Omics log R ratio data dimensions:\n  ",
                    paste(omic_dims, collapse = "\n  "))

        }

        message("Bulk LRR annotations: ",
                ifelse(add_bulk_lrr, slot(x, "bulk.data")[["label"]], add_bulk_lrr))
        message("White scale quantiles: ", paste(white_scale, collapse = " - "))
        message("Output file: ", ifelse(is.null(filename), "None", filename))
    }

    # Create list of heatmap objects
    list_ht <- lapply(names(x@omics), function(omic_name) {
        muscomic <- x@omics[[omic_name]]

        pdf(file = NULL)
        ComplexHeatmap::ht_opt(message = F)

        # Get chromosome info for features
        coord <- muscomic@coverage$logratio$coord.features
        chrom <- factor(coord[coord$keep, "CHROM"], levels = unique(coord[coord$keep, "CHROM"]))

        # Define color breaks for heatmap using white_scale argument
        col_breaks <- c(-5, muscomic@coverage$logratio$ref.logratio.perc[as.character(white_scale[[omic_name]])], 5)

        # Extract log R ratio matrix
        if (show_missing == TRUE) {
            mat <- matLogRatio(muscomic)
            cells.diff <- setdiff(all_cells, rownames(mat)) # identify missing cells
            if (length(cells.diff) > 0) {
                mat.na <- matrix(
                    data = NA,
                    nrow = length(cells.diff),
                    ncol = ncol(mat),
                    dimnames = list(cells.diff, colnames(mat))
                )
                mat <- rbind(mat, mat.na)[all_cells, ] # add empty rows to matrix
            } else {
                mat <- mat[all_cells, , drop = FALSE]
            }
        }
        if (show_missing == FALSE) {
            mat <- matLogRatio(muscomic)
            mat <- mat[common_cells, , drop = FALSE]
        }

        # Create heatmap
        ht <- ComplexHeatmap::Heatmap(
            mat,
            name = muscomic@label.omic,
            heatmap_legend_param = list(title = muscomic@label.omic),
            show_column_names = F,
            show_row_names = F,
            cluster_columns = F,
            column_split = chrom,
            column_title = paste(
                muscomic@label.omic,
                "coverage on",
                muscomic@coverage$logratio$label.features
            ),
            row_title_gp = grid::gpar(fontsize = 10),
            column_title_gp = grid::gpar(fontsize = 10),
            border_gp = grid::gpar(col = "black", lwd = 1),
            heatmap_height = unit(12, "cm"),
            heatmap_width = unit(18, "cm"),
            col = circlize::colorRamp2(col_breaks, c(
                "#00008E", "white", "white", "#630000"
            )),
            row_title_rot = 0,
            raster_device = "png",
            raster_quality = 3
        )

        # Add empty chromosome annotation
        emp <- ComplexHeatmap::anno_empty(border = FALSE, height = unit(2, "mm"))
        ht <- ComplexHeatmap::attach_annotation(ht,
                                                columnAnnotation(
                                                    emp_annot = emp,
                                                    name = paste0("chrom_", muscomic@label.omic)
                                                ),
                                                side = "top")
        # Rename chromosome annotation
        ht@top_annotation@anno_list[["emp_annot"]]@name <- paste0("chrom_", muscomic@label.omic)
        ht@top_annotation@anno_list[["emp_annot"]]@label <- paste0("chrom_", muscomic@label.omic)
        ht@top_annotation@anno_list[["emp_annot"]]@name_param[["label"]] <- paste0("chrom_", muscomic@label.omic)
        names(ht@top_annotation@anno_list)[1] <- paste0("chrom_", muscomic@label.omic)

        # Add bulk LRR data as annotation
        if (add_bulk_lrr) {
            # Retrieve bulk lrr values on features
            bulk_df <- getLogRatioBulk(muscomic, x@bulk.data$logratio)
            # Define color scale
            bulk_col <- list(circlize::colorRamp2(c(
                min(bulk_df$bulk.lrr, na.rm = TRUE),
                median(bulk_df$bulk.lrr, na.rm = TRUE),
                max(bulk_df$bulk.lrr, na.rm = TRUE)
            ), c("#00008E", "white", "#630000")))
            names(bulk_col) <- paste0("bulk_", muscomic@label.omic)
            # Create data frame for annotation
            bulk_df <- as.data.frame(bulk_df$bulk.lrr)
            colnames(bulk_df) <- paste0("bulk_", muscomic@label.omic)
            # Create bulk LRR annotation
            annot_bulk <- ComplexHeatmap::HeatmapAnnotation(
                df = bulk_df,
                col = bulk_col,
                annotation_label = x@bulk.data$label,
                annotation_name_side = "left",
                annotation_name_gp = grid::gpar(cex = 0.7)
            )
            # Attach annotation to heatmap
            ht <- ComplexHeatmap::attach_annotation(ht, annot_bulk, side = "bottom")
        }

        dev.off()
        return(ht)
    })


    # Create empty annotation to add number of cells (ncells)
    ha <- ComplexHeatmap::rowAnnotation(ncells = anno_empty(border = FALSE, width = unit(12, "mm")))
    # Combine heatmaps and ncells empty annotation
    ht_list <- ha + Reduce("+", list_ht)

    # Create layout for the list of heatmaps
    pdf(file = NULL)
    if (show_missing == TRUE) {
        if (is.null(clusters) & !is.null(partition)) {
            # 1. Cluster assignments from the muscadet object clustering with all cells
            clusters <- x@clustering$clusters[[as.character(partition)]]
        } else if (is.null(clusters) & is.null(partition)) {
            clusters <- x@cnacalling$clusters
        }
        n_cells <- table(clusters)[unique(clusters)]

        # Add supplementary annotations
        if (!is.null(row_annots)) {
            # Filter cells to match cells in heatmap
            row_annots <- lapply(row_annots, function(ha) {
                ha_cells <- names(ha@anno_list[["annot"]]@fun@var_env[["value"]])
                ha <- ha[which(ha_cells %in% all_cells)]
                ha
            })

            ht_list <- Reduce("+", row_annots) + ht_list

            ht_gap <- unit(c(rep(.3, length(row_annots)), .3, 1), "cm")

            annotation_legend_list <- lapply(row_annots, function(ha) {
                Legend(
                    labels = ha@anno_list[["annot"]]@fun@var_env[["color_mapping"]]@levels,
                    title = ha@anno_list[["annot"]]@label,
                    legend_gp = grid::gpar(fill = ha@anno_list[["annot"]]@fun@var_env[["color_mapping"]]@colors)
                )
            })

        } else {
            ht_gap <- unit(1, "cm")
            annotation_legend_list <- list()
        }

        # row_split: Handle the case of a single row (with averages = TRUE and a single cluster)
        if (length(all_cells) > 1) {
            row_split <- factor(clusters[all_cells], levels = unique(clusters))
        } else {
            row_split <- NULL
        }

        # Draw heatmap
        ht_all <- ComplexHeatmap::draw(
            ht_list,
            column_title = title,
            ht_gap = ht_gap,
            row_split = row_split,
            row_order = names(clusters),
            cluster_rows = F,
            annotation_legend_list = annotation_legend_list,
            merge_legend = TRUE
        )


    } else if (show_missing == FALSE) {

        # Add supplementary annotations
        if (!is.null(row_annots)) {
            # Filter cells to match cells in heatmap
            row_annots <- lapply(row_annots, function(ha) {
                ha_cells <- names(ha@anno_list[["annot"]]@fun@var_env[["value"]])
                ha <- ha[which(ha_cells %in% common_cells)]
                ha
            })

            ht_list <- Reduce("+", row_annots) + ht_list

            ht_gap <- unit(c(rep(.3, length(row_annots)), .3, 1), "cm")

            annotation_legend_list <- lapply(row_annots, function(ha) {
                Legend(
                    labels = ha@anno_list[["annot"]]@fun@var_env[["color_mapping"]]@levels,
                    title = ha@anno_list[["annot"]]@label,
                    legend_gp = grid::gpar(fill = ha@anno_list[["annot"]]@fun@var_env[["color_mapping"]]@colors)
                )
            })

        } else {
            ht_gap <- unit(1, "cm")
            annotation_legend_list <- list()
        }

        if (is.null(clusters) & !is.null(partition)) {

            if (x@clustering$params$method == "seurat") {
                # 2. Clustering seurat from the muscadet object clustering with common cells
                clusters <- x@clustering$clusters[[as.character(partition)]]
                n_cells <- table(clusters[common_cells])[unique(clusters[common_cells])]

                # row_split: Handle the case of a single row (with averages = TRUE and a single cluster)
                if (length(common_cells) > 1) {
                    row_split <- factor(clusters[common_cells], levels = unique(clusters))
                } else {
                    row_split <- NULL
                }

                # Draw heatmap
                ht_all <- ComplexHeatmap::draw(
                    ht_list,
                    column_title = title,
                    ht_gap = ht_gap,
                    row_split = row_split,
                    row_order = names(clusters)[names(clusters) %in% common_cells],
                    cluster_rows = FALSE,
                    annotation_legend_list = annotation_legend_list,
                    merge_legend = TRUE
                )

            } else if (x@clustering$params$method == "hclust") {
                # 2. Clustering hclust from the muscadet object clustering with common cells
                hc <- x@clustering$hclust # hclust object to print the dendrogram on the heatmap
                clusters <- x@clustering$clusters[[as.character(partition)]]
                n_cells <- table(clusters[common_cells])[unique(clusters[common_cells])]

                # Draw heatmap
                ht_all <- ComplexHeatmap::draw(
                    ht_list,
                    column_title = title,
                    ht_gap = ht_gap,
                    cluster_rows = hc,
                    row_split = as.integer(partition),
                    row_dend_reorder = FALSE,
                    annotation_legend_list = annotation_legend_list,
                    merge_legend = TRUE
                )
            }

        } else {
            # 3. Custom/Assigned cluster assignments vector
            if (is.null(clusters) & is.null(partition)) {
                clusters <- x@cnacalling$clusters
            }
            n_cells <- table(clusters[common_cells])

            # row_split: Handle the case of a single row (with averages = TRUE and a single cluster)
            if (length(common_cells) > 1) {
                row_split <- factor(clusters[common_cells], levels = unique(clusters))
            } else {
                row_split <- NULL
            }

            # Draw heatmap
            ht_all <- ComplexHeatmap::draw(
                ht_list,
                column_title = title,
                ht_gap = ht_gap,
                row_split = row_split,
                row_order = names(clusters)[names(clusters) %in% common_cells],
                cluster_rows = FALSE,
                annotation_legend_list = annotation_legend_list,
                merge_legend = TRUE
            )
        }
    }
    dev.off()

    # Save complete plot of heatmaps as PNG PDF or SVG
    if (!is.null(filename)) {
        if (grepl(".png", basename(filename))) {
            grDevices::png(
                filename = filename,
                width = ht_all@ht_list_param[["width"]],
                height = ht_all@ht_list_param[["height"]],
                units = "mm",
                res = png_res
            )
        } else if (grepl(".pdf", basename(filename))) {
            grDevices::pdf(
                file = filename,
                width = (ht_all@ht_list_param[["width"]]) / 25.4,
                # from mm to inches
                height = (ht_all@ht_list_param[["height"]]) / 25.4 # from mm to inches
            )
        } else if (grepl(".svg", basename(filename))) {
            grDevices::svg(
                filename = filename,
                width = (ht_all@ht_list_param[["width"]]) / 25.4,
                # from mm to inches
                height = (ht_all@ht_list_param[["height"]]) / 25.4 # from mm to inches
            )
        }
    } else {
        pdf(file = NULL)
    }

    # Print plot
    print(ht_all)

    # Add annotation: number of cells per cluster
    for (i in 1:length(n_cells)) {
        ComplexHeatmap::decorate_annotation("ncells", slice = i, envir = environment(), {
            grid::grid.rect(
                x = 1,
                width = unit(2, "mm"),
                gp = grid::gpar(fill = colors[i], col = NA),
                just = "right"
            )
            if (!averages) {
                grid::grid.text(
                    n_cells[i],
                    x = 0.7,
                    just = "right",
                    gp = grid::gpar(cex = 0.75)
                )
            }
        })
    }

    # Add annotation: chromosome numbers
    for (ht_name in unlist(lapply(list_ht, function(l) l@name))) {
        chrom_names <- unique(names(ht_all@ht_list[[ht_name]]@column_order_list))
        for (i in 1:length(chrom_names)) {
            ComplexHeatmap::decorate_annotation(paste0("chrom_", ht_name),
                                                slice = i,
                                                envir = environment(),
                                                {
                                                    grid::grid.text(chrom_names[i],
                                                                    just = "center",
                                                                    gp = grid::gpar(fontsize = 8.5))
                                                })
        }
    }

    # store plot object
    if (is.null(filename)) {
        plot.obj <- list(
            plot = grid::grid.grab(),
            width = ht_all@ht_list_param[["width"]],
            height = ht_all@ht_list_param[["height"]]
        )
    }

    dev.off()

    # return the plot object
    if (is.null(filename)) {
        return(plot.obj)
    }
}





#' Silhouette plot for `muscadet` object
#'
#' Generate a silhouette plot for a specified clustering partition within a
#' [`muscadet`][muscadet-class] object.
#'
#' @param x A [`muscadet`][muscadet-class] object containing clustering data (using
#'   [clusterMuscadet()]).
#'
#' @param partition Value specifying the clustering partition to plot (`numeric`
#'   or `character`). It should be either the resolution or the k number of
#'   cluster (k) used for clustering depending on the clustering method
#'   (`res_range` or `k_range` with [clusterMuscadet()]).
#'
#' @param colors Vector of colors for the cluster annotation (`character`
#'   vector). Default is `NULL`, which uses predefined colors.
#'
#' @param title Character string for the title of the plot (`character`
#'   string). If `NULL`, a default title is generated.
#'
#' @param annotations `TRUE` or `FALSE` (`logical`). Whether to add annotations
#'   per clusters. By default: `TRUE`.
#'
#' @return A ggplot object representing the silhouette plot.
#'
#' @include objects.R
#'
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' library("ggplot2")
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#' plotSil(exdata_muscadet, partition = 0.3)
#'
#' # Loop over partitions
#' for (p in names(exdata_muscadet$clustering$clusters)) {
#'     plot <- plotSil(exdata_muscadet, p)
#'     ggsave(paste0("plot_silhouette_", p, ".png"), plot)
#' }
#' }
#'
plotSil <- function(x,
                    partition,
                    colors = NULL,
                    title = NULL,
                    annotations = TRUE) {
    # Validate input: x must be a muscadet object
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Validate the muscadet object contains clustering results
    stopifnot(
        "The muscadet object `x` does not contain clustering data (use `clusterMuscadet()` to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Validate the clustering result for the specified partition
    stopifnot(
        "The muscadet object `x` must contain clustering results for the specified `partition`." =
            as.character(partition) %in% names(x@clustering$clusters)
    )

    # Validate partition with at least 2 clusters
    stopifnot(
        "The selected clustering `partition` contains only one cluster, silhouette scores cannot be computed." =
            length(unique(x@clustering$clusters[[as.character(partition)]])) != 1
    )


    # Set default color palette for clusters if not provided
    if (is.null(colors)) {
        colors <- rep(
            c("#FABC2A", "#7FC97F", "#EE6C4D", "#39ADBD", "#BEAED4", "#FEE672", "#F76F8E",
              "#487BEA", "#B67BE6", "#F38D68", "#7FD8BE", "#F2AFC0"), 3)
    }

    # Generate a default title if none is provided
    if (is.null(title)) {

        if (x$clustering$params$method == "seurat") partition_name <- "res"
        if (x$clustering$params$method == "hclust") partition_name <- "k"

        title <- paste0("Silhouette Plot (", partition_name, " = ", as.character(partition), ")")
    }

    # Extract silhouette data for the specified partition
    sil <- x@clustering$silhouette$sil.obj[[as.character(partition)]]
    df <- as.data.frame(sil[, 1:3], stringsAsFactors = TRUE)

    # Order data for plotting
    df <- df[order(df$cluster, -df$sil_width), ]
    df$name <- factor(rownames(df), levels = rownames(df))
    df$cluster <- as.factor(df$cluster)

    # Calculate average silhouette width per cluster
    avg_sil_width <- stats::aggregate(df$sil_width ~ df$cluster, data = df, mean)
    colnames(avg_sil_width) <- c("cluster", "sil_width")

    # Determine positions for cluster annotations
    cluster_counts <- table(df$cluster)
    cluster_positions <- sapply(unique(df$cluster), function(i) {
        # Reverse y-axis positions since it is flipped
        rev_y_position <- sum(cluster_counts) -
            (sum(cluster_counts[as.numeric(levels(df$cluster)[1:i])]) - cluster_counts[i] / 2)
        rev_y_position
    })

    # Create the ggplot object
    p <- ggplot(df, aes(
        x = .data$sil_width,
        y = .data$name,
        fill = .data$cluster
    )) +
        geom_bar(stat = "identity", width = 1) +
        scale_y_discrete(limits = rev) +
        theme_bw() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
            legend.position = "none"
        ) +
        scale_fill_manual(values = colors) +
        labs(
            x = "Silhouette Widths",
            y = "Cells",
            title = title,
            subtitle = paste0("Average Silhouette Width = ", round(mean(df$sil_width), 4))
        ) +
        ggplot2::xlim(c(min(0, min(df$sil_width) - 0.01), max(df$sil_width) + 0.25)) +
        geom_vline(
            xintercept = mean(df$sil_width),
            linetype = "dashed",
            color = "red"
        )

    if (annotations) {
        p <- p + labs(
            subtitle = paste0(
                "Average Silhouette Width = ",
                round(mean(df$sil_width), 4),
                "\n",
                "Annotation: cluster | number of cells | average silhouette width"
            )
        )

        # Add annotations for each cluster
        for (i in seq_along(cluster_positions)) {
            cluster <- unique(df$cluster)[i]
            avg_width <- round(avg_sil_width$sil_width[avg_sil_width$cluster == cluster], 4)
            count <- cluster_counts[i]
            p <- p + geom_text(
                label = paste(cluster, "|", count, "|", avg_width),
                x = max(df$sil_width) + 0.05,
                y = cluster_positions[i],
                inherit.aes = FALSE,
                hjust = 0,
                vjust = -0.5,
                check_overlap = TRUE
            )
        }
    }

    return(p)
}



#' Plot clustering validation indexes for a `muscadet` object
#'
#' It generates a plot of clustering validation indexes for a clustering
#' partition within a [`muscadet`][muscadet-class] object. The index values are
#' computed only using distances between common cells across omics in the
#' [`muscadet`][muscadet-class] object.
#'
#' @param x A [`muscadet`][muscadet-class] object containing clustering data
#'   (generated using [clusterMuscadet()]).
#'
#' @param index Character vector specifying one or more validation indexes to
#'   plot among `"silhouette"`, `"dunn2"`, `"daviesbouldin"`, `"pearsongamma"`,
#'   and `"c"`. If `NULL`, by default all available indexes are included. If
#'   multiple indexes are selected, the values are normalized for comparability.
#'
#' @param colors Vector of colors for each index in the plot (`character`
#'   vector). Default is `NULL`, which uses predefined colors for the indexes.
#'
#' @param title Character string for the title of the plot (`character` string).
#'   If `NULL`, a default title is generated.
#'
#' @return A ggplot object visualizing the clustering validation indexes across
#'   different clustering partitions (`res` resolution or `k` number of clusters
#'   depending on the used clustering method).
#'
#' @details The function computes several clustering validation indexes,
#'   including:
#'   \itemize{
#'     \item \strong{Silhouette}: Measures how similar an object is to its own
#'     cluster compared to others (see [cluster::silhouette()]). Average of
#'     individual silhouette widths.
#'     \item \strong{Dunn2}: The ratio of the smallest distance between
#'     observations in different clusters to the largest within-cluster distance
#'     (see [fpc::cluster.stats()] `$dunn2`). Minimum average dissimilarity
#'     between two cluster / maximum average within cluster dissimilarity.
#'     \item \strong{Davies-Bouldin}: Measures cluster compactness and
#'     separation (see [clusterSim::index.DB()]).
#'     \item \strong{Pearson's Gamma}: Evaluates the goodness of clustering
#'     based on correlation (see [fpc::cluster.stats()] `$pearsongamma`).
#'     Correlation between distances and a 0-1-vector where 0 means same
#'     cluster, 1 means different clusters. "Normalized gamma" in Halkidi et al.
#'     (2001).
#'     \item \strong{C Index} (Hubert & Levin C index): Measures the internal
#'     cluster quality compared to random data (see [clusterSim::index.C()]).
#'   }
#'
#'   If multiple indexes are selected, the values are normalized to fall between
#'   0 and 1. For indexes that are better when minimized ("pearsongamma" and
#'   "c"), their values are reversed for easier comparison. The partition for
#'   which the mean of indexes is maximal is highlighted with a dot.
#'
#' @include objects.R
#'
#' @import ggplot2
#' @importFrom stats as.dist
#' @importFrom cluster silhouette
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' # Plot all indexes
#' plotIndexes(exdata_muscadet)
#' ggsave("plot_indexes.png", width = 8, height = 4)
#'
#' # Plot a specific index
#' plotIndexes(exdata_muscadet, index = "silhouette")
#' ggsave("plot_indexes_sil.png", width = 7, height = 4)
#' }
#'
plotIndexes <- function(x,
                        index = "silhouette",
                        colors = NULL,
                        title = NULL) {
    # Validate input: x must be a muscadet object
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Validate the muscadet object contains clustering results
    stopifnot(
        "The muscadet object `x` does not contain clustering data (use `clusterMuscadet()` to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Define default indexes if none are provided
    if (is.null(index)) {
        index <- c("silhouette",
                   "dunn2",
                   "daviesbouldin",
                   "pearsongamma",
                   "c")
    }
    # Check for indexes correct names
    index <- match.arg(
        index,
        c("silhouette", "dunn2", "daviesbouldin", "pearsongamma", "c"),
        several.ok = TRUE
    )

    # Set default colors if not provided
    if (is.null(colors)) {
        colors <- c(
            "silhouette" = "brown",
            "dunn2" = "coral2",
            "daviesbouldin" = "tan2",
            "pearsongamma" = "turquoise4",
            "c" = "skyblue4"
        )
    }

    # Extract the clustering data
    partitions <- as.numeric(names(x@clustering$clusters))  # Partitions (res or k)
    dist <- stats::as.dist(x@clustering$dist)  # Distance matrix to dist object

    # Remove partitions with only one cluster
    n_clusters <- sapply(partitions, function(p) length(unique(x@clustering$clusters[[as.character(p)]])))
    partitions <- partitions[n_clusters > 1]

    stopifnot("All clustering partitions contain only one cluster, indexes cannot be computed." =
                  length(partitions) != 0)

    # Initialize a data frame to store index values
    df_indexes <- data.frame(partition = as.factor(partitions))

    # Compute selected indexes for each partition
    for (p in partitions) {
        # Extract cluster assignments for the current partition
        clusters <- x@clustering$clusters[[as.character(p)]]
        clusters <- clusters[labels(dist)] # Restrict to common cells in dist

        # Compute each index if selected
        if ("silhouette" %in% index) {
            df_indexes[df_indexes$partition == p, "silhouette"] <- summary(cluster::silhouette(as.integer(clusters), dist))[["avg.width"]]
        }
        if ("dunn2" %in% index) {
            if (!requireNamespace("fpc", quietly = TRUE)) {
                stop("Package 'fpc' is required to compute 'dunn2' index. Install it with install.packages('fpc').")
            }
            df_indexes[df_indexes$partition == p, "dunn2"] <- fpc::cluster.stats(dist, as.integer(clusters))$dunn2
        }
        if ("daviesbouldin" %in% index) {
            if (!requireNamespace("clusterSim", quietly = TRUE)) {
                stop("Package 'clusterSim' is required to compute 'daviesbouldin' index. Install it with install.packages('clusterSim').")
            }
            df_indexes[df_indexes$partition == p, "daviesbouldin"] <- clusterSim::index.DB(dist, as.integer(clusters))$DB
        }
        if ("pearsongamma" %in% index) {
            if (!requireNamespace("fpc", quietly = TRUE)) {
                stop("Package 'fpc' is required to compute 'pearsongamma' index. Install it with install.packages('fpc').")
            }
            df_indexes[df_indexes$partition == p, "pearsongamma"] <- fpc::cluster.stats(dist, as.integer(clusters))$pearsongamma
        }
        if ("c" %in% index) {
            if (!requireNamespace("clusterSim", quietly = TRUE)) {
                stop("Package 'clusterSim' is required to compute 'c' index. Install it with install.packages('clusterSim').")
            }
            df_indexes[df_indexes$partition == p, "c"] <- clusterSim::index.C(dist, as.integer(clusters))
        }
    }

    # If multiple indexes selected:
    if (length(index) > 1) {
        # Normalize indexes to values between 0 and 1
        df_indexes[, 2:ncol(df_indexes)] <- apply(df_indexes[, 2:ncol(df_indexes)], 2, function(x) {
            (x - min(x)) / diff(range(x))
        })
        # Reverse indexes that should be minimized, for comparability
        if ("pearsongamma" %in% index) {
            df_indexes <- mutate(df_indexes, pearsongamma = 1 - .data$pearsongamma)
        }
        if ("c" %in% index) {
            df_indexes <- mutate(df_indexes, c = 1 - .data$c)
        }
    }

    # Find optimal partition based on selected indexes
    if (ncol(df_indexes) > 2) {
        # Find the partition with the maximum mean index value
        opt_p <- df_indexes[which(rowMeans(df_indexes[, 2:ncol(df_indexes)]) == max(rowMeans(df_indexes[, 2:ncol(df_indexes)]))), "partition"]
    } else {
        if (colnames(df_indexes)[2] %in% c("pearsongamma", "c")) {
            # Find partition with the minimal value for indexes to minimize
            opt_p <- df_indexes[which(df_indexes[, 2] == min(df_indexes[, 2])), "partition"]
        } else {
            # Find partition with the maximal value for indexes to maximize
            opt_p <- df_indexes[which(df_indexes[, 2] == max(df_indexes[, 2])), "partition"]
        }
    }

    # Transform the data frame for plotting
    df_plot <- stats::reshape(
        df_indexes,
        varying = setdiff(names(df_indexes), "partition"),
        v.names = "Value",
        timevar = "Index",
        times = setdiff(names(df_indexes), "partition"),
        direction = "long"
    )
    row.names(df_plot) <- NULL
    df_plot$Index <- factor(df_plot$Index, levels = index) # Maintain order of indexes

    # Labels
    if (x@clustering$params$method == "seurat") x_lab <- "Clustering partition (resolution)"
    if (x@clustering$params$method == "hclust") x_lab <- "Clustering partition (number of clusters)"
    if (length(index) > 1) y_lab <- "Normalized Index Value"
    if (length(index) == 1) {
        if (index == "silhouette") y_lab <- "Average silhouette width"
        if (index == "dunn2") y_lab <- "Dunn2 index"
        if (index == "pearsongamma") y_lab <- "Pearson's Gamma index"
        if (index == "daviesbouldin") y_lab <- "Davies-Bouldin index"
        if (index == "c") y_lab <- "Hubert & Levin C index"
    }
    # Generate a default title if none is provided
    if (is.null(title) & length(index) > 1) {
        title <- "Clustering Validation Indexes Across Clustering Partitions"
    }
    if (is.null(title) & length(index) == 1) {
        if (index == "silhouette") title <- "Silhouette Score Across Clustering Partitions"
        if (index == "dunn2") y_lab <- "Dunn2 Index Across Clustering Partitions"
        if (index == "pearsongamma") y_lab <- "Pearson's Gamma Index Across Clustering Partitions"
        if (index == "daviesbouldin") y_lab <- "Davies-Bouldin Index Across Clustering Partitions"
        if (index == "c") y_lab <- "C Index Across Clustering Partitions"
    }

    # Generate the plot
    plot <- ggplot(df_plot,
                   aes(
                       x = .data$partition,
                       y = .data$Value,
                       color = .data$Index,
                       group = .data$Index
                   )) +
        geom_line(linewidth = 1) +
        geom_point(data = df_plot[which(df_plot$partition == opt_p), ]) +
        labs(x = x_lab, y = y_lab, title = title) +
        scale_x_discrete(breaks = unique(df_plot$partition)) +
        scale_color_manual(
            values = colors,
            labels = c(
                "silhouette" = "Silhouette",
                "dunn2" = "Dunn2",
                "pearsongamma" = "Pearson's Gamma",
                "daviesbouldin" = "Davies-Bouldin",
                "c" = "C Index"
            )
        ) +
        theme_classic() +
        theme(legend.position = "none")

    if (length(index) > 1) {
        plot <- plot +
            theme(legend.position = "right")
    }

    return(plot)
}



#' Plot CNA profiles from muscadet object
#'
#' This function generates a multi-panel plot of copy number alteration (CNA)
#' profiles from a [`muscadet`][muscadet-class] object, including: log R ratios
#' values, log odds ratio (or variant allele frequency), copy numbers, CNA
#' status and cell fractions.
#'
#' @param x A [`muscadet`][muscadet-class] object containing CNA calling data to
#'   be visualized (generated using [cnaCalling()]).
#' @param data Either a cluster identifier to plot data of a cluster or
#'   "allcells" to plot data on all cells.
#' @param title An optional title for the plot. Default is `NULL`.
#' @param allelic.type A character string indicating the allelic metric to plot:
#'   "lor" for log odds ratio or "vaf" for variant allele frequency. Default is
#'   "vaf".
#' @param point.cex Numeric vector of length 1 or 2 specifying the size of
#'   points in the plots. If a single value is provided, it will be replicated
#'   for both plots. Default is `c(0.4, 0.5)`.
#' @param chrom.colors A character vector of length 2 defining alternating
#'   chromosome colors. Default is `c("slategrey", "skyblue")`.
#' @param var.colors A character vector of length 2 for variant positions point
#'   colors depending of variant allele frequency in all cells. Use "none" to
#'   use the alternating chromosome colors (defined by `chrom.colors`). Default
#'   is `c("peachpuff2", "paleturquoise3")`.
#' @param cn.colors A character vector of length 2 for total copy number and
#'   minor allele copy number segment colors. Default is `c("black", "brown2")`.
#' @param cna.colors A vector of 3 colors for CNA states: gain, loss, and cnloh
#'   (or named vector where names are "gain", "loss", and "cnloh" and the values
#'   are their respective colors). Default is `c("gain" = "#EF6F6AFF", "loss" =
#'   "#6699CCFF", "cnloh" = "#44AA99FF")`.
#' @param cf.colors A character vector of length 3 for cellular fraction
#'   gradient (of 10 values): start color of the gradient, end color of the
#'   gradient, and color for normal diploid (depending on the ploidy). Default
#'   is `c("white", "steelblue", "bisque2")`.
#' @param dipLogR.color A character string for the diploid log R ratio line
#'   color. Default is "magenta4".
#' @param seg.color A character string for the color of segment medians. Default
#'   is "brown2".
#'
#' @return A multi-panel plot of CNA profiles is produced.
#'
#' @include objects.R
#'
#' @import graphics
#' @import grDevices
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' # Plot profile for first cluster
#' pdf("CNAprofile_allcells.pdf", width = 15, height = 7.5) # Save as PDF
#' plotProfile(exdata_muscadet, data = "1", title = "Example dataset - cluster 1")
#' dev.off()
#'
#' # Plot profile for all cells
#' pdf("CNAprofile_allcells.pdf", width = 15, height = 7.5) # Save as PDF
#' plotProfile(exdata_muscadet, data = "allcells", title = "Example dataset - all cells")
#' dev.off()
#' }
#'
plotProfile <- function(x,
                        data,
                        title = NULL,
                        allelic.type = "vaf",
                        point.cex = c(0.4, 0.5),
                        chrom.colors = c("slategrey", "skyblue"),
                        var.colors = c("peachpuff2", "paleturquoise3"),
                        cn.colors = c("grey20", "brown2"),
                        cna.colors = c(
                            "gain" = "#EF6F6AFF",
                            "loss" = "#6699CCFF",
                            "cnloh" = "#44AA99FF"
                        ),
                        cf.colors = c("white", "grey20", "bisque2"),
                        dipLogR.color = c("magenta4"),
                        seg.color = c("brown2")) {
    # Argument checks
    stopifnot(
        "Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"),

        "Input object `x` must contain CNA calling data." = !is.null(x@cnacalling$table),

        "Invalid `allelic.type'. Use allelic.type = \"lor\" or allelic.type = \"vaf\"." = allelic.type %in% c("lor", "vaf"),
        "`point.cex` must be numeric." = is.numeric(point.cex),
        "`chrom.colors` must be a character vector of length 2." = is.character(chrom.colors) &&
            length(chrom.colors) == 2,
        "`var.colors` must be a character vector of length 2 or `none`." = (is.character(var.colors) &&
                                                                                length(var.colors) == 2) ||
            all(var.colors == "none"),
        "`cn.colors` must be a character vector of length 2." = is.character(cn.colors) &&
            length(cn.colors) == 2,
        "`cf.colors `must be a character vector of length 3." = is.character(cf.colors) &&
            length(cf.colors) == 3,
        "`dipLogR.color` must be a single character value." = is.character(dipLogR.color) &&
            length(dipLogR.color) == 1,
        "`seg.color` must be a single character value." = is.character(seg.color) &&
            length(seg.color) == 1
    )
    if (length(point.cex) == 1) {
        point.cex <- rep(point.cex, 2)
    }


    # Extract data -------------------------------------------------------------
    if (data == "allcells") {
        pos <- x@cnacalling$positions.allcells
        segs <- x@cnacalling$segments.allcells
        ploidy <- x@cnacalling$ploidy
        dipLogR <- x@cnacalling$dipLogR.allcells
    } else {
        pos <- x@cnacalling$positions
        segs <- x@cnacalling$segments
        ploidy <- x@cnacalling$ploidy
        dipLogR <- x@cnacalling$dipLogR.clusters
    }
    stopifnot(
        "`data` must be \"allcells\" or a valid cluster identifier." =
            (data == "allcells" || data %in% unique(pos$cluster))
    )
    if (data != "allcells") {
        pos <- pos[which(pos$cluster == data), ]
        segs <- segs[which(segs$cluster == data), ]
    }

    # Adjust chromosomes levels to get only numeric chromosomes
    pos$chrom <- factor(pos$chrom, levels = unique(pos$chrom))
    segs$chrom <- factor(segs$chrom, levels = unique(segs$chrom))
    chromlevels <- levels(pos$chrom)
    levels(pos$chrom) <- 1:length(levels(pos$chrom))
    levels(segs$chrom) <- 1:length(levels(segs$chrom))
    pos$chrom <- as.numeric(pos$chrom)
    segs$chrom <- as.numeric(segs$chrom)

    segs <- dplyr::mutate(
        segs,
        gnl = dplyr::case_when(
            round(.data$tcn.em) - ploidy > 0 ~ 1,
            round(.data$tcn.em) - ploidy < 0 ~ -1,
            round(.data$tcn.em) - ploidy == 0 ~ 0
        ),
        loh = dplyr::case_when(round(.data$lcn.em) == 0 ~ T, round(.data$lcn.em) > 0 ~ F),
        state = dplyr::case_when(
            .data$gnl == 1 ~ "gain",
            .data$gnl == -1 ~ "loss",
            .data$gnl == 0 & .data$loh == TRUE ~ "cnloh",
            .data$gnl == 0 & .data$loh == FALSE ~ "neu"
        ),
        cna = dplyr::case_when(.data$gnl == 0 & .data$loh == FALSE ~ F, .data$gnl != 0 | .data$loh == TRUE ~ T),
        cna_state = dplyr::case_when(.data$cna == T ~ .data$state, .data$cna == F ~ NA)
    )

    # Chromosome alternating dual colors
    chrcol <- 1 + rep(segs$chrom - 2 * floor(segs$chrom / 2), segs$num.mark)

    # Chromosomes position boundaries
    chrbdry <- which(diff(pos$chrom) != 0)

    # Segment position boundaries
    segbdry <- cumsum(c(0, segs$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]


    # Layout params ------------------------------------------------------------
    def.par <- par(no.readonly = TRUE)
    layout(matrix(rep(1:5, c(9, 9, 6, 1, 1)), ncol = 1))
    par(
        mar = c(0.25, 3, 0.25, 1),
        mgp = c(1.75, 0.6, 0),
        oma = c(3, 1.5, 1.5, 0),
        xaxs = "i"
    )

    # 1- Plot the LRR data -----------------------------------------------------
    plot(
        pos$cnlr,
        pch = 16,
        cex = point.cex[1],
        col = chrom.colors[chrcol],
        ylab = "Log R ratio",
        cex.lab = 1.5,
        xaxt = "n",
        yaxt = "n"
    )
    # Add y axis labels with rotation for horizontal labels
    axis(2, las = 2)
    # Add chromosomes boundaries
    abline(v = chrbdry, lwd = 0.25)
    # # Add LRR median
    # abline(h=median(pos$cnlr, na.rm=TRUE), col="green3")
    # Add diploid LRR
    abline(h = dipLogR, col = dipLogR.color)
    # Add LRR segment medians
    segments(
        segstart,
        segs$cnlr.median,
        segend,
        segs$cnlr.median,
        lwd = 1.75,
        col = seg.color
    )

    # 2- Plot the LOR data -----------------------------------------------------
    if (any(var.colors == "none")) {
        cols <- chrom.colors[chrcol]
    } else {
        cols <- var.colors[pos$colVAR]
    }
    if (allelic.type == "lor") {
        plot(
            pos$valor,
            pch = 16,
            cex = point.cex[2],
            col = cols,
            ylab = "Log odds ratio",
            cex.lab = 1.5,
            ylim = c(-4, 4),
            xaxt = "n",
            yaxt = "n"
        )
        # Add y axis labels with rotation for horizontal labels
        axis(2, las = 2)
        # Add chromosomes boundaries
        abline(v = chrbdry, lwd = 0.25)
        # Add segments medians
        segments(segstart,
                 sqrt(abs(segs$mafR)),
                 segend,
                 sqrt(abs(segs$mafR)),
                 lwd = 1.75,
                 col = seg.color)
        segments(segstart,-sqrt(abs(segs$mafR)),
                 segend,-sqrt(abs(segs$mafR)),
                 lwd = 1.75,
                 col = seg.color)
    }
    if (allelic.type == "vaf") {
        pos[which(pos$signal == "coverage"), "vafT"] <- NA
        plot(
            pos$vafT,
            pch = 16,
            cex = point.cex[2],
            col = cols,
            ylab = "Variant allele frequency",
            cex.lab = 1.5,
            ylim = c(0, 1),
            xaxt = "n",
            yaxt = "n"
        )
        # Add y axis labels with rotation for horizontal labels
        axis(2, las = 2)
        # Add chromosomes boundaries
        abline(v = chrbdry, lwd = 0.25)
        # Add segments medians
        segments(
            segstart,
            segs$vafT.median,
            segend,
            segs$vafT.median,
            lwd = 1.75,
            col = seg.color
        )
        segments(
            segstart,
            1 - segs$vafT.median,
            segend,
            1 - segs$vafT.median,
            lwd = 1.75,
            col = seg.color
        )
    }

    # 3- Plot the estimated copy numbers and cf --------------------------------

    # Transform to tcn to log scale over 10
    tcn.i <- which(segs$tcn.em > 10 & !is.na(segs$tcn.em))
    if (length(tcn.i) > 0) {
        segs$tcn.em[tcn.i] <- 9 + log10(segs$tcn.em[tcn.i])
    }

    # Transform to lcn to log scale over 5
    lcn.i <- which(segs$lcn.em > 5)
    if (length(lcn.i) > 0) {
        segs$lcn.em[lcn.i] <- 5 + log10(segs$lcn.em[lcn.i])
    }

    plot(
        c(0, nrow(pos)),
        c(0, max(segs$tcn.em, na.rm = T)),
        type = "n",
        ylab = "Copy number",
        cex.lab = 1.5,
        xaxt = "n",
        yaxt = "n"
    )
    # Add y axis labels with rotation for horizontal labels
    axis(2, las = 2)
    # Add chromosomes boundaries
    abline(v = chrbdry, lwd = 0.25)
    # Add lcn
    segments(segstart,
             segs$lcn.em,
             segend,
             segs$lcn.em,
             lwd = 1.75,
             col = cn.colors[2])
    # Add tcn
    segments(segstart,
             segs$tcn.em,
             segend,
             segs$tcn.em,
             lwd = 1.75,
             col = cn.colors[1])

    # Add gnl status
    plot(
        c(0, nrow(pos)),
        0:1,
        type = "n",
        ylab = "",
        xaxt = "n",
        yaxt = "n"
    )
    mtext(
        "status",
        side = 2,
        at = 0.5,
        line = 1,
        las = 2,
        cex = 0.8
    )
    cnacol <- cna.colors[segs$cna_state]
    rect(segstart, 0, segend, 1, col = cnacol, border = NA)

    # Add cf
    plot(
        c(0, nrow(pos)),
        0:1,
        type = "n",
        ylab = "",
        xaxt = "n",
        yaxt = "n"
    )
    mtext(
        "cf",
        side = 2,
        at = 0.5,
        line = 1,
        las = 2,
        cex = 0.8
    )

    cfpalette <- colorRampPalette(c(cf.colors[1], cf.colors[2]))(10)
    cfcol <- cfpalette[round(10 * segs$cf.em)]
    # cfcol[segs$tcn.em == ploidy &
    #           segs$lcn.em == ploidy / 2] <- cf.colors[3]
    rect(segstart, 0, segend, 1, col = cfcol, border = NA)

    # Add chromosome ticks on x-axis -------------------------------------------

    # Number positions per chromosomes
    nn <- cumsum(table(pos$chrom))
    # Ticks
    axis(
        labels = chromlevels,
        side = 1,
        at = (nn + c(0, nn[-length(nn)])) / 2,
        cex = 1.5
    )
    # Labels
    mtext(side = 1,
          line = 1.75,
          "Chromosomes",
          cex = 1)

    # Add title ----------------------------------------------------------------
    mtext(
        title,
        side = 3,
        line = 0,
        outer = TRUE,
        cex = 1.1
    )

    # Reset layout
    par(def.par)
}




#' Plot CNA segments across clusters from a muscadet object
#'
#' This function visualizes copy number alteration (CNA) segments across
#' clusters based on data stored in a [`muscadet`][muscadet-class] object. It
#' displays CNAs for each clusters and scales the y-axis based on the proportion
#' of cells in each cluster.
#'
#' @param x A [`muscadet`][muscadet-class] object containing CNA calling data to be
#'   visualized (generated using [cnaCalling()]).
#' @param title An optional title for the plot. Default is `NULL`.
#' @param cna.colors A vector of 3 colors for CNA states: gain, loss, and cnloh
#'   (or named vector where names are "gain", "loss", and "cnloh" and the values
#'   are their respective colors). Default is `c("gain" = "#EF6F6AFF", "loss" =
#'   "#6699CCFF", "cnloh" = "#44AA99FF")`.
#' @param cf.gradient Logical. If `TRUE` adds a alpha transparency gradient on
#'   CNA color block based on the cell fraction. Default is `TRUE`.
#'
#' @return A ggplot object representing the CNA segments plot.
#'
#' @include objects.R
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom stats complete.cases
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("ggplot2")
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' # Plot CNA segments
#' p <- plotCNA(exdata_muscadet, title = "Copy Number Alterations in Example Data")
#' p
#' ggsave(
#'     filename = file.path("CNAplot.png"),
#'     plot = p,
#'     width = 3000,
#'     height = 800,
#'     units = "px"
#' )
#' }
#'
plotCNA <- function(x,
                    title = NULL,
                    cna.colors = c(
                        "gain" = "#EF6F6AFF",
                        "loss" = "#6699CCFF",
                        "cnloh" = "#44AA99FF"
                    ),
                    cf.gradient = TRUE
                    ) {
    # Argument checks
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Extract data table
    data <- x@cnacalling$table

    # Extract number of cells per cluster
    ncells <- x@cnacalling$ncells

    # Keep only autosomes to display
    data <- data[!data$chrom %in% c("X", "Y", "M"), ]

    # Check for colors
    cna_states <- c("gain", "loss", "cnloh")
    if(!is.null(names(cna.colors))) {
        if(!all(cna_states %in% names(cna.colors))) {
            cna.colors <- cna.colors[1:3]
            names(cna.colors) <- cna_states
        }
    }
    if (is.null(names(cna.colors)) & length(cna.colors) <= 3) {
        names(cna.colors) <- cna_states[1:length(cna.colors)]
    } else if (length(cna.colors) > 3) {
        cna.colors <- cna.colors[1:3]
        if (is.null(names(cna.colors))) {
            names(cna.colors) <- cna_states
        }
    }

    # Extract genome
    if (x@genome == "hg38") {
        genome_chrom <- hg38_chrom
    }
    if (x@genome == "hg19") {
        genome_chrom <- hg19_chrom
    }
    if (x@genome == "mm10") {
        genome_chrom <- mm10_chrom
    }
    chromSizes <- as.data.frame(genome_chrom)

    # Keep only autosomes to display
    chromSizes$seqnames <- as.character(chromSizes$seqnames)
    chromSizes <- chromSizes[!chromSizes$seqnames %in% c("X", "Y", "M"), ]

    # Chromosomes start coordinates
    chromStarts <- data.frame(
        chrom = factor(chromSizes$seqnames, levels = unique(chromSizes$seqnames)),
        chrom.start = c(0, cumsum(as.numeric(
            chromSizes$width
        )))[-(nrow(chromSizes) + 1)]
    )

    # Ordered chromosomes names
    chromNames <- factor(unique(chromSizes$seqnames),
                         levels = unique(chromSizes$seqnames))


    # Add chromosomes start and end coordinates on x axis
    df <- dplyr::left_join(data, chromStarts, by = "chrom") %>%
        dplyr::mutate(
            start.x = .data$start + .data$chrom.start,
            end.x = .data$end + .data$chrom.start
        )

    # Remove consensus segs that have no data in a cluster
    df <- df[complete.cases(df$cluster), ]

    # Compute proportions of cells per cluster for y axis
    prop_clus <- unique(df$prop.cluster[!is.na(df$prop.cluster)])
    prop_starts <- c(0, cumsum(prop_clus[-length(prop_clus)]))
    prop_ends <- cumsum(prop_clus)
    names(prop_starts) <- unique(df$cluster[!is.na(df$cluster)])
    names(prop_ends) <- unique(df$cluster[!is.na(df$cluster)])
    props <- unique(c(prop_starts, prop_ends))
    props_breaks <- props[-length(props)] + (diff(props) / 2)

    # Add clusters coordinates on y axis
    df <- df %>%
        dplyr::mutate(start.y = prop_starts[as.character(.data$cluster)], end.y = prop_ends[as.character(.data$cluster)])

    # Ordered levels for CNA states
    df$cna_state <- factor(df$cna_state, levels = cna_states[cna_states %in% unique(na.omit(data$cna_state))])

    # Construct plot
    cna_plot <- ggplot2::ggplot(
        df,
        aes(
            xmin = .data$start.x,
            xmax = .data$end.x,
            ymin = .data$start.y,
            ymax = .data$end.y,
            fill = .data$cna_state,
            alpha = if(cf.gradient) .data$cf.em else NULL
        )
    ) +
        geom_rect() +
        # lines between chromosomes
        geom_vline(
            xintercept = chromStarts$chrom.start,
            colour = "grey",
            linewidth = 0.2
        ) +
        # lines between clusters
        geom_hline(
            yintercept = props,
            colour = "black",
            linewidth = 0.2
        ) +
        # chromosomes labels placement in the center of each chromosome
        scale_x_continuous(
            expand = c(0, 0),
            breaks = chromStarts$chrom.start + (chromSizes$width / 2),
            labels = as.character(chromNames)
        ) +
        # cluster labels in the center
        scale_y_continuous(
            expand = c(0, 0),
            trans = "reverse",
            limits = c(1, 0),
            breaks = props_breaks,
            labels = paste0("cluster ", unique(df$cluster), "\n", ncells, " cells")
        ) +
        # set colors for the calls and remove the name
        scale_fill_manual(
            name = "",
            values = cna.colors,
            na.value = "white",
            na.translate = FALSE,
            drop = FALSE
        ) +
        ## Set axis labels
        labs(title = title) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.title.y = element_text(size = 10, face = "bold"),
            axis.text.x = element_text(size = 6, face = "bold"),
            axis.text.y = element_text(size = 8, face = "bold"),
            axis.ticks = element_blank(),
            legend.position = "top"
        )

    # Add optional layers
    if (cf.gradient) {
        cna_plot <- cna_plot +
            scale_alpha_continuous(
                name = "cell fraction",
                range = c(0, 1),
                limits = c(0, 1),
                breaks = c(0.2, 0.4, 0.6, 0.8, 1)
            ) +
            guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2))

    }
    return(cna_plot)
}


#' Plot UMAP Coordinates from a muscadet Object
#'
#' Visualize UMAP (Uniform Manifold Approximation and Projection) coordinates
#' stored in a [`muscadet`][muscadet-class] object. The function allows coloring
#' by clusters and optionally adding cluster labels.
#'
#' @param x A [`muscadet`][muscadet-class] object containing clustering data (using
#'   [clusterMuscadet()]).
#'
#' @inheritParams heatmapMuscadet
#'
#' @param lab.x Label for the x-axis (`character` string). Default is "UMAP 1".
#' @param lab.y Label for the y-axis (`character` string). Default is "UMAP 2".
#' @param add_clusters_labels Logical. If `TRUE`, adds the cluster names as
#'   text or boxed labels using [add_labels()]. Default is `FALSE`.
#' @param point.size Numeric. Size of the points in the UMAP plot passed to
#'   [ggplot2::geom_point()]. Default is `0.5`.
#' @param legend.point.size Numeric. Size of the points in the legend of the
#'   UMAP plot. Default is `3`.
#' @param ... Additional arguments passed to [add_labels()] providing an
#'   underlying geom for label names ([`geom_text`][ggplot2::geom_text()],
#'   [`geom_label`][ggplot2::geom_label()],
#'   [`geom_text_repel`][ggrepel::geom_text_repel()], or
#'   [`geom_label_repel`][ggrepel::geom_label_repel()]).
#'
#' @return A `ggplot` object.
#'
#' @include objects.R
#'
#' @import ggplot2
#' @importFrom SeuratObject Cells
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' p <- plotUMAP(exdata_muscadet, partition = 0.3, title = "UMAP copy-number clusters")
#' print(p)
#' }
#'
plotUMAP <- function(x,
                     partition = NULL,
                     clusters = NULL,
                     colors = NULL,
                     title = "",
                     lab.x = "UMAP 1",
                     lab.y = "UMAP 2",
                     add_clusters_labels = FALSE,
                     point.size = 0.5,
                     legend.point.size = 3,
                     ...) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Validate clustering results present in the muscadet object
    stopifnot("Input object `x` must contain clustering umap data (x@clustering$umap)." = !is.null(x@clustering$umap))

    # Set default color palette for clusters if not provided
    if (is.null(colors)) {
        colors <- rep(
            c("#FABC2A", "#7FC97F", "#EE6C4D", "#39ADBD", "#BEAED4", "#FEE672", "#F76F8E",
              "#487BEA", "#B67BE6", "#F38D68", "#7FD8BE", "#F2AFC0"), 3)
    }

    # Get common cells
    common_cells <- sort(Reduce(intersect, SeuratObject::Cells(x)))

    # Filter cells based on provided `clusters` argument
    if (!is.null(clusters)) {
        # Check `clusters` cells
        stopifnot(
            "The `clusters` argument must have names, corresponding to cell names within the muscadet object `x`." = !is.null(names(clusters))
        )

        stopifnot(
            "Names of `clusters` don't match cell names within the muscadet object `x`." = length(intersect(names(clusters), common_cells)) > 0
        )

        # cells not in `clusters`
        cells_filtered <- setdiff(common_cells, names(clusters))

        if (length(cells_filtered) > 0) {
            warning(
                paste(
                    "The `clusters` argument does not contain cluster assignments for all common cells.",
                    length(cells_filtered),
                    "cells in the muscadet object `x` are filtered out."
                )
            )
            # filter cells based on `clusters` cells
            common_cells <- common_cells[common_cells %in% names(clusters)]
        }

    } else if (!is.null(partition)) {
        stopifnot(
            "The muscadet object `x` must contain the clustering results for the provided `partition`." =
                as.character(partition) %in% names(x@clustering$clusters)
        )

        clusters <- x@clustering$clusters[[as.character(partition)]]

    } else if (is.null(partition) && is.null(clusters)) {

        stopifnot(
            "The muscadet object `x` must contain the assigned clusters (x@cnacalling$clusters) (use assignClusters())
            if `partition` and `clusters` arguments are NULL." =
                !is.null(x@cnacalling$clusters)
        )

        clusters <- x@cnacalling$clusters
    }


    # Combine clusters and umap data
    df <- as.data.frame(x@clustering$umap)[common_cells, ]
    df$cluster <- as.factor(clusters[common_cells])

    p <- ggplot(df, aes(x = .data$UMAP_1 , y = .data$UMAP_2, color = .data$cluster)) +
        geom_point(size = point.size) +
        scale_color_manual(name = "Clusters", values = colors) +
        labs(title = title, x = lab.x, y = lab.y) +
        coord_fixed() +
        guides(colour = guide_legend(override.aes = list(size = legend.point.size))) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)
        )

    if (add_clusters_labels) {
        p <- add_labels(p, labels = "cluster", ...)
    }

    return(p)
}




#' Add Labels to a ggplot Object
#'
#' This function adds labels at the median position of each group to a ggplot
#' object. Labels can be added either as plain text or as label boxes, with
#' optional repelling to avoid overlaps (using the
#' [`ggprepel`][ggrepel::ggrepel] package if installed).
#'
#' @param p A ggplot object with mapping and data (`ggplot`).
#' @param labels Column name (unquoted) indicating the group label to display.
#' @param color Color of the label text (`character`). Default is `"black"`.
#' @param repel Logical. If `TRUE`(default), overlapping labels are repelled
#'   using the [`ggprepel`][ggrepel::ggrepel] package.
#' @param label.box Logical. If `TRUE` it uses a boxed label (`geom_label`)
#'   instead of plain text (`geom_text`). Default is `FALSE`.
#' @param size Size of the label text (`numeric`). Default is `2`.
#'
#' @param ... Additional arguments passed to the corresponding underlying geom:
#' - [ggplot2::geom_text()] (`repel`= `FALSE` and `label.box` = `FALSE`)
#' - [ggplot2::geom_label()] (`repel`= `FALSE` and `label.box` = `TRUE`)
#' - [ggrepel::geom_text_repel()] (`repel`= `TRUE` and `label.box` = `FALSE`)
#' - [ggrepel::geom_label_repel()] (`repel`= `TRUE` and `label.box` = `TRUE`)
#'
#' @return A ggplot2 layer  ([`geom_text`][ggplot2::geom_text()],
#'   [`geom_label`][ggplot2::geom_label()],
#'   [`geom_text_repel`][ggrepel::geom_text_repel()], or
#'   [`geom_label_repel`][ggrepel::geom_label_repel()]).
#'
#' @details The function summarizes the data by computing the median x and y
#' positions for each label group.
#'
#' If \code{repel = TRUE}, it uses
#' [`geom_text_repel()`][ggrepel::geom_text_repel()] (\code{label.box = FALSE})
#' or [`geom_label_repel()`][ggrepel::geom_label_repel()] (\code{label.box =
#' TRUE}).
#'
#' If \code{repel = FALSE}, it uses [`geom_text()`][ggplot2::geom_text()] (\code{label.box
#' = FALSE}) or [`geom_label()`][ggplot2::geom_label()] (\code{label.box = TRUE}).
#'
#' If \code{repel = TRUE}, the [`ggprepel`][ggrepel::ggrepel] package must be
#' installed.
#'
#' @importFrom rlang enquo as_label
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(ggrepel)
#'
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
#'   geom_point()
#'
#' p2 <- add_labels(
#'   p,
#'   labels = cyl,
#'   repel = TRUE,
#'   label.box = FALSE,
#'   size = 5,
#'   min.segment.length = 0
#' )
#' p2
#'
#' p3 <- add_labels(
#'     p,
#'     labels = cyl,
#'     repel = TRUE,
#'     label.box = TRUE,
#'     size = 3,
#'     box.padding = 1
#' )
#' p3
#'
#' p4 <- add_labels(
#'     p,
#'     labels = cyl,
#'     repel = FALSE,
#'     label.box = TRUE,
#'     size = 4,
#' )
#' p4
#' }
#'
add_labels <- function(
        p,
        labels,
        color = "black",
        repel = TRUE,
        label.box = FALSE,
        size = 3,
        ...
) {
    labels <- rlang::enquo(labels)
    mapping <- p$mapping
    data <- p$data

    # Check for ggrepel if repel is TRUE
    if (repel && !requireNamespace("ggrepel", quietly = TRUE)) {
        warning("Package 'ggrepel' is required for `repel` = TRUE. Install it with install.packages('ggrepel').")
        repel <- FALSE
    }

    # Extract x and y from mapping
    x_col <- rlang::as_label(mapping$x)
    y_col <- rlang::as_label(mapping$y)

    # Summarize data
    summarized_data <- data %>%
        dplyr::group_by(.data[[labels]]) %>%
        dplyr::summarize(
            med_x = median(.data[[x_col]], na.rm = TRUE),
            med_y = median(.data[[y_col]], na.rm = TRUE),
            .groups = "drop"
        )

    # Choose appropriate geom
    if (repel) {
        geom_fun <- if (label.box) ggrepel::geom_label_repel else ggrepel::geom_text_repel

        # Build layer
        ggobj <- geom_fun(
            data = summarized_data,
            mapping = aes(
                x = .data$med_x,
                y = .data$med_y,
                label = !!labels
            ),
            color = color,
            size = size,
            inherit.aes = FALSE,
            show.legend = FALSE,
            ...
        )

    } else {
        geom_fun <- if (label.box) ggplot2::geom_label else ggplot2::geom_text

        # Build layer
        ggobj <- geom_fun(
            data = summarized_data,
            mapping = aes(
                x = .data$med_x,
                y = .data$med_y,
                label = !!labels
            ),
            color = color,
            size = size,
            inherit.aes = FALSE,
            show.legend = FALSE,
            ...
        )
    }

    p + ggobj
}



#' Create heatmap and distribution plots of the different steps of computing log
#' R ratios
#'
#' This function generates heatmap and distribution plots for tumor and
#' reference cells for any step of computing log R ratios matrices. The
#' input object corresponds to the output of [computeLogRatioATAC()] or
#' [computeLogRatioRNA()] with the argument `all_steps = TRUE`.
#'
#' @param obj A list provided as output by the [computeLogRatioATAC()] or
#'   [computeLogRatioRNA()] functions with the argument `all_steps = TRUE`. It
#'   includes tumor and reference matrices at each step of the computing of log
#'   R ratio matrices.
#' @param step The step within the `obj` list to use for plotting (`character`
#'   string). It must match one of the names in `obj`.
#' @param filename File path to save the output plot (`character` string). The
#'   file format is inferred from the extension (".png", ".pdf" or ".svg").
#' @param title Title of the plot (`character` string). If `NULL`, the title is
#'   automatically generated using the provided step argument and its
#'   corresponding value name (`obj[[step]]$name`).
#' @param col_quantiles A numeric vector of length 4, specifying the quantiles
#'   to use for the color breakpoints in the heatmap (`numeric`). Either
#'   `col_quantiles` or `col_breaks` must be provided, if both are provided
#'   `col_breaks` is used. Default is `c(0.1, 0.4, 0.6, 0.9)`.
#' @param col_breaks A numeric vector of length 4, specifying custom breakpoints
#'   for the color scale in the heatmap (`numeric`). Either `col_quantiles` or
#'   `col_breaks` must be provided, if both are provided `col_breaks` is used.
#'   Default is `NULL`.
#' @param colors A character vector of 4 colors used for the color scale of the
#'   heatmap (`character` vector). Default is `c("#00008E", "white", "white",
#'   "#630000")`.
#'
#' @return The function does not return any value but saves a
#'   heatmaps-histograms plot to the specified file.
#'
#' @include objects.R
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid unit grid.grab
#' @importFrom grDevices png pdf svg
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create muscomic objects
#' atac <- CreateMuscomicObject(
#'     type = "ATAC",
#'     mat_counts = exdata_mat_counts_atac_tumor,
#'     features = exdata_peaks
#' )
#' rna <- CreateMuscomicObject(
#'     type = "RNA",
#'     mat_counts = exdata_mat_counts_rna_tumor,
#'     features = exdata_genes
#' )
#' atac_ref <- CreateMuscomicObject(
#'     type = "ATAC",
#'     mat_counts = exdata_mat_counts_atac_ref,
#'     features = exdata_peaks
#' )
#' rna_ref <- CreateMuscomicObject(
#'     type = "RNA",
#'     mat_counts = exdata_mat_counts_rna_ref,
#'     features = exdata_genes
#' )
#'
#' # Create muscadet objects
#' muscadet <- CreateMuscadetObject(
#'     omics = list(atac, rna),
#'     bulk.lrr = exdata_bulk_lrr,
#'     bulk.label = "WGS",
#'     genome = "hg38"
#' )
#' muscadet_ref <- CreateMuscadetObject(
#'     omics = list(atac_ref, rna_ref),
#'     genome = "hg38"
#' )
#'
#' # Compute log R ratios with `all_steps = TRUE`
#' obj_atac_all <- computeLogRatioATAC(
#'     matTumor = matCounts(muscadet)$ATAC,
#'     matRef = matCounts(muscadet_ref)$ATAC,
#'     peaksCoord = coordFeatures(muscadet)$ATAC,
#'     genome = muscadet$genome,
#'     minReads = 1, # low value for example subsampled datasets
#'     minPeaks = 1, # low value for example subsampled datasets
#'     all_steps = TRUE
#' )
#' names(obj_atac_all)
#'
#' # Plot heatmap and distribution of values for Step01
#' heatmapStep(obj = obj_atac_all,
#'             step = "step01",
#'             filename = file.path(tempdir(), "step01.png"),
#'             title = "Example dataset - Step 01")
#'
#' # Plot heatmap and distribution of values for all steps
#' for (step in grep("step", names(obj_atac_all), value = TRUE)) {
#'     heatmapStep(
#'         obj_atac_all,
#'         step,
#'         filename = file.path(tempdir(), paste0("ATAC_", step, ".pdf")),
#'         title = paste("ATAC -", step)
#'     )
#' }
#' }
#'
heatmapStep <- function(obj,
                        step,
                        filename,
                        title = NULL,
                        col_quantiles = c(0.1, 0.4, 0.6, 0.9),
                        col_breaks = NULL,
                        colors = c("#00008E", "white", "white", "#630000")) {

    if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required for this function. Install it with install.packages('patchwork').")
    }

    # Argument checks

    # obj
    stopifnot("The `obj` argument must be a list." = is.list(obj))

    # step
    if (!(step %in% names(obj))) {
        stop(paste(
            "The specified `step` (`",
            step,
            "`) is not found in `obj`.",
            sep = ""
        ))
    }

    # filename extension
    stopifnot("The `filename` argument must end with either .png, .pdf or .svg." = grepl(".(png|pdf|svg)$", filename))

    # col_quantiles and col_breaks
    stopifnot("Either `col_quantiles` or `col_breaks` must be provided." =
                  !(is.null(col_quantiles) && is.null(col_breaks)))
    if (!is.null(col_quantiles)) {
        stopifnot("`col_quantiles` must be a numeric vector of length 4." =
                      is.numeric(col_quantiles) && length(col_quantiles) == 4)
        stopifnot("Values in `col_quantiles` must be between 0 and 1." =
                      all(col_quantiles >= 0 & col_quantiles <= 1))
    }
    if (!is.null(col_breaks)) {
        stopifnot("`col_breaks` must be a numeric vector of length 4." =
                      is.numeric(col_breaks) && length(col_breaks) == 4)
        stopifnot("`col_breaks` must contain at least two distinct values." =
                      length(unique(col_breaks)) >= 2)
    }
    # if (!is.null(col_quantiles) && !is.null(col_breaks)) {
    #     message("Both `col_quantiles` and `col_breaks` are provided. Only `col_breaks` is used.")
    # }

    # colors
    if (!is.character(colors) || length(colors) != 4) {
        stop("The `colors` argument should be a character vector of length 4.")
    }

    # Extract the tumor and reference matrices
    matTumor <- as.matrix(obj[[step]]$matTumor)
    matRef <- as.matrix(obj[[step]]$matRef)
    name <- obj[[step]]$name

    # Generate title if NULL
    if (is.null(title)) {
        title <- paste(step, "-", name)
    }

    # Chromosome factor from coordinates
    coord <- obj$coord
    chrom <- coord[which(coord$id %in% colnames(matTumor)), "CHROM"]
    chrom <- factor(chrom, levels = unique(chrom))
    mat <- rbind(matTumor, matRef)

    # Combine matrices and calculate breaks for color scale
    if(!is.null(col_breaks)) {
        # Color scale function
        col_fun <- circlize::colorRamp2(col_breaks, colors)
    } else {
        # define breaks based on quantiles
        col_breaks <- quantile(mat, col_quantiles)
        stopifnot("The breaks defined by `col_quantiles` must contain at least two distinct values." =
                      length(unique(col_breaks)) >= 2)
        # Color scale function
        col_fun <- circlize::colorRamp2(col_breaks, colors)
    }

    # Calculate bin_width
    bin_width <- (max(mat) - min(mat)) / 1000

    ComplexHeatmap::ht_opt(message = F)

    # Create heatmaps
    ht_Tum <- ComplexHeatmap::Heatmap(
        matTumor,
        name = "Tumor",
        heatmap_legend_param = list(title = name),
        row_title = "Tumor cells",
        row_title_gp = grid::gpar(fontsize = 12),
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        column_split = chrom,
        column_title_gp = grid::gpar(fontsize = 10),
        border_gp = grid::gpar(col = "black", lwd = 1),
        heatmap_height = unit(12, "cm"),
        heatmap_width = unit(18, "cm"),
        col = col_fun,
        raster_device = "png",
        raster_quality = 3
    )

    ht_Ref <- ComplexHeatmap::Heatmap(
        matRef,
        name = "Ref",
        heatmap_legend_param = list(title = name),
        row_title = "Reference cells",
        row_title_gp = grid::gpar(fontsize = 12),
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        column_split = chrom,
        column_title_gp = grid::gpar(fontsize = 10),
        border_gp = grid::gpar(col = "black", lwd = 1),
        heatmap_height = unit(12, "cm"),
        heatmap_width = unit(18, "cm"),
        col = col_fun,
        raster_device = "png",
        raster_quality = 3
    )

    # Draw heatmaps
    pdf(file = NULL)  # Temporarily create the plot device for heatmap images
    ht_Tum_2 <- ComplexHeatmap::draw(ht_Tum, column_title = paste(name, "in tumor cells"))
    ht_Tum_grob <- grid.grab()  # Capture heatmap as grob
    ht_Ref_2 <- ComplexHeatmap::draw(ht_Ref, column_title = paste(name, "in reference cells"))
    ht_Ref_grob <- grid.grab()  # Capture heatmap as grob
    dev.off()

    # Histogram calculations and data preparation
    data_Tum <- as.data.frame(as.table(matTumor))
    data_Ref <- as.data.frame(as.table(matRef))
    colnames(data_Tum) <- colnames(data_Ref) <- c("Row", "Column", "Value")

    # Create bins based on the bin_width and compute frequency for gradient bar
    data_Tum$Value_bin <- cut(
        data_Tum$Value,
        breaks = seq(floor(min(data_Tum$Value)), ceiling(max(data_Tum$Value)), by = bin_width),
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
    )

    gradient_data_Tum <- data.frame(
        x_min = seq(min(matTumor), max(matTumor), length.out = 500)[-500],
        x_max = seq(min(matTumor), max(matTumor), length.out = 500)[-1],
        y_min = rep(-(max(
            table(data_Tum$Value_bin)
        ) * 0.02), 499),
        y_max = rep(-(max(
            table(data_Tum$Value_bin)
        ) * 0.1), 499)
    )
    gradient_data_Tum$fill <- col_fun((gradient_data_Tum$x_min + gradient_data_Tum$x_max) / 2)

    data_Ref$Value_bin <- cut(
        data_Ref$Value,
        breaks = seq(floor(min(data_Ref$Value)), ceiling(max(data_Ref$Value)), by = bin_width),
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
    )

    gradient_data_Ref <- data.frame(
        x_min = seq(min(matRef), max(matRef), length.out = 500)[-500],
        x_max = seq(min(matRef), max(matRef), length.out = 500)[-1],
        y_min = rep(-(max(
            table(data_Ref$Value_bin)
        ) * 0.02), 499),
        y_max = rep(-(max(
            table(data_Ref$Value_bin)
        ) * 0.1), 499)
    )
    gradient_data_Ref$fill <- col_fun((gradient_data_Ref$x_min + gradient_data_Ref$x_max) / 2)

    # Create histograms
    hist_Tum <- ggplot(data_Tum, aes(x = .data$Value)) +
        geom_histogram(
            aes(y = after_stat(.data$count)),
            binwidth = bin_width,
            color = "black",
            fill = "black"
        ) +
        geom_rect(
            data = gradient_data_Tum,
            aes(
                xmin = .data$x_min,
                xmax = .data$x_max,
                ymin = .data$y_min,
                ymax = .data$y_max,
                fill = .data$fill
            ),
            inherit.aes = FALSE
        ) +
        scale_fill_identity() +
        labs(
            title = paste(name, "distribution in tumor cells"),
            x = "Value",
            y = "Frequency"
        ) +
        theme_classic()

    hist_Ref <- ggplot(data_Ref, aes(x = .data$Value)) +
        geom_histogram(
            aes(y = after_stat(.data$count)),
            binwidth = bin_width,
            color = "black",
            fill = "black"
        ) +
        geom_rect(
            data = gradient_data_Ref,
            aes(
                xmin = .data$x_min,
                xmax = .data$x_max,
                ymin = .data$y_min,
                ymax = .data$y_max,
                fill = .data$fill
            ),
            inherit.aes = FALSE
        ) +
        scale_fill_identity() +
        labs(
            title = paste(name, "distribution in reference cells"),
            x = "Value",
            y = "Frequency"
        ) +
        theme_classic()

    # Combine heatmap and histogram plots
    final_plot <- patchwork::wrap_plots(c(list(ht_Tum_grob, hist_Tum), list(ht_Ref_grob, hist_Ref)), ncol = 2) +
        patchwork::plot_layout(ncol = 2, widths = c(2, 1)) +
        patchwork::plot_annotation(title = title,
                        theme = theme(plot.title = element_text(size = 16)))

    # Output to file
    if (grepl(".png", basename(filename))) {
        grDevices::png(
            filename = filename,
            width = ht_Tum_2@ht_list_param[["width"]] * 1.75,
            height = ht_Tum_2@ht_list_param[["height"]] * 2.1,
            units = "mm",
            res = 300
        )
        print(final_plot)
        dev.off()

    } else if (grepl(".pdf", basename(filename))) {
        grDevices::pdf(
            file = filename,
            width = (ht_Tum_2@ht_list_param[["width"]] * 1.75) / 25.4, # in inches
            height = (ht_Tum_2@ht_list_param[["height"]] * 2.1) / 25.4  # in inches
        )
        print(final_plot)
        dev.off()

    } else if (grepl(".svg", basename(filename))) {
        grDevices::svg(
            filename = filename,
            width = (ht_Tum_2@ht_list_param[["width"]] * 1.75) / 25.4, # in inches
            height = (ht_Tum_2@ht_list_param[["height"]] * 2.1) / 25.4  # in inches
        )
        print(final_plot)
        dev.off()
    }
}





