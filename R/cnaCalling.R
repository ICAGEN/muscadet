#' Assign a cluster assignment to a `muscadet` object
#'
#' Add the user selected cluster assignments to cells in a
#' [`muscadet`][muscadet-class] object. This function allows the user to choose the
#' cluster assignments they consider to fit the data and their requirements, or
#' cluster assignments based on other data and methods.
#'
#' @param x A [`muscadet`][muscadet-class] object (`muscadet`).
#'
#' @param partition Value specifying the clustering partition to choose from the
#'   [`muscadet`][muscadet-class] object (`numeric` or `character`). It should
#'   be either the resolution or the k number of cluster (k) used for clustering
#'   depending on the clustering method (`res_range` or `k_range` with
#'   [clusterMuscadet()]). Should be provided if `clusters` is `NULL`.
#'
#' @param clusters A custom named vector of cluster assignments (`vector`). The
#'   vector names must match cell names in the muscadet object `x`, at least
#'   cluster assignments for all common cells must be provided if
#'   `redo_imputation` is set to true, otherwise, all cells within the muscadet
#'   object `x` must be provided. Should be provided if `partition` is `NULL`.
#'
#' @param mapping Optional named vector specifying how to remap cluster values
#'   (`vector`). The names of the vector correspond to the original cluster
#'   values, and the values are the remapped cluster values. For example, `c("1"
#'   = 1, "2" = 1, "3" = 2, "4" = 3)` would merge clusters 1 and 2 into 1,
#'   cluster 3 into 2, and cluster 4 into 3.
#'
#' @param redo_imputation Logical. If `TRUE` (default), reruns the imputation
#'   process to assign clusters to cells with missing data. This ensures that
#'   imputed clusters are updated if the clustering has changed due to remapping
#'   or to the use of custom clusters.
#'
#' @inheritParams imputeClusters
#'
#' @details
#' - The clusters can be taken directly from the [`muscadet`][muscadet-class]
#' object clustering results with setting the `partition` argument (e.g.
#' `muscadet_obj$clustering$clusters[["0.8"]]` for res=`0.8`).
#' - A custom vector of cluster assignments
#' can be attributed using the `clusters` argument.
#' - Either way, the clusters assignments can be rearranged using the `mapping`
#' argument.
#'
#' @return A [`muscadet`][muscadet-class] object updated with the user chosen
#'   cluster assignments in `muscadet_obj$cnacalling$clusters`.
#'
#' @include objects.R
#'
#' @importFrom SeuratObject Cells
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' # Select clustering result for partition = 0.6
#' exdata_muscadet<- assignClusters(exdata_muscadet, partition = 0.6)
#' table(exdata_muscadet$cnacalling$clusters)
#'
#' # Assign custom clusters
#' set.seed(42)
#' cell_names <- Reduce(union, SeuratObject::Cells(exdata_muscadet))
#' n1 <- sample(1:length(cell_names), 1)
#' n2 <- length(cell_names) - n1
#' custom_clusters <- setNames(c(rep.int(1, n1), rep.int(2, n2)), cell_names)
#' table(custom_clusters)
#' exdata_muscadet <- assignClusters(exdata_muscadet, clusters = custom_clusters)
#' table(exdata_muscadet$cnacalling$clusters)
#'
#' # Assign clusters with remapping
#' # example to remap from partition=0.8 with merging of clusters 2 and 3
#' clusters <- exdata_muscadet$clustering$clusters[["0.8"]]
#' table(clusters) # 3 clusters
#' mapping <- c("1" = 1, "2" = 2, "3" = 2) # remap to 2 clusters
#'
#' exdata_muscadet <- assignClusters(exdata_muscadet, clusters = clusters, mapping = mapping)
#' table(exdata_muscadet$cnacalling$clusters)
#' # check original and remapped clusters
#' table(clusters, exdata_muscadet$cnacalling$clusters)
#'
#' exdata_muscadet <- assignClusters(exdata_muscadet, partition = 0.8, mapping = mapping)
#' table(exdata_muscadet$cnacalling$clusters)
#' # check original and remapped clusters
#' table(exdata_muscadet$clustering$clusters[["0.8"]],
#'       exdata_muscadet$cnacalling$clusters)
#'
#' # Visualize clusters on heatmap
#' heatmapMuscadet(
#'     exdata_muscadet,
#'     partition = 0.8,
#'     filename = file.path("heatmap_muscadet_res0.8.png"),
#'     title = "Example dataset | res=0.8"
#' )
#' heatmapMuscadet(
#'     exdata_muscadet,
#'     clusters = exdata_muscadet$cnacalling$clusters,
#'     filename = file.path("heatmap_muscadet_custom_res0.8.png"),
#'     title = "Example dataset | rearranged clusters from res=0.8"
#' )
#' }
#'
assignClusters <- function(x,
                           partition = NULL,
                           clusters = NULL,
                           mapping = NULL,
                           redo_imputation = TRUE,
                           knn_imp = 10) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object `x` must be of class `muscadet`." = inherits(x, "muscadet"))

    # Validate that either partition or clusters is provided, but not both
    stopifnot(
        "Either `partition` or `clusters` must be provided, but not both."
        = xor(!is.null(partition), !is.null(clusters))
    )

    # If partition is provided, validate that clustering has been performed
    if (!is.null(partition)) {
        stopifnot(
            "Clustering results are not available in the muscadet object `x`.
            Perform clustering first using `clusterMuscadet()`."
            = !is.null(x@clustering),
            "Clustering results for the chosen `partition` are not available.
            Use `clusterMuscadet()` with the `res_range`/`k_range` argument to compute the partition."
            = as.character(partition) %in% names(x@clustering$clusters)
        )
    }

    # If clusters is provided, validate its format
    if (!is.null(clusters)) {
        if (is.factor(clusters) && !is.null(names(clusters))) {
            clusters <- stats::setNames(as.vector(clusters), names(clusters))
        } else {
            stopifnot("`clusters` must be a named vector or factor."
                      = is.vector(clusters) && !is.null(names(clusters)))
        }
    }

    # Reset CNA calling results to avoid discordant outputs
    x@cnacalling <- list()

    ########

    # Helper function to impute clusters if needed
    .check_and_impute_clusters <- function(clusters_vec) {

        common_cells <- sort(Reduce(intersect, Cells(x)))
        stopifnot("`clusters` must contain at least all common cells when `redo_imputation = TRUE`." =
                      all(common_cells %in% names(clusters_vec)))

        # Remove cells missing data that needs to go through cluster imputation (cleaner)
        clusters_vec <- clusters_vec[intersect(names(clusters_vec), common_cells)]

        mat_list <- matLogRatio(x)

        imputeClusters(mat_list, clusters_vec, knn_imp = knn_imp)
    }

    # Apply mapping if provided
    if (!is.null(mapping)) {
        # Validate mapping
        stopifnot(
            "`mapping` must be a named vector." = is.vector(mapping) &&
                !is.null(names(mapping)),
            "All cluster values must have corresponding mappings in `mapping`."
            = all(as.character(unique(clusters)) %in% names(mapping))
        )

        # Get clusters from partition
        if (!is.null(partition)) {
            clusters <- x@clustering$clusters[[as.character(partition)]]
        }

        # Remap clusters
        remapped_clusters <- setNames(mapping[as.character(clusters)], names(clusters))

        stopifnot(
            "`clusters` must contain all cell names within the muscadet object `x`."
            = all(names(remapped_clusters) %in% Reduce(union, Cells(x)))
        )

        # Rerun cluster imputation (following remapping)
        if (redo_imputation) {
            remapped_clusters <- .check_and_impute_clusters(remapped_clusters)
        }

        x@cnacalling[["clusters"]] <- remapped_clusters

    } else {
        # Assign clusters without remapping - custom clusters
        if (!is.null(clusters)) {

            stopifnot(
                "`clusters` must contain all cell names within the muscadet object `x`."
                = all(names(clusters) %in% Reduce(union, Cells(x)))
            )

            # Rerun cluster imputation (in case custom clusters have been modified)
            if (redo_imputation) {
                clusters <- .check_and_impute_clusters(clusters)
            }
            x@cnacalling[["clusters"]] <- clusters

        } else if (!is.null(partition)) {
            # Assign clusters without remapping - clusters from partition
            x@cnacalling[["clusters"]] <- x@clustering$clusters[[as.character(partition)]]
        }
    }

    return(x)
}


#' Aggregate counts for `muscadet` objects
#'
#' This function combines allelic (counts at variant positions, either common
#' SNPs or individual-specific heterozygous positions) and coverage counts
#' (counts on features) from all omics per cluster for both sample and
#' reference. The resulting merged data is stored in the `cnacalling` slot of
#' the sample [`muscadet`][muscadet-class] object.
#'
#' @param x A [`muscadet`][muscadet-class] object containing sample data
#'   (`muscadet`). This object must include clustering assignments in the
#'   `cnacalling$clusters` slot.
#' @param reference A [`muscadet`][muscadet-class] object containing reference
#'   data (`muscadet`).
#' @param nor.het A logical value to specify if normal reference allele counts
#'   are modified to: total normal depth counts divided by 2, to force these
#'   positions to be heterozygous in the normal reference in allelic data (e.g.
#'   when reference cells are from a different individual but
#'   individual-specific heterozygous positions based on matched bulk sequencing
#'   data are used and thereby assumed to be heterozygous) before combining coverage and
#'   allelic data. Default is `TRUE`.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return
#' A modified [`muscadet`][muscadet-class] object corresponding to the `x`
#' [`muscadet`][muscadet-class] object, with updated `cnacalling` slot containing:
#' \itemize{
#'   \item \code{allelic.counts}: Processed allelic counts on variant positions, for all omics.
#'   \item \code{coverage.counts}: Processed coverage counts merged with the reference.
#'   \item \code{combined.counts}: Combined data for allelic and coverage counts.
#' }
#' Abbreviations:
#' - RD = Reference allele read depth
#' - AD = Alternative allele read depth
#' - DP = Total read depth
#' - TUM = tumor sample
#' - NOR = normal reference
#'
#' @importFrom data.table data.table rbindlist :=
#' @importFrom Matrix sparseMatrix colSums t
#' @importFrom dplyr left_join
#'
#' @export
#'
#' @examples
#' # Load example muscadet object
#' # data("exdata_muscadet")
#' # data("exdata_muscadet_ref")
#'
#' # Aggregate counts from all omics from both sample and reference
#' exdata_muscadet <- aggregateCounts(exdata_muscadet, exdata_muscadet_ref)
#'
aggregateCounts <- function(x,
                            reference,
                            nor.het = TRUE,
                            quiet = FALSE) {


    # Validate input: x and reference must be muscadet objects
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))
    stopifnot("Input object 'reference' must be of class 'muscadet'." = inherits(reference, "muscadet"))

    # Ensure clustering data is available in the sample muscadet object
    stopifnot(
        "Cluster assignments not found in the 'x' muscadet object. Use 'assignClusters()' to add them."
        = !is.null(x@cnacalling$clusters)
    )

    # Ensure allele counts data is available in the sample and reference muscadet object
    stopifnot(
        "Input object 'x' must contain allele counts data in the 'allelic' slot of omics. Use 'addAlleleCounts()' to add them."
        = all(unlist(lapply(x@omics, function(omic) {
            !is.null(omic@allelic$mat.RD)
        })))
    )
    stopifnot(
        "Input object 'reference' must contain allele counts data in the 'allelic' slot of omics. Use 'addAlleleCounts()' to add them."
        = all(unlist(lapply(reference@omics, function(omic) {
            !is.null(omic@allelic$mat.RD)
        })))
    )

    # Reset CNA calling results to avoid discordant outputs
    x@cnacalling <- list(clusters = x@cnacalling$clusters)

    if (!quiet) {
        table <- table(x@cnacalling$clusters)
        message(paste("Clusters used:", paste(names(table), paste0("(", table, " cells", ")"), collapse = ", ")))
    }

    # Helper functions ---------------------------------------------------------

    ## Aggregate counts from cells to clusters counts in sparse matrices ---
    aggregate_sparse_counts <- function(mat, cells, clusters) {

        # Get integer row indices directly
        row_idx <- match(names(clusters), cells)
        row_idx <- row_idx[!is.na(row_idx)]  # keep valid rows

        # Cluster assignments in correct order
        cl <- clusters[cells[row_idx]]

        # Unique cluster labels
        ucl <- unique(cl)

        # Build cluster x cell indicator matrix
        I <- Matrix::sparseMatrix(
            i = match(cl, ucl),
            j = seq_along(cl),
            x = 1L,
            dims = c(length(ucl), length(cl))
        )

        # Multiplication
        out <- I %*% mat[row_idx, , drop = FALSE]

        rownames(out) <- ucl

        return(out)
    }

    ## Extract per-cluster allelic counts for sample cells ---
    extract_allelic_for_sample <- function(x) {

        # Bind variables (for data.table use)
        RD.omic <- AD.omic <- DP.omic <- NULL
        cluster <- omic <- CHROM <- POS <- DP <- NULL

        # Extract clusters
        clusters <- x@cnacalling$clusters

        # Sum across cells per cluster and per omic
        res <- data.table::rbindlist(lapply(names(x@omics), function(omic_name) {
            allelic <- x@omics[[omic_name]]@allelic
            RD <- allelic$mat.RD
            AD <- allelic$mat.AD
            cells <- rownames(RD)
            coord <- allelic$coord.vars

            # Keep only tumor cells with cluster assignments
            keep_cells <- cells %in% names(clusters)

            # Keep only loci with information
            keep_loci <- which(colSums(RD != 0) > 0 | colSums(AD != 0) > 0)

            RD <- RD[keep_cells, keep_loci, drop = FALSE]
            AD <- AD[keep_cells, keep_loci, drop = FALSE]
            coord <- coord[keep_loci, ]

            # aggregate sparse
            agg_RD <- aggregate_sparse_counts(RD, cells[keep_cells], clusters)
            agg_AD <- aggregate_sparse_counts(AD, cells[keep_cells], clusters)
            agg_DP <- agg_RD + agg_AD

            # convert to long format
            ncl <- nrow(agg_DP)
            nloci <- ncol(agg_DP)

            dt <- data.table::data.table(
                cluster = rep(rownames(agg_DP), each = nloci),
                id = rep(coord$id, times = ncl),
                CHROM = rep(coord$CHROM, times = ncl),
                POS = rep(coord$POS, times = ncl),
                omic = omic_name,
                RD.omic = as.vector(Matrix::t(agg_RD)),
                AD.omic = as.vector(Matrix::t(agg_AD)),
                DP.omic = as.vector(Matrix::t(agg_DP))
            )
            dt[DP.omic > 0]
            return(dt)
        }))

        # Sum across omics
        res <- res[, list(
            RD = sum(RD.omic, na.rm = TRUE),
            AD = sum(AD.omic, na.rm = TRUE),
            DP = sum(DP.omic, na.rm = TRUE)
        ), by = list(cluster, id, omic, CHROM, POS)][DP > 0]

        # Sort data
        res <- res[order(res$CHROM, res$POS, res$cluster), ]

        return(as.data.frame(res))
    }

    ## Extract allelic counts for reference cells ---
    extract_allelic_for_ref <- function(reference) {

        # Bind variables (for data.table use)
        RD.omic <- AD.omic <- DP.omic <- NULL
        omic <- CHROM <- POS <- DP <- NULL

        # Sum across cells per omic
        res <- data.table::rbindlist(lapply(names(reference@omics), function(omic_name) {

            allelic <- reference@omics[[omic_name]]@allelic
            coord <- allelic$coord.vars

            data.table::data.table(
                id = coord$id,
                CHROM = coord$CHROM,
                POS = coord$POS,
                omic = omic_name,
                RD.omic = colSums(allelic$mat.RD),
                AD.omic = colSums(allelic$mat.AD),
                DP.omic = colSums(allelic$mat.RD + allelic$mat.AD)
            )
        }))

        # Sum across omics
        res <- res[, list(
            RD = sum(RD.omic, na.rm = TRUE),
            AD = sum(AD.omic, na.rm = TRUE),
            DP = sum(DP.omic, na.rm = TRUE)
        ), by = list(id, omic, CHROM, POS)]

        # Sort data
        res <- res[order(res$CHROM, res$POS), ]

        return(as.data.frame(res))
    }

    ## Extract per-cluster coverage counts for sample cells ---
    extract_coverage_for_sample <- function(x) {

        # Bind variables (for data.table use)
        DP.omic <- cluster <- omic <- CHROM <- POS <- DP <- NULL

        # Extract clusters
        clusters <- x@cnacalling$clusters

        # Sum across cells per cluster and per omic
        res <- data.table::rbindlist(lapply(names(x@omics), function(omic_name) {

            coverage <- x@omics[[omic_name]]@coverage
            mat <- coverage$counts$mat
            cells <- rownames(mat)
            coord <- coverage$counts$coord.features

            ## Keep only tumor cells with clusters
            keep_cells <- cells %in% names(clusters)

            ## Remove features with zero coverage across cells
            keep_loci <- which(colSums(mat != 0) > 0)

            mat <- mat[keep_cells, keep_loci, drop = FALSE]
            coord <- coord[keep_loci, ]

            ## Aggregate counts by cluster
            agg_DP <- aggregate_sparse_counts(mat, cells[keep_cells], clusters)

            ncl <- nrow(agg_DP)
            nloci <- ncol(agg_DP)

            ## Build long table
            dt <- data.table::data.table(
                cluster = rep(rownames(agg_DP), each = nloci),
                id = rep(coord$id, times = ncl),
                CHROM = rep(coord$CHROM, times = ncl),
                POS = rep(round((coord$start + coord$end)/2), times = ncl),
                omic = omic_name,
                DP.omic  = as.vector(Matrix::t(agg_DP))   # row-major flattening
            )

            dt[DP.omic > 0]
            return(dt)

        }))

        # Sum across omics
        res <- res[, list(
            DP = sum(DP.omic, na.rm = TRUE)
        ), by = list(cluster, id, omic, CHROM, POS)][DP > 0]

        # Sort data
        res <- res[order(res$CHROM, res$POS, res$cluster), ]

        return(as.data.frame(res))
    }

    ## Extract coverage counts for reference cells ---
    extract_coverage_for_ref <- function(reference) {

        # Bind variables (for data.table use)
        DP.omic <- cluster <- omic <- CHROM <- POS <- DP <- NULL

        # Sum across cells per omic
        res <- data.table::rbindlist(lapply(names(reference@omics), function(omic_name) {

            coverage <- reference@omics[[omic_name]]@coverage
            mat <- coverage$counts$mat
            coord <- coverage$counts$coord.features

            # filter features
            coord <- coord[coord$id %in% colnames(mat), ]
            mat <- mat[, coord$id]

            dt <- data.table::data.table(
                id = coord$id,
                CHROM = coord$CHROM,
                POS = round((coord$start + coord$end) / 2),
                omic = omic_name,
                DP.omic = colSums(mat)
            )
            return(dt)
        }))

        # Sum across omics
        res <- res[, list(DP = sum(DP.omic, na.rm = TRUE)), by = list( id, omic, CHROM, POS)]

        # Sort data
        res <- res[order(res$CHROM, res$POS), ]

        return(as.data.frame(res))
    }


    # Allelic Data Processing --------------------------------------------------
    if (!quiet) {
        message("Allelic data processing...")
    }

    # Extract allele counts for all omics in both sample and reference
    allele_sample <- extract_allelic_for_sample(x)
    allele_ref <- extract_allelic_for_ref(reference)

    # Join sample and reference allelic data
    var <- dplyr::left_join(
        allele_sample,
        allele_ref,
        by = c("CHROM", "POS", "id", "omic"),
        suffix = c(".TUM", ".NOR")
    )
    var <- var[order(var$CHROM, var$POS, var$cluster), ]
    var$signal <- "allelic"
    var <- var[, c("id", "CHROM", "POS", "cluster", "signal", "omic", "RD.NOR", "AD.NOR", "DP.NOR", "RD.TUM", "AD.TUM", "DP.TUM")]

    x@cnacalling[["allelic.counts"]] <- as.data.frame(var)

    # Auto set normal variant positions as heterozygous before combining
    if(nor.het == TRUE) {
        var[, "RD.NOR"] <- round(var[, "DP.NOR"] / 2, 0)
    }


    # Coverage Data Processing -------------------------------------------------
    if (!quiet) {
        message("Coverage data processing...")
    }

    # Extract coverage counts for all omics in both sample and reference
    coverage_sample <- extract_coverage_for_sample(x)
    coverage_ref <- extract_coverage_for_ref(reference)

    # Join sample and reference coverage data
    cov <- dplyr::left_join(
        coverage_sample,
        coverage_ref,
        by = c("CHROM", "POS", "id", "omic"),
        suffix = c(".TUM", ".NOR")
    )
    cov <- cov[order(cov$CHROM, cov$POS, cov$cluster), ]
    cov$signal <- "coverage"
    cov$RD.NOR <- cov$DP.NOR
    cov$RD.TUM <- cov$DP.TUM
    cov <- cov[, c("id", "CHROM", "POS", "cluster", "signal", "omic", "RD.NOR", "DP.NOR", "RD.TUM", "DP.TUM")]

    x@cnacalling[["coverage.counts"]] <- as.data.frame(cov)


    # Combine Allelic and Coverage Data ----------------------------------------
    if (!quiet) {
        message("Combining allelic and coverage data...")
    }

    # Format data
    cols <- c("CHROM", "POS", "DP.NOR", "RD.NOR", "DP.TUM", "RD.TUM", "cluster", "signal", "omic", "id")
    var <- var[, cols]
    var <- unique(var)
    cov <- cov[, cols]
    cov <- unique(cov)
    colnames(var) <- c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD", "cluster", "signal", "omic", "id")
    colnames(cov) <- colnames(var)

    # Make sure the levels match for binding
    levels(var$Chromosome) <- union(levels(var$Chromosome), levels(cov$Chromosome))
    levels(cov$Chromosome) <- union(levels(var$Chromosome), levels(cov$Chromosome))

    # Combine both data
    combined <- rbind(var, cov)
    combined <- combined[order(combined$Chromosome, combined$Position), ]
    rownames(combined) <- NULL

    x@cnacalling[["combined.counts"]] <- as.data.frame(combined)

    return(x)
}


#' Copy Number Alteration (CNA) Calling from muscadet object
#'
#' Performs copy number alteration (CNA) analysis on a
#' [`muscadet`][muscadet-class] object by processing allelic and coverage counts
#' across clusters and evaluating cell fractions and copy numbers.
#'
#' @param x A [`muscadet`][muscadet-class] object. Must contain:
#'   - Clustering assignments in the `cnacalling$clusters` slot
#'   (use [assignClusters()]).
#'   - Combined allelic and coverage counts per cluster in the
#'   `cnacalling$combined.counts` slot (use [aggregateCounts()]).
#' @param omics.coverage A vector of omics names to select for coverage log R
#'   ratio data. RECOMMENDED: select "ATAC" when ATAC and RNA omics are
#'   available, the ATAC coverage (DNA) signal is less noisy than RNA signal. By
#'   default, `NULL` selects all available data.
#' @param depthmin.a.clusters Minimum allelic depth per clusters in tumor cells
#'   (default: 30).
#' @param depthmin.c.clusters Minimum coverage depth per clusters in tumor cells
#'   (default: 50).
#' @param depthmin.a.allcells Minimum allelic depth for all tumor cells
#'   (default: 30).
#' @param depthmin.c.allcells Minimum coverage depth for all tumor cells
#'   (default: 50).
#' @param depthmin.c.nor Minimum coverage depth for normal sample (default: 1000).
#' @param depthmax.nor Optional. Maximum depth for normal sample (default: `NULL`).
#' @param het.thresh VAF (Variant Allele Frequency) threshold to call variant
#'   positions heterozygous for [preProcSample2()] (default: 0.25).
#' @param snp.nbhd Window size for selecting SNP loci to reduce serial
#'   correlation for [preProcSample2()] (default: 250).
#' @param hetscale Logical value indicating whether log odds ratio (logOR)
#'   should be scaled to give more weight in the test statistics for
#'   segmentation and clustering [preProcSample2()]. (default: `TRUE`)
#' @param cval1 Critical value for segmentation for [preProcSample2()] (default:
#'   25).
#' @param cval2 Critical value for segmentation for [facets::procSample()]
#'   (default: 150).
#' @param min.nhet Minimum number of heterozygous positions in a segment for
#'   [facets::procSample()] and [facets::emcncf()] (default: 5).
#' @param clonal.thresh Threshold of minimum cell proportion to label a
#'   segment as clonal (default: 0.9).
#' @param cf.thresh Numeric threshold to set the minimum cell fraction (CF) allowed. CNA segments
#'   with a CF below this value are considered neutral (i.e., 2:1). The default
#'   is `0.5`. Can be set to `0` or `NULL` to disable this threshold and retain all segments
#'   regardless of CF.
#' @param dist.breakpoints Minimum distance between breakpoints to define
#'   distinct segments (see [getSegConsensus()]) (default: 1e6).
#' @param minoverlap Minimum distance between a cluster specific segment and a
#'   consensus segment for them to overlap (see [annotateSegments()]) (default:
#'   1e6).
#' @param ploidy Specifies ploidy assumption: `"auto"`, `"median"`, or numeric
#'   value (default: `"auto"`).
#' @param dipLogR Optional. Numeric value specifying an expected log ratio for
#'   diploid regions to use for the analysis. If `NULL`, by default, the diploid
#'   log ratio is estimated automatically from the data.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return
#' A modified [`muscadet`][muscadet-class] object with added CNA analysis results in
#' the `cnacalling` slot, including: filtered counts and positions, segmentation
#' data for clusters and all cells, consensus segments across clusters based on
#' breakpoints, diploid log R ratio, purity and ploidy.
#'
#' Details of the cnacalling slot:
#' \itemize{
#'   \item \code{combined.counts.filtered}: Filtered counts per clusters.
#'   \item \code{combined.counts.allcells}: Counts summed for all cells (no
#'   cluster distinction).
#'   \item \code{combined.counts.allcells.filtered}: Filtered counts summed for
#'   all cells (no cluster distinction).
#'
#'   \item \code{positions}: Data frame of positions from the per cluster
#'   analysis. Positions in rows and associated data in columns: `chrom`,
#'   `maploc` (position), `rCountT` (read count in tumor), `rCountN` (read count
#'   in normal), `vafT` (variant allele frequency in tumor), `vafN` (variant
#'   allele frequency in normal), `cluster` (cluster id), `signal` (whether the
#'   counts come from coverage or allelic data), `het` (heterozygous status),
#'   `keep` (whether to keep position), `gcpct` (GC percentage), `gcbias` (GC
#'   bias correction), `cnlr` (log R ratio), `valor` (log odds ratio), `lorvar`
#'   (variance of log odds ratio), `seg0`, `seg_ori` (segment original id within
#'   each cluster), `seg` (segment id), `segclust` (cluster of segments id),
#'   `vafT.allcells` (vairiant allele frequency in all tumor cells), `colVAR`
#'   (integer for allelic position color depending on `vafT.allcells`).
#'
#'   \item \code{segments}: Data frame of segments from the per cluster
#'   analysis. Segments in rows and associated data in columns: `chrom`, `seg`
#'   (segment id), `num.mark` (number of positions in segment), `nhet` (number
#'   of heterezygous positions in segment), `cnlr.median` (segment log R ratio
#'   median), `mafR` (segment square of expected log odds ratio), vafT.median
#'   (segment variant allele frequency median), `cluster` (cluster id),
#'   `seg_ori` (segment original id within each cluster), `segclust` (cluster of
#'   segments id), `cnlr.median.clust` (segment cluster log R ratio median),
#'   `mafR.clust` (segment cluster square of expected log odds ratio), `cf`
#'   (cell fraction), `tcn` (total copy number), `lcn` (lower copy number),
#'   `start`, `end`, `cf.em` (cell fraction computed with EM algorithm),
#'   `tcn.em`, (total copy number computed with EM algorithm), `lcn.em` (lower
#'   copy number computed with EM algorithm).
#'
#'   \item \code{positions.allcells}: Same as `positions` but from the all cells analysis.
#'   \item \code{segments.allcells}: Same as `segments` but from the all cells analysis.
#'
#'   \item \code{consensus.segs}: Data frame of unique consensus segments across
#'   clusters, with the `cna` (`logical`) and `cna_clonal` (`logical`)
#'   information.
#'
#'   \item \code{table}: Data frame of consensus segments across clusters with
#'   associated information per cluster in columns: `chrom`, `start`, `end`, `id`,
#'   `cluster`, `cf.em` (cell fraction computed with EM algorithm), `tcn.em`
#'   (total copy number computed with EM algorithm), `lcn.em` (lower copy number
#'   computed with EM algorithm), `ncells` (number of cells in cluster),
#'   `prop.cluster` (proportion of cells per cluster), `gnl` (gain;neutral;loss :
#'   1;0;-1), `loh` (loss of heterozygosity status), `state` (state of segments),
#'   `cna` (whether the segment is a CNA), `cna_state` (state of CNA segments),
#'   `prop.tot` (proportion of cells with the same state per segment),
#'   `state_clonal` (state of the segment if its `prop.tot` is above
#'   `clonal.thresh`), `cna_clonal` (whether the segment is a clonal CNA),
#'   `cna_clonal_state` (state of clonal CNA segments).
#'
#'   \item \code{ncells}: Vector of number of cells per cluster.
#'   \item \code{dipLogR.clusters}: Diploid log R ratio estimated during the per cluster analysis.
#'   \item \code{dipLogR.allcells}: Diploid log R ratio estimated during the all cells analysis.
#'   \item \code{purity.clusters}: Purity estimated during the per cluster analysis.
#'   \item \code{purity.allcells}: Purity estimated during the all cells analysis.
#'   \item \code{ploidy.clusters}: Ploidy estimated during the per cluster analysis.
#'   \item \code{ploidy.allcells}: Ploidy estimated during the all cells analysis.
#'   \item \code{ploidy}: Ploidy used for the CNA analysis.
#' }
#'
#' @seealso [assignClusters()], [aggregateCounts()],[preProcSample2()]
#'
#' @source This function uses several functions from the
#'   [`facets`][facets::facets-package] package, including:
#'   - [facets::clustersegs()]
#'   - [facets::emcncf()]
#'   - [facets::findDiploidLogR()]
#'   - [facets::fitcncf()]
#'   - [facets::procSample()]
#'   - [facets::procSnps()]
#'   - adapted function [preProcSample2()]
#'
#'   Seshan VE, Shen R (2021). _facets: Cellular Fraction and Copy Numbers from
#'   Tumor Sequencing_. R package version 0.6.2,
#'   [https://github.com/mskcc/facets](https://github.com/mskcc/facets).
#'
#' @references
#' \describe{
#'   \item{[`facets`][facets::facets-package] package}{Shen R, Seshan VE. FACETS:
#'   allele-specific copy number and clonal heterogeneity analysis tool for
#'   high-throughput DNA sequencing. Nucleic Acids Res. 2016 Sep 19;44(16):e131.
#'   doi: [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
#'   27270079; PMCID: PMC5027494.}
#' }
#'
#' @include objects.R
#'
#' @import dplyr
#' @import facets
#' @import pctGCdata
#' @importFrom GenomicRanges GRanges
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' library("facets")
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' exdata_muscadet <- cnaCalling(
#'     exdata_muscadet,
#'     omics.coverage = "ATAC",
#'     depthmin.a.clusters = 3, # set low thresholds for example data
#'     depthmin.c.clusters = 5,
#'     depthmin.a.allcells = 3,
#'     depthmin.c.allcells = 5,
#'     depthmin.c.nor = 0
#' )
#'
cnaCalling <- function(
        x,
        omics.coverage = NULL,
        depthmin.a.clusters = 30,
        depthmin.c.clusters = 50,
        depthmin.a.allcells = 30,
        depthmin.c.allcells = 50,
        depthmin.c.nor = 1000,
        depthmax.nor = NULL,
        het.thresh = 0.25,
        snp.nbhd = 250,
        hetscale = TRUE,
        cval1 = 25,
        cval2 = 150,
        min.nhet = 5,
        clonal.thresh = 0.9,
        cf.thresh = 0.5,
        dist.breakpoints = 1e6,
        minoverlap = 1e6,
        ploidy = "auto",
        dipLogR = NULL,
        quiet = FALSE
) {

    # ==========================================================================
    # Condition messages if quiet = FALSE
    msg <- function(...) if (!quiet) message(...)

    # ==========================================================================
    # Argument validation
    stopifnot("Input object `x` must be of class 'muscadet'." = inherits(x, "muscadet"))
    stopifnot(
        is.numeric(depthmin.a.clusters),
        is.numeric(depthmin.c.clusters),
        is.numeric(depthmin.a.allcells),
        is.numeric(depthmin.c.allcells),
        is.numeric(depthmin.c.nor),
        is.numeric(het.thresh),
        is.numeric(snp.nbhd),
        is.numeric(min.nhet),
        is.numeric(cval1),
        is.numeric(cval2),
        is.numeric(clonal.thresh),
        is.numeric(dist.breakpoints)
    )
    stopifnot(
        "The `ploidy` argument must be either numeric or 'median' or 'auto'." =
            is.character(ploidy) && ploidy %in% c("auto", "median") || is.numeric(ploidy)
    )

    # ==========================================================================
    # Get internal functions from facets package (to avoid using `:::`)
    clustersegs <- utils::getFromNamespace("clustersegs", "facets")
    findDiploidLogR <- utils::getFromNamespace("findDiploidLogR", "facets")
    fitcncf <- utils::getFromNamespace("fitcncf", "facets")
    emcncf <- utils::getFromNamespace("emcncf", "facets")
    procSnps <- utils::getFromNamespace("procSnps", "facets")

    # ==========================================================================
    # INPUTS -------------------------------------------------------------------

    # Extract required inputs
    gbuild <- x$genome
    rcmat <- x@cnacalling$combined.counts
    chromlevels <- levels(rcmat$Chromosome)
    clusters <- as.character(x@cnacalling$clusters)
    ncells <- table(clusters)

    # Determine chromosome number for X based on genome build
    nX <- switch(
        gbuild,
        "hg19" = 23,
        "hg38" = 23,
        "mm10" = 20,
        stop("Unsupported genome: ", gbuild)
    )

    if (is.null(rcmat) || is.null(clusters)) {
        stop("The input muscadet object `x` must contain valid 'x@cnacalling$clusters' and 'x@cnacalling$combined.counts'.")
    }

    # Set depthmax.nor to maximum value if NULL
    if (is.null(depthmax.nor)) depthmax.nor <- max(rcmat$NOR.DP, na.rm = TRUE) + 1

    # Reset CNA calling results to avoid discordant outputs
    x@cnacalling <- list(
        clusters = x@cnacalling$clusters,
        allelic.counts = x@cnacalling$allelic.counts,
        coverage.counts = x@cnacalling$coverage.counts,
        combined.counts = x@cnacalling$combined.counts
    )
    # ==========================================================================


    # ==========================================================================
    # 1. PER CLUSTER ANALYSIS --------------------------------------------------
    # ==========================================================================

    msg("- Analysis per cluster -")

    msg("Initial number of positions: ", nrow(rcmat))

    allelic_raw  <- dplyr::filter(rcmat, .data$signal == "allelic")
    coverage_raw <- dplyr::filter(rcmat, .data$signal == "coverage")

    msg("Initial number of allelic positions: ", nrow(allelic_raw))
    msg("Initial number of coverage positions: ", nrow(coverage_raw))

    # ==========================================================================
    # Process allelic counts

    if (length(unique(allelic_raw$omic)) > 1) {

        msg("Integrating omics...")

        cols <- colnames(allelic_raw)[1:8]
        allelic_raw <- allelic_raw %>%
            dplyr::group_by(.data$Chromosome, .data$Position, .data$cluster, .data$signal) %>%
            dplyr::summarise(
                NOR.DP = sum(.data$NOR.DP, na.rm = TRUE),
                NOR.RD = sum(.data$NOR.RD, na.rm = TRUE),
                TUM.DP = sum(.data$TUM.DP, na.rm = TRUE),
                TUM.RD = sum(.data$TUM.RD, na.rm = TRUE),
                .groups = "drop"
            )
        allelic_raw <- allelic_raw[, cols] # reorder columns
    } else {
        allelic_raw <- allelic_raw %>% dplyr::select(-"omic", -"id")
    }

    # Filter allelic counts
    msg("Filtering allelic positions: tumor depth >= ", depthmin.a.clusters, " reads")

    allelic_rcmat <- allelic_raw %>%
        dplyr::filter(.data$TUM.DP >= depthmin.a.clusters) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)

    allelic_rcmat <- na.omit(allelic_rcmat)

    # ==========================================================================
    # Process coverage counts

    # Optionally select omics
    if (!is.null(omics.coverage)) {

        msg("Selecting coverage data from omic(s): ", omics.coverage)

        stopifnot(
            "The `omics.coverage` argument must match one or more omic name" =
                omics.coverage %in% unique(coverage_raw$omic)
        )

        coverage_raw <- coverage_raw %>% dplyr::filter(.data$omic %in% omics.coverage)
    }

    # Filter coverage counts
    msg("Filtering coverage positions: tumor depth >= ", depthmin.c.clusters, " reads")
    msg("Filtering coverage positions: normal depth >= ", depthmin.c.nor, " reads")

    coverage_rcmat <- coverage_raw %>%
        dplyr::filter(
            .data$TUM.DP >= depthmin.c.clusters,
            .data$NOR.DP >= depthmin.c.nor
        ) %>%
        dplyr::select(-c("omic", "id")) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)

    coverage_rcmat <- na.omit(coverage_rcmat)

    # ==========================================================================
    # Recombine counts

    msg("Allelic positions kept: ", nrow(allelic_rcmat))
    msg("Coverage positions kept: ", nrow(coverage_rcmat))

    rcmat_filtered <- dplyr::bind_rows(allelic_rcmat, coverage_rcmat) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)

    msg("Final number of positions: ", nrow(rcmat_filtered))

    x@cnacalling[["combined.counts.filtered"]] <- rcmat_filtered


    # ==========================================================================
    # GET SEGMENTS AND ASSOCIATED DATA -----------------------------------------

    msg("Performing segmentation per cluster...")

    oo_list <- list()
    for (i in unique(clusters)) {
        rcmat_clus <- dplyr::filter(rcmat_filtered, .data$cluster == i)
        if (length(unique(rcmat_clus$signal)) == 2) {
            xx <- preProcSample2(
                rcmat_clus,
                het.thresh = het.thresh,
                cval = cval1,
                snp.nbhd = snp.nbhd,
                gbuild = gbuild,
                ndepth = 1,
                ndepthmax = depthmax.nor
            )
            oo <- facets::procSample(xx,
                                     cval = cval2,
                                     min.nhet = min.nhet,
                                     dipLogR = dipLogR)
            oo[["clusters"]] <- as.character(i)
            oo_list[[as.character(i)]] <- oo
        }
    }

    # ==========================================================================
    # Filtered clusters without enough data
    processed_clusters <- names(oo_list)
    initial_clusters <- as.character(unique(clusters))
    skipped_clusters <- setdiff(initial_clusters, processed_clusters)

    if (length(skipped_clusters) >= 1) {
        msg("Clusters without enough data are filtered out: ", paste(skipped_clusters, collapse = ", "))

        # Update clusters and ncells
        clusters <- clusters[which(clusters %in% processed_clusters)]
        ncells <- table(clusters)

        x@cnacalling[["clusters"]] <- clusters
        x@cnacalling[["skipped_clusters"]] <- skipped_clusters
    }
    # ==========================================================================

    # Add cluster column to output "out" table
    oo_list <- lapply(oo_list, function(x) {
        x$out$cluster <- x[["clusters"]]
        return(x)
    })
    # Bind "out" table of segments from all clusters
    out_full <- lapply(oo_list, function(x) {
        x$out[ , c(1:6, 13)]
    })
    out_full <- do.call(rbind, out_full)
    # Rename segments with unique number
    out_full$seg_ori <- out_full$seg
    out_full$seg <- 1:nrow(out_full)

    # Bind "jointseg" table of segment data from all clusters
    jseg_full <- lapply(oo_list, function(x) {
        x$jointseg[, 1:17]
    })
    jseg_full <- do.call(rbind, jseg_full)
    # Add the new unique segment number
    jseg_full$seg_ori <- jseg_full$seg
    jseg_full <- dplyr::select(jseg_full, colnames(jseg_full)[colnames(jseg_full) != "seg"])
    jseg_full <- dplyr::left_join(jseg_full,
                                  out_full[, c("seg", "cluster", "seg_ori")],
                                  by =c("cluster", "seg_ori"))

    # Cluster segments and add the segclust IDs to the jointseg complete data table
    outclust_full <- tryCatch(
        expr = {
            clustersegs(out_full, jseg_full, min.nhet = min.nhet)
        },
        error = function(e) {
            message("Error using facets::clustersegs(), consider modifying the `min.nhet` argument.")
            message(e)
            stop()
        }
    )

    jseg_full <- dplyr::left_join(jseg_full,
                                  outclust_full[, c("seg", "segclust")],
                                  by = "seg")

    # ==========================================================================
    # FIND DIPLOID Log ratio ---------------------------------------------------

    oo_full <- findDiploidLogR(outclust_full, jseg_full$cnlr)

    if (is.null(dipLogR)) {
        msg("Finding diploid log R ratio on clusters...")

        dipLogR <- oo_full$dipLogR

    }
    msg("Diploid log R ratio = ", round(dipLogR, 3))


    # ==========================================================================
    # COMPUTE CF, TCN, LCN -----------------------------------------------------

    msg("Computing cell fractions and copy numbers on clusters...")

    # Compute CF Cell fraction and CN Copy Number (L for lower, T for total)
    outclust_full <- fitcncf(outclust_full, dipLogR, nX)
    oo_full <- c(
        list(
            jointseg = jseg_full,
            out = outclust_full,
            nX = nX,
            chromlevels = chromlevels,
            clusters = unique(outclust_full$cluster)
        ),
        oo_full[-1]
    )
    # Compute CF Cell fraction and CN Copy Number with EM algorithm
    fit_full <- emcncf(oo_full, min.nhet = min.nhet)
    # Add CF and CN computed with EM to table
    out_full <- cbind(oo_full$out, fit_full$cncf[, c("start", "end", "cf.em", "tcn.em", "lcn.em")])

    # Compute median of vafT per segment
    jseg_full <- jseg_full %>%
        dplyr::group_by(.data$seg) %>%
        dplyr::mutate(vafT.median = median(.data$vafT[.data$signal == "allelic"], na.rm = TRUE))
    out_full <- dplyr::left_join(out_full, unique(jseg_full[, c("seg", "vafT.median")]), by="seg")

    jseg_full <- jseg_full[, colnames(jseg_full) != "vafT.median"]
    out_full <- out_full[, colnames(out_full)[c(1:6, 20, 7:19)]]

    # Reorder
    jseg_full <- jseg_full[order(jseg_full$chrom, jseg_full$maploc), ]
    out_full <- out_full[order(out_full$chrom, out_full$start), ]

    # Rename chr X instead of 23
    out_full$chrom[out_full$chrom == 23] <- "X"
    jseg_full$chrom[jseg_full$chrom == 23] <- "X"


    # ==========================================================================
    # CF THRESHOLD -------------------------------------------------------------

    if (!is.null(cf.thresh)) {
        # Correct CN to diploid / neutral is CF below threshold
        out_full <- dplyr::mutate(
            out_full,
            tcn.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 2, .default = .data$tcn.em),
            lcn.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 1, .default = .data$lcn.em),
            cf.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 1, .default = .data$cf.em))
    }


    # ==========================================================================
    # 2. ALL CELLS ANALYSIS ----------------------------------------------------
    # ==========================================================================

    msg("- Analysis on all cells -")

    # ==========================================================================
    # Process allelic counts

    msg("Aggregating allelic counts of all cells...")

    cols <- colnames(allelic_raw)[1:8]
    allelic_raw_allcells <- allelic_raw %>%
        dplyr::group_by(.data$Chromosome, .data$Position, .data$signal) %>%
        dplyr::summarise(
            NOR.DP = sum(.data$NOR.DP, na.rm = TRUE),
            NOR.RD = sum(.data$NOR.RD, na.rm = TRUE),
            TUM.DP = sum(.data$TUM.DP, na.rm = TRUE),
            TUM.RD = sum(.data$TUM.RD, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(cluster = "allcells")
    allelic_raw_allcells <- allelic_raw_allcells[, cols] # reorder columns

    msg("Filtering allelic positions: tumor depth >= ", depthmin.a.allcells, " reads")

    allelic_rcmat_allcells <- allelic_raw_allcells %>%
        dplyr::filter(.data$TUM.DP >= depthmin.a.allcells) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)

    allelic_rcmat_allcells <- na.omit(allelic_rcmat_allcells)

    msg("Allelic positions kept: ", nrow(allelic_rcmat_allcells))

    # ==========================================================================
    # Process coverage counts

    msg("Aggregating coverage counts of all cells...")

    cols <- colnames(coverage_raw)[1:8]
    coverage_raw_allcells <- coverage_raw %>%
        dplyr::group_by(.data$Chromosome, .data$Position, .data$signal) %>%
        dplyr::summarise(
            NOR.DP = sum(.data$NOR.DP, na.rm = TRUE),
            NOR.RD = sum(.data$NOR.RD, na.rm = TRUE),
            TUM.DP = sum(.data$TUM.DP, na.rm = TRUE),
            TUM.RD = sum(.data$TUM.RD, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(cluster = "allcells")
    coverage_raw_allcells <- coverage_raw_allcells[, cols] # reorder columns

    msg("Filtering coverage positions: tumor depth >= ", depthmin.c.allcells, " reads")
    msg("Filtering coverage positions: normal depth >= ", depthmin.c.nor, " reads")

    coverage_rcmat_allcells <- coverage_raw_allcells %>%
        dplyr::filter(
            .data$TUM.DP >= depthmin.c.allcells,
            .data$NOR.DP >= depthmin.c.nor
        ) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)


    coverage_rcmat_allcells <- na.omit(coverage_rcmat_allcells)

    msg("Coverage positions kept: ", nrow(coverage_rcmat_allcells))

    # ==========================================================================
    # Recombine counts

    rcmat_allcells_filtered <- dplyr::bind_rows(allelic_rcmat_allcells, coverage_rcmat_allcells) %>%
        dplyr::arrange(.data$Chromosome, .data$Position, .data$cluster)

    msg("Final number of positions: ", nrow(rcmat_allcells_filtered))

    x@cnacalling[["combined.counts.allcells.filtered"]] <- rcmat_allcells_filtered


    # ==========================================================================
    # GET SEGMENTS AND ASSOCIATED DATA -----------------------------------------

    msg("Performing segmentation on all cells...")

    xx_allcells <- preProcSample2(
        rcmat_allcells_filtered,
        cval = cval1,
        snp.nbhd = snp.nbhd,
        gbuild = gbuild,
        ndepth = 1,
        ndepthmax = depthmax.nor)
    oo_allcells <- facets::procSample(xx_allcells,
                                      cval = cval2,
                                      min.nhet = min.nhet,
                                      dipLogR = NULL)

    # ==========================================================================
    # COMPUTE CF, TCN, LCN -----------------------------------------------------

    msg("Computing cell fractions and copy numbers on all cells...")

    # Compute CF Cell fraction and CN Copy Number with EM algorithm
    fit_allcells <- emcncf(oo_allcells, min.nhet=min.nhet)
    # Add CF and CN computed with EM to table
    out_allcells <- cbind(oo_allcells$out, fit_allcells$cncf[, c("start", "end", "cf.em", "tcn.em", "lcn.em")])
    jseg_allcells <- oo_allcells$jointseg

    # Compute median of vafT per segment
    jseg_allcells <- jseg_allcells %>%
        dplyr::group_by(.data$seg) %>%
        dplyr::mutate(vafT.median = median(.data$vafT[.data$signal == "allelic"], na.rm = TRUE))
    out_allcells <- dplyr::left_join(out_allcells, unique(jseg_allcells[, c("seg", "vafT.median")]), by="seg")

    jseg_allcells <- jseg_allcells[, colnames(jseg_allcells) != "vafT.median"]
    out_allcells <- out_allcells[, colnames(out_allcells)[c(1:6, 18, 7:17)]]

    # Rename chr X instead of 23
    out_allcells$chrom[out_allcells$chrom == 23] <- "X"
    jseg_allcells$chrom[jseg_allcells$chrom == 23] <- "X"

    # ==========================================================================
    # CF THRESHOLD -------------------------------------------------------------

    if (!is.null(cf.thresh)) {
        # Correct CN to diploid / neutral is CF below threshold
        out_allcells <- dplyr::mutate(
            out_allcells,
            tcn.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 2, .default = .data$tcn.em),
            lcn.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 1, .default = .data$lcn.em),
            cf.em = dplyr::case_when(!is.na(.data$cf.em) & .data$cf.em < cf.thresh ~ 1, .default = .data$cf.em))
    }

    # ==========================================================================
    # GET VARIANT ALLELE FREQUENCIES -------------------------------------------

    # Get vafT from allcells to set colors for variant positions depending on all cells frequency
    pmat_allcells_filtered <- procSnps(
        rcmat_allcells_filtered,
        snp.nbhd = snp.nbhd,
        nX = nX,
        ndepth = 1,
        ndepthmax = depthmax.nor)
    colnames(pmat_allcells_filtered)[7:8] <- c("cluster", "signal")  # Rename columns

    # Create new column colVAR depending on VAF in tumor cells
    pmat_allcells_filtered$colVAR <- rep(NA, nrow(pmat_allcells_filtered))
    pmat_allcells_filtered$colVAR[which(pmat_allcells_filtered$vafT >= 0.5 &
                                            pmat_allcells_filtered$signal == "allelic")] <- 1
    pmat_allcells_filtered$colVAR[which(pmat_allcells_filtered$vafT < 0.5 &
                                            pmat_allcells_filtered$signal == "allelic")] <- 2
    colnames(pmat_allcells_filtered)[colnames(pmat_allcells_filtered) == "vafT"] <- "vafT.allcells"
    # Rename chr X instead of 23
    pmat_allcells_filtered$chrom[pmat_allcells_filtered$chrom == 23] <- "X"

    # Add new columns to tables of data
    jseg_full <- dplyr::left_join(jseg_full,
                                  pmat_allcells_filtered[, c("chrom", "maploc", "signal", "vafT.allcells", "colVAR")],
                                  by = c("chrom", "maploc", "signal"))
    jseg_allcells <- dplyr::left_join(jseg_allcells,
                                      pmat_allcells_filtered[, c("chrom", "maploc", "signal", "vafT.allcells", "colVAR")],
                                      by = c("chrom", "maploc", "signal"))

    # ==========================================================================
    # SAVE OBJECTS -------------------------------------------------------------

    x@cnacalling[["positions"]] <- jseg_full
    x@cnacalling[["segments"]] <- out_full

    x@cnacalling[["positions.allcells"]] <- jseg_allcells
    x@cnacalling[["segments.allcells"]] <- out_allcells


    # ==========================================================================
    # 3. GET CONSENSUS SEGMENTS ------------------------------------------------
    # ==========================================================================

    msg("- Consensus segments accross clusters -")
    msg("Finding consensus segments...")

    consensus_segs <- getSegConsensus(out_full,
                                      ncells = ncells,
                                      dist.breakpoints = dist.breakpoints)

    # ==========================================================================
    # DEFINE PLOIDY ------------------------------------------------------------

    if(ploidy == "median") { ploidy <- round(median(out_full$tcn.em)) }
    if (ploidy == "auto") {
        # define ploidy based on the number of segments being 2:1 or 4:2
        diploid <- length(which(out_full$tcn.em == 2 & out_full$lcn.em == 1))
        tetraploid <- length(which(out_full$tcn.em == 4 & out_full$lcn.em == 2))
        if (diploid >= tetraploid) ploidy <- 2
        if (diploid < tetraploid) ploidy <- 4
    }

    # ==========================================================================
    # ANNOTATE SEGMENTS --------------------------------------------------------

    table <- annotateSegments(
        x = out_full,
        consensus_segs = consensus_segs,
        ncells,
        ploidy = ploidy,
        minoverlap = minoverlap,
        clonal.thresh = clonal.thresh
    )

    # Add CNA information to consensus segments
    consensus_segs <- dplyr::left_join(consensus_segs,
                                       table[, c("chrom", "start", "end", "cna", "cna_clonal")],
                                       by = c("chrom", "start", "end")) %>%
        dplyr::group_by(.data$chrom, .data$start) %>%
        dplyr::mutate(cna = any(.data$cna)) %>%
        dplyr::mutate(cna_clonal = any(.data$cna_clonal)) %>%
        unique()

    msg(paste(nrow(consensus_segs), "consensus segments identified,", sum(consensus_segs$cna), "CNA segments identified"))

    # ==========================================================================
    # SAVE OBJECTS -------------------------------------------------------------

    x@cnacalling[["consensus.segs"]] <- consensus_segs
    x@cnacalling[["table"]] <- table
    x@cnacalling[["ncells"]] <- ncells

    x@cnacalling[["dipLogR.clusters"]] <- dipLogR
    x@cnacalling[["dipLogR.allcells"]] <- fit_allcells$dipLogR
    x@cnacalling[["purity.clusters"]] <- fit_full$purity
    x@cnacalling[["purity.allcells"]] <- fit_allcells$purity
    x@cnacalling[["ploidy.clusters"]] <- fit_full$ploidy
    x@cnacalling[["ploidy.allcells"]] <- fit_allcells$ploidy
    x@cnacalling[["ploidy"]] <- ploidy

    msg("Done.")

    return(x)
}


#' Process read count matrix and segmentation
#'
#' Function adapted from [facets::preProcSample()] to process read count matrix
#' and generates a segmentation tree. The modifications from the original
#' function includes:
#' - Incorporation of `cluster` and `signal` columns into the final result.
#' - Change in the log ratio (`cnlr`) for allelic data, it is computed as the
#' log ratio (`cnlr`) mean between the previous and next coverage positions.
#'
#' @param rcmat A data frame with 8 required columns: `Chrom`, `Pos`, `NOR.DP`,
#'   `NOR.RD`, `TUM.DP`, `TUM.RD`, `cluster`, and `signal` (`data.frame`).
#' @param het.thresh VAF (Variant Allele Frequency) threshold to call variant positions
#'   heterozygous (`numeric`). Default: 0.25.
#' @param snp.nbhd Window size for selecting loci to reduce serial correlation
#'   (`numeric`). Default: 250.
#' @param cval Critical value for segmentation (`numeric`). Default: 25.
#' @param gbuild Genome build used for alignment. One of \code{"hg19"},
#'   \code{"hg38"}, or \code{"mm10"} (`character`). Default: \code{"hg19"}.
#' @param hetscale Logical value indicating whether log odds ratio (logOR)
#'   should be scaled to give more weight in the test statistics for
#'   segmentation and clustering (`logical`). Usually only 10% of snps are hets
#'   and hetscale gives the logOR contribution to T-square as 0.25/proportion of
#'   hets. Default: \code{TRUE}.
#' @param ndepth Minimum depth in normal reference to keep (`numeric`). Default:
#'   5.
#' @param ndepthmax Maximum normal coverage threshold for filtering loci
#'   (`numeric`). Default: 1000.
#'
#' @return A list containing:
#' \describe{
#'   \item{pmat}{Read counts and other elements of all the loci.}
#'   \item{gbuild}{Genome build used for the analysis.}
#'   \item{nX}{Chromosome number for X (e.g., 23 for human, 20 for mouse).}
#'   \item{clusters}{Unique clusters from the processed data.}
#'   \item{seg.tree}{Segmentation tree for each chromosome.}
#'   \item{jointseg}{Segmented variant positions data.}
#'   \item{hscl}{Scaling factor for logOR data.}
#' }
#' @details
#' The function processes variant positions data to generate a segmentation
#' tree. It uses \code{\link[facets]{procSnps}} to compute initial values,
#' adjusts the log ratio for allelic signals, and computes segmentation using
#' \code{\link[facets]{segsnps}}.
#'
#' SNPs (or variants) in a genome are not evenly spaced, and loci are sampled
#' within the specified window (\code{snp.nbhd}) to reduce serial correlation.
#'
#' @seealso \code{\link[facets]{preProcSample}}, \code{\link[facets]{procSnps}},
#'   \code{\link[facets]{segsnps}}
#'
#' @source This function is derived from the [facets::preProcSample()] function,
#' with the modifications to fit its use in muscadet.
#'
#' Seshan VE, Shen R (2021). _facets: Cellular Fraction and Copy Numbers from Tumor Sequencing_.
#' R package version 0.6.2, [https://github.com/mskcc/facets](https://github.com/mskcc/facets).
#'
#' @references
#' \describe{
#'   \item{[`facets`][facets::facets-package] package}{Shen R, Seshan VE. FACETS:
#'   allele-specific copy number and clonal heterogeneity analysis tool for
#'   high-throughput DNA sequencing. Nucleic Acids Res. 2016 Sep 19;44(16):e131.
#'   doi: [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
#'   27270079; PMCID: PMC5027494.}
#' }
#'
#' @include objects.R
#'
#' @import dplyr
#' @import facets
#' @import pctGCdata
#'
#' @export
#'
#' @examples
#' library(facets)
#'
#' # Load example muscadet object
#' # data("exdata_muscadet")
#'
#' counts <- exdata_muscadet$cnacalling$combined.counts
#' counts <- counts[complete.cases(counts), 1:8]
#' counts_clus <- counts[which(counts$cluster == 1),]
#' result <- preProcSample2(counts_clus)
#'
preProcSample2 <- function(
        rcmat,
        het.thresh = 0.25,
        snp.nbhd = 250,
        cval = 25,
        gbuild = "hg38",
        hetscale = TRUE,
        ndepth = 1,
        ndepthmax = 1000
) {

    # rcmat correct format
    rcmat <- as.data.frame(rcmat)

    cols <- c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD", "cluster", "signal")
    stopifnot(
        "The `rcmat` must contain the required columns in the correct order." =
            all(colnames(rcmat) == cols)
    )

    gbuild <- match.arg(gbuild, c("hg19", "hg38", "mm10"))

    # Export gbuild CG data from pctGC data correctly
    # NOTE:
    # The function pctGCdata::getGCpct() attempts to retrieve the GC content data
    # using `get(gbuild, pos = "package:pctGCdata")`:
    # https://github.com/veseshan/pctGCdata/blob/d2d4fafd10595750977f9bf5ebfeaed6e78da781/R/getGCpct.R#L5
    # This approach requires the 'pctGCdata' package to be attached (via library(pctGCdata)),
    # which is not standard practice and can cause errors when the package is only loaded via namespace.
    # To work around this, we manually retrieve the GC data object and pass it directly to `facets::counts2logROR()`
    # via the `ugcpct` argument while setting `gbuild = "udef"`.
    # This bypasses the environment-dependent `get()` call and avoids requiring
    # the 'pctGCdata' package to be attached globally.
    gbuild_obj <- getExportedValue("pctGCdata", gbuild)

    # Determine chromosome number for X based on genome build
    nX <- switch(
        gbuild,
        "hg19" = 23,
        "hg38" = 23,
        "mm10" = 20,
        stop("Unsupported genome: ", gbuild)
    )

    procSnps <- utils::getFromNamespace("procSnps", "facets")

    # Step 1: Process SNPs and retain cluster & signal columns
    pmat_1 <- procSnps(
        rcmat,
        ndepth = ndepth,
        het.thresh = het.thresh,
        snp.nbhd = snp.nbhd,
        nX = nX,
        unmatched = FALSE,
        ndepthmax = ndepthmax
    )
    colnames(pmat_1)[7:8] <- c("cluster", "signal")  # Rename columns

    # Step 2: Process SNPs with only required columns to get het & keep columns
    pmat_2 <- procSnps(
        rcmat[, 1:6],
        ndepth = ndepth,
        het.thresh = het.thresh,
        snp.nbhd = snp.nbhd,
        nX = nX,
        unmatched = FALSE,
        ndepthmax = ndepthmax
    )

    # Combine both results
    pmat <- cbind(pmat_1, pmat_2[, c("het", "keep")])

    # Step 3: Compute logR and logOR with GC correction for logR
    counts2logROR <- utils::getFromNamespace("counts2logROR", "facets")
    dmat <- counts2logROR(pmat[pmat$rCountT > 0, ], gbuild = "udef", ugcpct = gbuild_obj, unmatched = FALSE)

    # Step 4: Exclude log ratios for allelic data, replace by mean of neighbors
    dmat <- dmat %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$chrom) %>% # Group by chromosome
        dplyr::arrange(.data$maploc, .by_group = TRUE) %>% # Ensure rows are ordered by maploc
        dplyr::mutate(
            cnlr = dplyr::if_else(
                .data$signal == "allelic",
                # Compute mean of neighbors
                (dplyr::lag(.data$cnlr, default = 0) + dplyr::lead(.data$cnlr, default = 0)) /
                    (ifelse(is.na(lag(.data$cnlr, default = NA)), 0, 1) + ifelse(is.na(lead(.data$cnlr, default = NA)), 0, 1)),
                .data$cnlr # Keep existing values for coverage rows
            )
        ) %>%
        dplyr::ungroup()
    dmat <- as.data.frame(dmat)

    # Step 5: Segment SNPs
    segsnps <- utils::getFromNamespace("segsnps", "facets")
    seg_results <- segsnps(dmat,
                           cval = cval,
                           hetscale = hetscale,
                           delta = 0)

    # Make sure the joint segmentation data frame is ordered
    seg_results$jointseg <- seg_results$jointseg[order(seg_results$jointseg[, "maploc"]), ]
    seg_results$jointseg <- seg_results$jointseg[order(seg_results$jointseg[, "chrom"]), ]

    # Prepare output
    result <- list(
        pmat = pmat,
        gbuild = gbuild,
        nX = nX,
        clusters = unique(pmat$cluster)
    )
    if(!is.null(pmat$cluster)) result$clusters <- unique(pmat$cluster)

    # Combine output with segmentation results
    return(c(result, seg_results))
}


#' Get consensus segments across clusters
#'
#' This function processes segments data, in which clusters have different
#' segment breakpoints, to identify consensus segments across all clusters. It
#' groups breakpoints within a specified genomic distance (`dist.breakpoints`)
#' and calculates representative breakpoints for each group based on cluster
#' size and median coordinates.
#'
#' @param x A data frame containing cluster segments information with the
#'   following required columns:
#'   - `chrom`: Chromosome name (`factor` or `character`).
#'   - `start`: Start position of the segment (`numeric`).
#'   - `end`: End position of the segment (`numeric`).
#'   - `cluster`: Cluster identifier (`numeric`).
#' @param ncells A named vector specifying the number of cells per cluster
#'   (`numeric` vector). The names must match the cluster identifiers in the
#'   `cluster` column of `x`.
#' @param dist.breakpoints A numeric value specifying the minimum genomic
#'   distance between adjacent breakpoints to be grouped into the same consensus
#'   segment (`numeric` value). Default: `1e6`.
#'
#' @return A data frame containing consensus segments with the following columns:
#'   - `chrom`: Chromosome name.
#'   - `start`: Start position of the consensus segment.
#'   - `end`: End position of the consensus segment.
#'
#' @import dplyr
#' @importFrom GenomicRanges GRanges findOverlaps reduce
#' @importFrom rlang .data
#'
#' @export
#'
#' @seealso [muscadet::annotateSegments()], [muscadet::cnaCalling()]
#'
#' @examples
#' # Example data frame
#' segs <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2", "chr2"),
#'   start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
#'   end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
#'   cluster = c("1", "2", "1", "2")
#' )
#'
#' # Generate consensus segments
#' consensus_segs <- getSegConsensus(segs,
#'                                   ncells = c("1" = 50, "2" = 30),
#'                                   dist.breakpoints = 1e6)
#' print(consensus_segs)
#'
getSegConsensus <- function(x, ncells, dist.breakpoints = 1e6) {
    # Check if x is a data frame
    if (!is.data.frame(x)) {
        stop("'x' must be a data frame.")
    }

    # Check for required columns in x
    required_columns <- c("chrom", "start", "end", "cluster")
    missing_columns <- setdiff(required_columns, colnames(x))
    if (length(missing_columns) > 0) {
        stop(
            "The data frame 'x' is missing the following required columns: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    # Check if dist.breakpoints is numeric and positive
    if (!is.numeric(dist.breakpoints) ||
        length(dist.breakpoints) != 1 || dist.breakpoints <= 0) {
        stop("'dist.breakpoints' must be a single positive numeric value.")
    }

    # Check for missing values in required columns
    if (anyNA(x[, required_columns])) {
        stop("The data frame 'x' contains missing values in required columns.")
    }

    # Get start and end breakpoints
    x$chrom <- factor(x$chrom, levels = unique(x$chrom))
    breakpoints_start <- dplyr::arrange(x, .data$chrom, .data$start)[, c("chrom", "start", "start", "cluster")]
    colnames(breakpoints_start)[c(2, 3)] <- c("start", "end")
    breakpoints_end <- dplyr::arrange(x, .data$chrom, .data$end)[, c("chrom", "end", "end", "cluster")]
    colnames(breakpoints_end)[c(2, 3)] <- c("start", "end")

    # Merge all breakpoints in one GRanges object
    breakpoints <- GenomicRanges::GRanges(arrange(
        rbind(breakpoints_start, breakpoints_end),
        .data$chrom,
        .data$start,
        .data$end
    ))

    # Find clusters of breakpoints (named "group" here to distinguish from the clusters of cells)
    hits <- GenomicRanges::findOverlaps(
        breakpoints,
        GenomicRanges::reduce(breakpoints, min.gapwidth = dist.breakpoints)
    )
    breakpoints <- as.data.frame(breakpoints)
    breakpoints$group <- hits@to

    # Order clusters of cells from the smallest to the largest
    breakpoints$cluster <- factor(breakpoints$cluster,
                                  levels = names(ncells)[order(ncells)],
                                  ordered = T)

    # In case of unique breakpoint in a chromosome (unique segment < dist.breakpoints):
    # separation into 2 distinct new groups
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$seqnames, .data$cluster) %>%
        dplyr::mutate(ngroupperchr = length(unique(.data$group))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$seqnames, .data$start) %>%
        dplyr::mutate(group = dplyr::case_when(
            .data$ngroupperchr == 1 ~ max(breakpoints$group) + cur_group_id(),
            .default = .data$group
        )) %>%
        dplyr::ungroup()

    # If there are several breakpoints in one group of breakpoints per cluster:
    # compute median to get a unique breakpoint
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$group, .data$cluster) %>%
        dplyr::mutate(start = round(median(.data$start)), end = round(median(.data$end))) %>%
        dplyr::ungroup() %>%
        unique()

    # In each group of breakpoints: select the breakpoint coordinates from the largest cluster of cells
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$group) %>%
        dplyr::mutate(start = .data$start[.data$cluster == max(.data$cluster)],
                      end = .data$end[.data$cluster == max(.data$cluster)]) %>%
        dplyr::select(-"cluster", -"ngroupperchr") %>%
        unique()

    # Transform breakpoints back into segments
    colnames(breakpoints)[1] <- "chrom"
    consensus_segs <- breakpoints %>%
        dplyr::group_by(.data$chrom) %>%
        dplyr::mutate(start = if (n() == 1) {
            .data$start # Single breakpoint: start remains as is
        } else {
            c(.data$start[1:(length(.data$start) - 1)], NA) # Multiple breakpoints: create segments
        }, end = if (n() == 1) {
            .data$end # Single breakpoint: end remains as is
        } else {
            c(.data$end[2:length(.data$end)], NA) # Multiple breakpoints: create segments
        }) %>%
        dplyr::select(-"group", -"width", -"strand") %>%
        dplyr::filter(!is.na(.data$start) & !is.na(.data$end)) # Remove rows with NA start or end

    return(as.data.frame(consensus_segs))
}



#' Annotate consensus segments across cluster with CNA states
#'
#' This function annotates consensus segments across clusters using
#' segment-level data from each cluster. It assigns states (i.e., gain, loss,
#' cnloh, neutral), computes proportions of cells per cluster, and defines
#' whether alterations are clonal.
#'
#' @param x A data frame containing cluster segments information with the
#'   following required columns:
#'   - `chrom`: Chromosome name (`factor` or `character`).
#'   - `start`: Start position of the segment (`numeric`).
#'   - `end`: End position of the segment (`numeric`).
#'   - `cluster`: Cluster identifier (`numeric` or `integer`).
#'   - `tcn.em`: Total copy number (EM algorithm) (`numeric`).
#'   - `lcn.em`: Lower copy number (EM algorithm) (`numeric`).
#' @param consensus_segs A `data.frame` of unique consensus segments.
#' @param ploidy An `integer` specifying the expected ploidy (default is `2`).
#'   Used to determine whether a segment is classified as gain, loss, or neutral.
#' @param minoverlap A non-negative integer specifying the minimum overlap
#'   between a cluster specific segment and a consensus segment (see
#'   [IRanges::findOverlaps()]) (`integer`). Default: `1e6`.
#' @param clonal.thresh Minimum cell proportion to label a
#'   segment as clonal. Default: `1e6`.
#'
#' @inheritParams getSegConsensus
#'
#' @return A `data.frame` with the consensus segments annotated per cluster,
#'   including CNA state, clonal status, and proportion of cells per cluster and
#'   across clusters:
#'   - `ncells`: Number of cells in cluster.
#'   - `prop.cluster`: Proportion of cells in cluster relatively to the total number of cells.
#'   - `gnl`: GNL (gain ; neutral ; loss) status as an integer value (1 ; 0 ; -1).
#'   - `loh`: Loss of heterozygosity status (logical).
#'   - `state`: State of segment (gain ; loss ; neu ; cnloh).
#'   - `cna`: Segment CNA status (logical).
#'   - `cna_state`: State of CNA segment (gain ; loss ; cnloh).
#'   - `prop.tot`: Proportion of cells for the state across clusters.
#'   - `state_clonal`: Clonal state of segment (gain ; loss ; neu ; cnloh).
#'   - `cna_clonal`: Segment clonal CNA status (logical).
#'   - `cna_clonal_state`: Clonal state of CNA segment (gain ; loss ; cnloh).
#'
#' @seealso [getSegConsensus()], [cnaCalling()].
#'
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' #' # Example data frame
#' segs <- data.frame(
#'     chrom = c("chr1", "chr1", "chr2", "chr2"),
#'     start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
#'     end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
#'     cluster = c("1", "2", "1", "2"),
#'     cf.em = c(1, 1, 1, 1),
#'     tcn.em = c(1, 2, 3, 3),
#'     lcn.em = c(0, 1, 1, 1)
#' )
#'
#' consensus_segs <- getSegConsensus(
#'     x = segs,
#'     ncells = c("1" = 50, "2" = 30)
#' )
#'
#' # Annotate consensus segments
#' table <- annotateSegments(
#'     x = segs,
#'     consensus_segs = consensus_segs,
#'     ncells = c("1" = 50, "2" = 30)
#' )
#'
#' print(table)
#'
annotateSegments <- function(
        x,
        consensus_segs,
        ncells,
        ploidy = 2,
        minoverlap = 1e6,
        clonal.thresh = 0.9
) {

    segs <- x[, c("chrom", "start", "end", "cluster", "cf.em", "tcn.em", "lcn.em")]

    # Annotate consensus segments data for each cluster
    out.segs <- Map(function(cluster) {

        segs_clus <- segs[segs$cluster == cluster, ]
        segsGR <- GenomicRanges::GRanges(segs_clus[, c("chrom", "start", "end")])
        consensGR <- GenomicRanges::GRanges(consensus_segs[, c("chrom", "start", "end")])

        idx_segs <- GenomicRanges::findOverlaps(
            consensGR,
            segsGR,
            minoverlap = minoverlap,
            select = "first"
        )
        valid_idx <- !is.na(idx_segs)

        cbind(consensus_segs[valid_idx, ],
              data.frame(id = which(valid_idx)),
              segs_clus[idx_segs[valid_idx], -c(1:3)])
    },
    unique(sort(segs$cluster))
    )

    out.segs <- Reduce(rbind, out.segs)

    # Add number of cells and cluster proportions
    out.segs$ncells <- ncells[as.character(out.segs$cluster)]
    prop.ncells <- ncells / sum(ncells)
    out.segs$prop.cluster <- prop.ncells[as.character(out.segs$cluster)]

    # Add gain-neutral-loss (gnl) status, state gain-neu-loss-cnloh
    out.segs <- dplyr::mutate(
        out.segs,
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
        cna = dplyr::case_when(
            .data$state == "neu" ~ FALSE,
            .data$state %in% c("gain", "loss", "cnloh") ~ TRUE,
            .default = NA
        ),
        cna_state = ifelse(.data$cna, .data$state, NA)
    )

    # Compute proportions across all clusters with same state
    # group by segment id and by state (to separate the same segment if this one has different states across clusters)
    out.segs <- out.segs %>%
        dplyr::group_by(.data$id, .data$state) %>%
        dplyr::mutate(prop.tot = sum(.data$prop.cluster)) %>%
        dplyr::ungroup()


    # Define clonal status
    out.segs <- dplyr::mutate(
        out.segs,
        state_clonal = dplyr::case_when(.data$prop.tot >= clonal.thresh ~ .data$state),
        cna_clonal = dplyr::case_when(
            .data$cna == TRUE & .data$prop.tot >= clonal.thresh ~ TRUE,
            .data$cna == TRUE & .data$prop.tot < clonal.thresh ~ FALSE,
            .data$cna == FALSE ~ FALSE
        ),
        cna_clonal_state = dplyr::case_when(.data$cna_clonal == TRUE ~ .data$state)
    )

    return(as.data.frame(out.segs))
}


