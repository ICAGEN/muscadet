#' Assign a cluster assignment to a `muscadet` object
#'
#' Add the user selected cluster assignments to cells in a
#' \code{\link{muscadet}} object. This function allows the user to choose the
#' cluster assignments they consider to fit the data and their requirements, or
#' cluster assignments based on other data and methods.
#'
#' @param x A \code{\link{muscadet}} object (`muscadet`).
#'
#' @param partition Value specifying the clustering partition to choose from the
#'   muscadet object (`numeric` or `character`). It should be either the resolution or the
#'   k number of cluster (k) used for clustering depending on the clustering
#'   method (`res_range` or `k_range` with [muscadet::clusterMuscadet()]).
#'   Should be provided if `clusters` is `NULL`.
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
#' - The clusters can be taken directly from the `muscadet` object clustering
#' results with setting the `parition` argument (e.g.
#' `muscadet_obj$clustering$clusters[["0.8"]]` for res=`0.8`).
#' - A custom vector of cluster assignments
#' can be attributed using the `clusters` argument.
#' - Either way, the clusters assignments can be rearranged using the `mapping`
#' argument.
#'
#' @return A \code{\link{muscadet}} object updated with the user chosen cluster
#' assignments in `muscadet_obj$cnacalling$clusters`.
#'
#' @include objects.R
#' @importFrom SeuratObject Cells
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' data(muscadet_obj)
#'
#' # Select clustering result for partition = 0.6
#' muscadet_obj <- assignClusters(muscadet_obj, partition = 0.6)
#' table(muscadet_obj$cnacalling$clusters)
#'
#' # Assign custom clusters
#' set.seed(42)
#' cell_names <- Reduce(union, SeuratObject::Cells(muscadet_obj))
#' n1 <- sample(1:length(cell_names), 1)
#' n2 <- length(cell_names) - n1
#' custom_clusters <- setNames(c(rep.int(1, n1), rep.int(2, n2)), cell_names)
#' table(custom_clusters)
#' muscadet_obj <- assignClusters(muscadet_obj, clusters = custom_clusters)
#' table(muscadet_obj$cnacalling$clusters)
#'
#' # Assign clusters with remapping
#' # example to remap from partition=0.8 with merging of clusters 2 and 3
#' clusters <- muscadet_obj$clustering$clusters[["0.8"]]
#' table(clusters) # 3 clusters
#' mapping <- c("1" = 1, "2" = 2, "3" = 2) # remap to 2 clusters
#'
#' muscadet_obj <- assignClusters(muscadet_obj, clusters = clusters, mapping = mapping)
#' table(muscadet_obj$cnacalling$clusters)
#' # check original and remapped clusters
#' table(clusters, muscadet_obj$cnacalling$clusters)
#'
#' muscadet_obj <- assignClusters(muscadet_obj, partition = 0.8, mapping = mapping)
#' table(muscadet_obj$cnacalling$clusters)
#' # check original and remapped clusters
#' table(muscadet_obj$clustering$clusters[["0.8"]],
#'       muscadet_obj$cnacalling$clusters)
#'
#' # Visualize clusters on heatmap
#' heatmapMuscadet(
#'     muscadet_obj,
#'     partition = 0.8,
#'     filename = file.path("heatmap_muscadet_res0.8.png"),
#'     title = "Example sample | res=0.8"
#' )
#' heatmapMuscadet(
#'     muscadet_obj,
#'     clusters = muscadet_obj$cnacalling$clusters,
#'     filename = file.path("heatmap_muscadet_custom_res0.8.png"),
#'     title = "Example sample | rearranged clusters from res=0.8"
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
        "Either `partition` or `clusters` must be provided, but not both." = xor(!is.null(partition), !is.null(clusters))
    )

    # If partition is provided, validate that clustering has been performed
    if (!is.null(partition)) {
        stopifnot(
            "Clustering results are not available in the muscadet object `x`. Perform clustering first using `clusterMuscadet()`."
            = !is.null(x@clustering),
            "Clustering results for the chosen `partition` are not available. Use `clusterMuscadet()` with the `res_range`/`k_range` argument to compute the partition."
            = as.character(partition) %in% names(x@clustering$clusters)
        )
    }

    # If clusters is provided, validate its format
    if (!is.null(clusters)) {
        if (is.factor(clusters) && !is.null(names(clusters))) {
            clusters <- stats::setNames(as.vector(clusters), names(clusters))
        } else {
            stopifnot("`clusters` must be a named vector or factor." = is.vector(clusters) && !is.null(names(clusters)))
        }
    }

    # Apply mapping if provided
    if (!is.null(mapping)) {
        # Validate mapping
        stopifnot(
            "`mapping` must be a named vector." = is.vector(mapping) &&
                !is.null(names(mapping)),
            "All cluster values must have corresponding mappings in `mapping`." = all(as.character(unique(clusters)) %in% names(mapping))
        )

        # Get clusters from partition
        if (!is.null(partition)) {
            clusters <- x@clustering$clusters[[as.character(partition)]]
        }

        # Remap clusters
        remapped_clusters <- setNames(mapping[as.character(clusters)], names(clusters))

        # Rerun cluster imputation (following remapping)
        if (redo_imputation) {

            mat_list <- lapply(muscadet::matLogRatio(x), t)
            common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames)))

            # Validate clusters cells
            stopifnot(
                "`clusters` must contain at least all common cells when `redo_imputation = TRUE`." = all(common_cells %in% names(remapped_clusters))
            )

            # Remove cells missing data that needs to go through cluster imputation (cleaner)
            remapped_clusters <- remapped_clusters[intersect(names(remapped_clusters), common_cells)]
            # Impute clusters
            remapped_clusters <- imputeClusters(mat_list, remapped_clusters, knn_imp = knn_imp)

        } else {
            # Validate clusters cells
            stopifnot(
                "`clusters` must contain all cell names within the muscadet object `x`." = all(names(remapped_clusters) %in% Reduce(union, Cells(x)))
            )
        }


        x@cnacalling[["clusters"]] <- remapped_clusters

    } else {
        # Assign clusters without remapping - custom clusters
        if (!is.null(clusters)) {

            # Rerun cluster imputation (in case custom clusters have been modified)
            if (redo_imputation) {

                mat_list <- lapply(muscadet::matLogRatio(x), t)
                common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames)))

                # Validate clusters cells
                stopifnot(
                    "`clusters` must contain at least all common cells when `redo_imputation = TRUE`." = all(common_cells %in% names(clusters))
                )

                # Remove cells missing data that needs to go through cluster imputation (cleaner)
                clusters <- clusters[intersect(names(clusters), common_cells)]
                # Impute clusters
                clusters <- imputeClusters(mat_list, clusters, knn_imp = knn_imp)

            } else {
                # Validate clusters cells
                stopifnot(
                    "`clusters` must contain all cell names within the muscadet object `x`." = all(names(clusters) %in% Reduce(union, Cells(x)))
                )
            }
            x@cnacalling[["clusters"]] <- clusters

        } else if (!is.null(partition)) {
            # Assign clusters without remapping - clusters from partition
            x@cnacalling[["clusters"]] <- x@clustering$clusters[[as.character(partition)]]
        }
    }

    return(x)
}


#' Add allele counts to a `muscadet` object
#'
#' This function adds allele counts data to the `omics` of a
#' \code{\link{muscadet}} object. The data frames in the `allele_counts` list
#' are assigned to the `allelic` slots of `omics`.
#'
#' @param x A \code{\link{muscadet}} object.
#' @param allele_counts A list of data frames where each data frame contains
#'   allele counts for a specific omic (`list`). The list must have the same
#'   length and order as the number of omics in the `x` object. Each data frames
#'   must contain the following columns : `cell`, `id`, `CHROM`, `POS`, `REF`,
#'   `ALT`, `RD`, `AD`, `DP`, `GT`. See [allele_counts] for details.
#'
#' @return
#' A modified \code{\link{muscadet}} object with updated allele counts in the
#' `allelic` slot of each `muscomic` object in the `omics` slot.
#'
#' @note
#' As the allele counts are not used for computing log R
#' ratios and cell clustering, they are not mandatory at the creation of
#' `muscomic` and `muscadet` objects. The allele counts data can be added to
#' objects later with this `addAlleleCounts` function, before using the
#' [muscadet::mergeCounts()] function.
#'
#' This function is also useful to add allele counts for individual-specific
#' variant positions to a common `muscadet` object, for example for the
#' reference `muscadet` object: a common `muscadet` object with reference
#' coverage data can be stored as a unique object to use as reference for
#' computing log R ratios for all samples, and then it can updated with allele
#' counts at individual-specific variant positions (e.g. found by bulk
#' sequencing) before Copy Number Alterations (CNAs) calling.
#'
#' @seealso [muscadet::CreateMuscomicObject()], [muscadet::mergeCounts()]
#'
#' @importFrom stringr str_remove
#' @importFrom gtools mixedsort
#'
#' @export
#'
#' @examples
#' # Load example muscadet objects
#' data(muscadet_obj)
#' data(muscadet_obj_ref)
#'
#' # Add allele counts data frames to muscadet objects
#' muscadet_obj <- addAlleleCounts(
#'     muscadet_obj,
#'     allele_counts = list(allele_counts_atac_tumor, allele_counts_rna_tumor))
#' muscadet_obj_ref <- addAlleleCounts(
#'     muscadet_obj_ref,
#'     allele_counts = list(allele_counts_atac_ref, allele_counts_rna_ref))
#'
addAlleleCounts <- function(x, allele_counts) {

    # Check that x is a muscadet object
    stopifnot("The argument 'x' must be a muscadet object." = inherits(x, "muscadet"))

    # Check that allele_counts is a list of dataframes
    stopifnot("The argument 'allele_counts' must be a list." = inherits(allele_counts, "list"))
    stopifnot("The argument 'allele_counts' must be a list of data frames." = all(unlist(
        lapply(allele_counts, function(x)
            inherits(x, "data.frame"))
    )))

    # Ensure allele_counts list matches the number of omics in x
    if (length(allele_counts) != length(x@omics)) {
        stop("The 'allele_counts' list must have the same length as the number of omics in 'x'.")
    }

    # Assign names to allele_counts based on omics names in x
    names(allele_counts) <- names(x@omics)
    # Extract types from omics in x
    type.omic <- unlist(lapply(x@omics, function(x) x$type))

    # Format allele count data frames
    allele_counts <- lapply(allele_counts, function(allele_df) {
        # Remove "chr" if necessary
        allele_df$CHROM <- stringr::str_remove(allele_df$CHROM, "chr")
        # Ordered chromosomes
        allele_df$CHROM <- ordered(allele_df$CHROM, levels = gtools::mixedsort(unique(allele_df$CHROM)))
        # Reorder data
        allele_df <- allele_df[order(allele_df[, "POS"]), ]
        allele_df <- allele_df[order(allele_df[, "CHROM"]), ]
        return(allele_df)

    })

    # Add allele counts to each omic in the muscadet object
    for (i in names(x@omics)) {
        slot(x@omics[[i]], "allelic") <- list(table.counts = data.frame(omic = type.omic[[i]], allele_counts[[i]]))
    }

    return(x)
}


#' Merge counts for `muscadet` objects
#'
#' This function combines allelic (counts at variant positions, either common
#' SNPs or individual-specific heterozygous positions) and coverage counts
#' (counts on features) from all omics per cluster for both sample and
#' reference. The resulting merged data is stored in the `cnacalling` slot of
#' the sample `muscadet` object.
#'
#' @param x A \code{\link{muscadet}} object containing sample data (`muscadet`).
#'   This object must include clustering assignments in the
#'   `cnacalling$clusters` slot.
#' @param reference A \code{\link{muscadet}} object containing reference data
#'   (`muscadet`).
#' @param nor.het A logical value to specify if normal reference allele counts
#'   are modified to: total normal depth counts divided by 2, to force these
#'   positions to be heterozygous in the normal reference in allelic data (e.g.
#'   when heterozygous positions are retrieve based on matched bulk sequencing
#'   data, and are thereby assumed to be heterozygous) before combining coverage
#'   and allelic data. Default is `TRUE`.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return
#' A modified \code{\link{muscadet}} object corresponding to the `x` muscadet object,
#' with updated `cnacalling` slot containing:
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
#' - omic = omic specific (`omic` column)
#' - all = for all omics
#'
#' @import dplyr
#' @importFrom data.table rbindlist
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' # Load example muscadet objects
#' data(muscadet_obj)
#' data(muscadet_obj_ref)
#'
#' # Merge counts from all omics from both sample and reference
#' muscadet_obj <- mergeCounts(muscadet_obj, muscadet_obj_ref)
#'
mergeCounts <- function(x,
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
            !is.null(omic@allelic$table.counts)
        })))
    )
    stopifnot(
        "Input object 'reference' must contain allele counts data in the 'allelic' slot of omics. Use 'addAlleleCounts()' to add them."
        = all(unlist(lapply(reference@omics, function(omic) {
            !is.null(omic@allelic$table.counts)
        })))
    )

    # Allelic Data Processing --------------------------------------------------
    if (!quiet) {
        message("Allelic data processing...")
    }
    # Extract allelic counts for all omics in both sample and reference
    allele_sample <- data.table::rbindlist(lapply(x@omics, function(omic) {
        omic@allelic$table.counts
    }), use.names = TRUE)
    allele_ref <- data.table::rbindlist(lapply(reference@omics, function(omic) {
        omic@allelic$table.counts
    }), use.names = TRUE)


    # Compute per-omic and total counts for the sample
    # Filter for cells with valid cluster assignments
    allele_sample <- allele_sample[allele_sample$cell %in% names(x@cnacalling$clusters),]
    # Add cluster assignments to table
    allele_sample$cluster <- x@cnacalling$clusters[match(allele_sample$cell, names(x@cnacalling$clusters))]
    # Compute per-omi counts
    summed_sample_omic <- allele_sample %>%
        dplyr::group_by(.data$cluster, .data$id, .data$omic)%>%
        dplyr::summarise(
            RD.omic = sum(.data$RD, na.rm = TRUE),
            AD.omic = sum(.data$AD, na.rm = TRUE),
            DP.omic = sum(.data$DP, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(AF.omic = round(.data$AD.omic / .data$DP.omic, 2))
    # Compute total counts
    summed_sample_all <- allele_sample %>%
        dplyr::group_by(.data$cluster, .data$id)%>%
        dplyr::summarise(
            RD.all = sum(.data$RD, na.rm = T),
            AD.all = sum(.data$AD, na.rm = T),
            DP.all = sum(.data$DP, na.rm = T),
            .groups = "drop"
        ) %>%
        dplyr::mutate(AF.all = round(.data$AD.all / .data$DP.all, 2))
    # Join data
    allele_sample_full <- allele_sample %>%
        dplyr::select(!c("cell", "RD", "AD", "DP")) %>%
        dplyr::left_join(.data, summed_sample_omic, by = c("cluster", "id", "omic")) %>%
        dplyr::left_join(.data, summed_sample_all, by = c("cluster", "id")) %>%
        subset(.data, !duplicated(.data))

    # Compute per-omic and total counts for the reference
    # Compute per-omi counts
    summed_ref_omic <- allele_ref %>%
        dplyr::group_by(.data$id, .data$omic)%>%
        dplyr::summarise(
            RD.omic = sum(.data$RD, na.rm = TRUE),
            AD.omic = sum(.data$AD, na.rm = TRUE),
            DP.omic = sum(.data$DP, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(AF.omic = round(.data$AD.omic / .data$DP.omic, 2))
    # Compute total counts
    summed_ref_all <- allele_ref %>%
        dplyr::group_by(.data$id)%>%
        dplyr::summarise(
            RD.all = sum(.data$RD, na.rm = T),
            AD.all = sum(.data$AD, na.rm = T),
            DP.all = sum(.data$DP, na.rm = T),
            .groups = "drop"
        ) %>%
        dplyr::mutate(AF.all = round(.data$AD.all / .data$DP.all, 2))
    # Join data
    allele_ref_full <- allele_ref %>%
        dplyr::select(!c("cell", "RD", "AD", "DP")) %>%
        dplyr::left_join(.data, summed_ref_omic, by = c("id", "omic")) %>%
        dplyr::left_join(.data, summed_ref_all, by = c("id")) %>%
        subset(.data, !duplicated(.data))

    # Merge sample and reference allelic data
    var <- dplyr::left_join(
        allele_sample_full,
        allele_ref_full,
        by = c("CHROM", "POS", "REF", "ALT", "id", "GT", "omic"),
        suffix = c(".TUM", ".NOR")
    ) %>%
        dplyr::arrange(.data$CHROM, .data$POS, .data$cluster) %>%
        dplyr::mutate(signal = "allelic")

    x@cnacalling[["allelic.counts"]] <- as.data.frame(var)


    # Auto set normal variant positions as heterozygous before combining
    if(nor.het == TRUE) {
        var[, "RD.all.NOR"] <- round(var[, "DP.all.NOR"] / 2, 0)
    }


    # Coverage Data Processing -------------------------------------------------
    if (!quiet) {
        message("Coverage data processing...")
    }
    # Extract coverage counts for all omics in both sample and reference
    coverage_sample <- data.table::rbindlist(lapply(x@omics, function(omic) {
        omic@coverage$table.counts
    }), use.names = TRUE)
    coverage_ref <- data.table::rbindlist(lapply(reference@omics, function(omic) {
        omic@coverage$table.counts
    }), use.names = TRUE)

    # Compute total counts for the sample
    # Filter for cells with valid cluster assignments
    coverage_sample <- coverage_sample[coverage_sample$cell %in% names(x@cnacalling$clusters),]
    # Add cluster assignments to table
    coverage_sample$cluster <- x@cnacalling$clusters[match(coverage_sample$cell, names(x@cnacalling$clusters))]
    summed_sample <- coverage_sample %>%
        dplyr::group_by(.data$cluster, .data$id) %>%
        dplyr::summarise(DP = sum(.data$DP),
                         .groups = "drop")
    coverage_sample_full <- coverage_sample %>%
        dplyr::select(!c("cell", "cluster", "DP")) %>%
        dplyr::left_join(.data, summed_sample, by = "id") %>%
        subset(.data, !duplicated(.data))

    # Compute total counts for the reference
    summed_ref <- coverage_ref %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarise(DP = sum(.data$DP),
                         .groups = "drop")
    coverage_ref_full <- coverage_ref %>%
        dplyr::select(!c("cell", "DP")) %>%
        dplyr::left_join(.data, summed_ref, by = "id") %>%
        subset(.data, !duplicated(.data))

    # Merge sample and reference coverage data
    cov <- dplyr::left_join(
        coverage_sample_full,
        coverage_ref_full,
        by = c("CHROM", "POS", "id", "omic"),
        suffix = c(".TUM", ".NOR")
    ) %>%
        dplyr::arrange(.data$CHROM, .data$POS, .data$cluster) %>%
        dplyr::select("omic", "id", "CHROM", "POS", "DP.NOR", "DP.TUM", "cluster") %>%
        dplyr::mutate(signal = "coverage")

    x@cnacalling[["coverage.counts"]] <- as.data.frame(cov)


    # Combine Allelic and Coverage Data ----------------------------------------
    if (!quiet) {
        message("Combining allelic and coverage data...")
    }

    # Format data
    var <- var[, c("CHROM", "POS", "DP.all.NOR", "RD.all.NOR", "DP.all.TUM", "RD.all.TUM", "cluster", "signal", "omic", "id")]
    var <- unique(var)
    cov <- cov[, c("CHROM", "POS", "DP.NOR", "DP.NOR", "DP.TUM", "DP.TUM", "cluster", "signal", "omic", "id")]
    cov <- unique(cov)

    colnames(var) <- c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD", "cluster", "signal", "omic", "id")
    colnames(cov) <- colnames(var)

    # Make sure the levels match for binding
    levels(var$Chromosome) <- union(levels(var$Chromosome), levels(cov$Chromosome))
    levels(cov$Chromosome) <- union(levels(var$Chromosome), levels(cov$Chromosome))

    # Combine both data
    combined <- rbind(var, cov) %>%
        arrange(.data$Chromosome, .data$Position)

    x@cnacalling[["combined.counts"]] <- as.data.frame(combined)

    return(x)
}


#' Save a data frame as a VCF File
#'
#' This function saves a given data frame in [variant calling format
#' (VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) with necessary
#' headers.
#'
#' @param data A data frame containing variant data (`data.frame`). The first
#'   column should represent the chromosome.
#' @param file A character string specifying the output file path with the
#'   `.vcf` or `.vcf.gz` extension (`character string`). If the file has the
#'   `.vcf.gz` extension it will be compressed.
#' @param header An optional custom header for the VCF format (`character`
#'   vector). By default, `NULL` sets a minimal header.
#'
#' @details
#' The function ensures that the first column is named `#CHROM` and writes standard VCF headers.
#' It then appends the variant data in a tab-separated format.
#'
#' @return The function does not return a value but writes a VCF file to the specified path.
#'
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' vcf_data <- data.frame(
#'   CHROM = c("chr1", "chr2"),
#'   POS = c(12345, 67890),
#'   ID = c(".", "."),
#'   REF = c("A", "T"),
#'   ALT = c("G", "C"),
#'   QUAL = c(".", "."),
#'   FILTER = c("PASS", "PASS"),
#'   INFO = c("AF=0.5", "AF=0.3"),
#'   FORMAT = c("GT", "GT"),
#'   SAMPLE = c("0/1", "1/1")
#' )
#' save_as_vcf(vcf_data, "output.vcf")
#' }
#'
#' @export
save_as_vcf <- function(data, file, header = NULL) {
    # Check if input is a data frame
    stopifnot("The `data` must be a data frame." = is.data.frame(data))

    # Check if file extension
    stopifnot("The output file must have a .vcf or .vcf.gz extension." = grepl("\\.vcf(\\.gz)?$", file, ignore.case = TRUE))

    # Ensure the first column is named "#CHROM" (VCF standard)
    colnames(data)[1] <- "#CHROM"

    # VCF headers lines
    if (is.null(header)) {
        header <- c(
            "##fileformat=VCFv4.1",
            "##INFO=<ID=AF,Number=.,Type=Integer>",
            "##FORMAT=<ID=GT,Number=1,Type=String>"
        )
    }

    # Determine if compression is needed
    compress <- grepl("\\.vcf.gz$", file, ignore.case = TRUE)

    # If compressing, write to a temporary .vcf file first
    if (compress) {
        temp_vcf <- tempfile(fileext = ".vcf")
    } else {
        temp_vcf <- file
    }

    # Write VCF headers
    writeLines(header, temp_vcf)
    # Append variant data in tab-separated format without quotes
    suppressWarnings(
        utils::write.table(
            data,
            temp_vcf,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            row.names = FALSE,
            col.names = TRUE
        )
    )

    # Compress if needed
    if (compress) {
        if (!requireNamespace("Rsamtools", quietly = TRUE)) {
            stop("The 'Rsamtools' package is required for bgzip compression and tabix indexing. Please install it.")
        }

        # Compress with bgzip
        Rsamtools::bgzip(temp_vcf, dest = file, overwrite = TRUE)
        unlink(temp_vcf)  # Remove temporary file

        # Index the VCF with tabix
        Rsamtools::indexTabix(file, format = "vcf")
    }
}



#' Process allele counts results from SCReadCounts format
#'
#' This function reformats SCReadCounts output into a standardized VCF-like data
#' frame and optionally merges unphased/phased genotype information from a
#' VCF-format data frame.
#'
#' @param data A data frame of SCReadCounts output table (`data.frame`).
#' @param vcf An optional VCF-format data frame with unphased/phased genotype
#'   information to add it in output (`data.frame`).
#' @param samplename An optional sample name to prefix cell barcodes
#'   (`character`).
#'
#' @return A cleaned data frame with VCF-like standardized columns:
#' \describe{
#'   \item{`cell`}{Barcodes of cells (`character`).}
#'   \item{`id`}{Variant unique identifier defined as
#'   CHROM_POS_REF_ALT, e.g. "1_920949_C_G" (`character`).}
#'   \item{`CHROM`}{Chromosome in integer format, e.g. 15 (X and Y chromosomes
#'   are not included) (`integer`).}
#'   \item{`POS`}{Position of the variant (1-base positions) (`integer`).}
#'   \item{`REF`}{Reference allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`ALT`}{Alternative allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`RD`}{Reference allele depth/count (`integer`).}
#'   \item{`AD`}{Alternative allele depth/count (`integer`).}
#'   \item{`DP`}{Total depth/count (`integer`).}
#'   \item{`GT`}{Genotype (if `vcf` provided): "0/1" or "1/0" if unphased; "0|1"
#'   or "1|0" if phased. (`character`).}
#' }
#'
#' @note
#' - [SCReadCounts website](https://horvathlab.github.io/NGS/SCReadCounts/)
#' - [Variant Call Format (VCF) format](https://en.wikipedia.org/wiki/Variant_Call_Format)
#'
#' @references
#' [SCReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/)
#'
#' Liu, H., Bousounis, P., Movassagh, M., Edwards, N., and Horvath, A.
#' SCReadCounts: estimation of cell-level SNVs expression from scRNA-seq data.
#' BMC Genomics 22, 689 (2021).
#' doi: [10.1186/s12864-021-07974-8](https://doi.org/10.1186/s12864-021-07974-8)
#'
#' @examples
#'
#' # Example data
#' sc_data <- data.frame(
#'   ReadGroup = c("cell1", "cell2"),
#'   CHROM = c(1, 1),
#'   POS = c(1001, 1002),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T"),
#'   SNVCount = c(3, 5),
#'   RefCount = c(7, 5),
#'   GoodReads = c(10, 10)
#' )
#'
#' vcf_data <- data.frame(
#'   CHROM = c(1, 1),
#'   POS = c(1001, 1002),
#'   ID = c(".", "."),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T"),
#'   QUAL = c(".", "."),
#'   FILTER = c(".", "."),
#'   INFO = c(".", "."),
#'   FORMAT = c("GT", "GT"),
#'   sample1 = c("0|1", "1|0"),
#'   stringsAsFactors = FALSE
#' )
#'
#' data <- process_allele(sc_data, vcf = vcf_data, samplename = "sampleA")
#'
#' @export
process_allele <- function(data, vcf = NULL, samplename = NULL) {

    # Check for inputs
    stopifnot("The `data` must be a data frame." = is.data.frame(data))

    required_cols <- c("CHROM", "POS", "REF", "ALT", "ReadGroup", "SNVCount", "RefCount", "GoodReads")
    missing <- setdiff(required_cols, colnames(data))
    if (length(missing) > 0) stop("Missing required columns in `data`: ", paste(missing, collapse = ", "))

    # Create ID column
    data$id <- paste(data$CHROM, data$POS, data$REF, data$ALT, sep = "_")

    if (!is.null(vcf)) {
        stopifnot("The `vcf` must be a data frame." = is.data.frame(vcf))
        # Rename the column 10 to standard (instead of the sample name)
        colnames(vcf)[10] <- "GT"
        # Create ID column
        vcf$id <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = "_")
        # Add genotype information (possibly phased genotype)
        df <- left_join(data, vcf[, c("id", "GT")], by = "id")
    } else {
        df <- data
        df$GT <- "" # temporary column
    }

    # Chromosomes as factor
    df$CHROM <- factor(df$CHROM, sort(unique(df$CHROM)))

    # Rename columns
    colnames(df)[colnames(df) == "ReadGroup"] <- "cell"
    colnames(df)[colnames(df) == "SNVCount"] <- "AD"
    colnames(df)[colnames(df) == "RefCount"] <- "RD"
    colnames(df)[colnames(df) == "GoodReads"] <- "DP"

    # Select final columns
    df <- df[, c("cell", "id", "CHROM", "POS", "REF", "ALT", "RD", "AD", "DP", "GT")]
    # Remove temporary GT column if vcf not provided
    if (is.null(vcf)) {
        df$GT <- NULL
    }

    # Order by position
    df <- df[order(df$CHROM, df$POS), ]

    # Add sample name to the cell barcode if samplename provided
    if (!is.null(samplename)) {
        df$cell <- paste(samplename, df$cell, sep = "_")
    }

    # Make sure unique data only remain
    df <- unique(df)

    return(df)
}
