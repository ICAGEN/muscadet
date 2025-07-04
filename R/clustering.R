#' Multi Omics Integration and Clustering on a `muscadet` Object
#'
#' Performs integration of multi omics and clustering of cells based on log
#' ratio data contained in a [`muscadet`][muscadet-class] object.
#'
#' Two methods are available for integration and clustering of common cells
#' between omics:
#' - Method `seurat` uses nearest neighbors for integration followed by
#' graph-based clustering.
#' - Method `hclust` uses Similarity Network Fusion (SNF) for integration
#' followed by hierarchical clustering.
#'
#' Then, clusters are imputed for cells missing data in at least one omic, by
#' similarity using nearest neighbor cells.
#'
#' Finally, silhouette widths are computed on the integrated distance matrix to
#' help identify the optimal clustering partition.
#'
#' @param x A [`muscadet`][muscadet-class] object containing omics data and
#'   previously computed log R ratio matrices (`muscadet`).
#' @param method The clustering method to apply (`character` string). One of
#'   `"seurat"` or `"hclust"`. For  medthod `"seurat"`, arguments for
#'   [cluster_seurat()] should be provided, and for method `"hclust"`, arguments
#'   for [cluster_hclust()] should be provided. Note: The `"seurat"` method can
#'   only be applied for a maximum of 2 omics. Default is `"seurat"`.
#' @param omics Optional character vector specifying omic names to use for
#'   clustering. Must match names of available omics in the `x` muscadet object.
#'   If `NULL` (default), all available omics are used.
#' @param knn_imp Number of k nearest neighbors to use for imputing cluster
#'   assignments of cells missing in one or more omics. Only relevant for more
#'   than one omic.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @inheritDotParams cluster_seurat res_range dims_list algorithm knn_seurat knn_range_seurat
#' @inheritDotParams cluster_hclust k_range dist_method hclust_method weights
#'
#'
#' @return The input [`muscadet`][muscadet-class] object with its `clustering`
#'   slot updated. This slot contains:
#' \describe{
#'   \item{params}{List of parameters used for clustering (`list`).}
#'   \item{...}{Output objects depending on the method. See [cluster_seurat()] or [cluster_hclust()].}
#'   \item{clusters}{A named list of cluster partitions (named vectors of cluster
#'   labels) for all cells (imputed clusters assignments for non-common cells),
#'   for each value in `k_range` or `res_range` (`list`).}
#'   \item{silhouette}{A list of silhouette objects and widths for each cluster partition (`list`).}
#'   \item{partition.opt}{Name of the optimal cluster partition based on maximum
#'   average silhouette width.}
#' }
#'
#' @seealso
#' Methodology and functionality:
#'
#' - [muscadet-class]
#' - [cluster_seurat()] for graph-based clustering using [`Seurat`][Seurat::Seurat-package].
#' - [cluster_hclust()] for hierarchical clustering of SNF-fused distances.
#' - [weightedSNF()] for weighted Similarity Network Fusion (SNF).
#' - [imputeClusters()] for imputing cluster labels across omics.
#'
#' Visualization:
#' - [heatmapMuscadet()] to plot clustering result as heatmap.
#' - [plotSil()] to plot silhouette widths.
#' - [plotIndexes()] to plot several normalized cluster validation indexes.
#'
#' Select clusters to continue with CNA calling:
#' - [assignClusters()] to assign final cluster assignments in the `muscadet`
#' object after cluster partition validation.
#'
#' @include objects.R
#'
#' @importFrom methods slot slot<-
#' @importFrom cluster silhouette
#' @importFrom stats as.dist
#'
#' @export
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Perform clustering with "seurat" method
#' muscadet_obj <- clusterMuscadet(
#'   x = muscadet_obj,
#'   method = "seurat",
#'   res_range = c(0.5, 0.8),
#'   dims_list = list(1:8, 1:8),
#'   knn_seurat = 10, # adapted to low number of cells in example data
#'   knn_range_seurat = 30 # adapted to low number of cells in example data
#' )
#'
#' # Perform clustering with "hclust" method
#' muscadet_obj <- clusterMuscadet(
#'   x = muscadet_obj,
#'   k_range = 2:4,
#'   method = "hclust",
#'   dist_method = "euclidean",
#'   hclust_method = "ward.D",
#'   weights = c(1, 1)
#' )
#'
#' # Retrieve cluster assignments
#' clusters <- muscadet_obj$clustering$clusters
#' lapply(clusters, table)
#'
clusterMuscadet <- function(x,
                            method = c("seurat", "hclust"),
                            omics = NULL,
                            knn_imp = 10,
                            quiet = FALSE,
                            ...) {

    # Validate inputs
    stopifnot("Input object `x` must be of class 'muscadet'." = inherits(x, "muscadet"))
    method <- match.arg(method)

    if (is.null(omics)) {
        omics <- names(x$omics)
    } else {
        stopifnot(
            "Some specified `omics` are not present in the muscadet object `x`." =
                all(omics %in% names(x$omics))
        )
    }

    # Empty clustering slot
    if(length(slot(x, "clustering")) > 0) {
        slot(x, "clustering") <- list()
    }
    # Empty cna calling slot to avoid inconsistent results
    if(length(slot(x, "cnacalling")) > 0) {
        slot(x, "cnacalling") <- list()
    }

    # Check that log ratio matrices exists
    stopifnot(
        "At least one log ratio matrix is missing for one omic." =
            !any(unlist(lapply(
                matLogRatio(x), is.null
            )))
    )
    # Extract and transpose log ratio matrices for each omic
    mat_list <- lapply(muscadet::matLogRatio(x), t)

    # Clustering ---------------------------------------------------------------

    # Method 1: Using Seurat functions (1 or 2 omics only)
    if(method == "seurat") {

        if (!quiet) {
            message("Clustering method: 'seurat'")
        }

        # Clustering
        res_clust <- cluster_seurat(mat_list[omics], ...)
    }

    # Method 2: Using hclust
    if(method == "hclust") {

        if (!quiet) {
            message("Clustering method: 'hclust'")
        }

        # Clustering
        res_clust <- cluster_hclust(mat_list[omics], ...)
    }

    slot(x, "clustering") <- res_clust

    params <- c(list(omics = omics), res_clust$params, list(knn_imp = knn_imp))
    slot(x, "clustering")[["params"]] <- params

    # Clusters imputation ------------------------------------------------------

    # Retrieve list of clusters
    clusters <- res_clust$clusters

    if (!quiet) {
        message("Imputing clusters...")
    }

    # Impute clusters assignments to missing cells based on nearest neighbors (if several omics)
    if (length(mat_list) > 1) {
        clusters_imp <- lapply(clusters, function(clus) {
            imputeClusters(mat_list, clus, knn_imp = knn_imp)
        })
        slot(x, "clustering")[["clusters"]] <- clusters_imp
    } else {
        slot(x, "clustering")[["clusters"]] <- clusters
    }

    # Find optimal partition ---------------------------------------------------

    if (!quiet) {
        message("Computing Silhouette scores...")
    }

    # Retrieve (integrated) distance matrix
    dist <- stats::as.dist(res_clust$dist)

    # Compute Silhouette
    partitions_list <- setNames(as.list(names(clusters)), names(clusters))
    sil <- lapply(partitions_list, function(partition) {
        cl <- as.vector(clusters[[as.character(partition)]])
        if (length(unique(cl)) > 1) {
            cluster::silhouette(as.vector(clusters[[as.character(partition)]]), dist)
        }
    })

    # Extract Silhouette widths
    sil_widths <- lapply(sil, function(x)
        if (!is.null(x)) {
            summary(x)$avg.width
        })
    sil_widths_cl <- lapply(sil, function(x)
        if (!is.null(x)) {
            summary(x)$clus.avg.widths
        })

    # Find optimal partition based on Silhouette
    sil_widths_vec <- unlist(sil_widths)
    if (!is.null(sil_widths_vec)) {
        partition_opt <- names(sil_widths_vec[sil_widths_vec == max(sil_widths_vec)])
    } else {
        partition_opt <- NULL
    }

    slot(x, "clustering")[["silhouette"]] <- list(sil.obj = sil,
                                                  sil.w.avg = sil_widths,
                                                  sil.w.avg.clusters = sil_widths_cl)

    slot(x, "clustering")[["partition.opt"]] <- partition_opt

    if (!quiet) {
        message("Done.")
    }

    # Return the updated muscadet object
    return(x)
}

#' Multi Omics Clustering using Seurat Multi Modal Graph-based Clustering
#'
#' Performs graph-based clustering of cells using
#' [`Seurat`][Seurat::Seurat-package], based on one or two log R ratio matrices
#' (`mat_list`), including shared nearest neighbors (SNN) graph construction on
#' selected dimensions from PCA (`dims_list`), to identify clusters of cells for
#' each specified resolution (`res_range`).
#'
#' - For two omics: multimodal integration is performed using
#' [Seurat::FindMultiModalNeighbors()] (weighted shared nearest neighbors graph).
#' Only common cells between omics are used.
#' - For a single omic: [Seurat::FindNeighbors()] (shared nearest neighbors
#' graph) is used.
#'
#' @param mat_list A named list of log R ratio matrices (cells x features), one
#'   per omic layer (`list`).
#' @param res_range A numeric non-negative vector specifying the resolution
#'   values to use for [Seurat::FindClusters()] (`numeric` vector). Default is
#'   `c(0.1, 0.2, 0.3, 0.4, 0.5)`.
#' @param dims_list A list of vectors of PC dimensions to use for each omic
#'   (`list`). Must match the length of `mat_list` (e.g., list(1:8) for 1 omic ;
#'   list(1:8, 1:8) for 2 omics). Default is the first 8 dimensions for each
#'   provided omic.
#' @param algorithm Integer specifying the algorithm for modularity optimization
#'   by [Seurat::FindClusters()] (`1` = original Louvain algorithm; `2` = Louvain
#'   algorithm with multilevel refinement; `3` = SLM algorithm; `4` = Leiden
#'   algorithm). Leiden requires the leidenalg python. Default is `1`.
#' @param knn_seurat Integer specifying the number of nearest neighbors used for
#'   graph construction with [Seurat::Seurat()] functions [Seurat::FindNeighbors()]
#'   (`k.param`) or [Seurat::FindMultiModalNeighbors()] (`k.nn`) (`integer`).
#'   Default is `20`.
#' @param knn_range_seurat Integer specifying the approximate number of nearest
#'   neighbors to compute for [Seurat::FindMultiModalNeighbors()] (`knn.range`)
#'   (`integer`). Default is `200`.
#' @param max_dim Integer specifying the maximum number of principal components
#'   to be used for PCA computation with [stats::prcomp()] (`integer`). Default
#'   is `200`.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{params}{List of parameters used for clustering (`list`).}
#'   \item{pcs}{List of principal components summaries for each omic (`list` of
#'   [`summary.prcomp`][stats::summary.prcomp]).}
#'   \item{nn}{Nearest neighbors object ([`Neighbor`][SeuratObject::Neighbor-class]).}
#'   \item{graph}{Shared nearest neighbors graph ([`Graph`][SeuratObject::Graph-class]).}
#'   \item{dist}{Distance matrix derived from the graph (`matrix`).}
#'   \item{umap}{UMAP coordinates (`matrix`).}
#'   \item{clusters}{A named list of clustering results (vectors of cluster
#'   labels) for each value in `res_range` (`list`).}
#' }
#'
#' @seealso [Weighted Nearest Neighbor Analysis Vignette from
#' Seurat](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)
#'
#' @include objects.R
#'
#' @importFrom SeuratObject CreateAssay5Object as.sparse
#' @importFrom Seurat CreateSeuratObject FindMultiModalNeighbors FindNeighbors RunUMAP FindClusters
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Format input
#' # transpose matrices to: cells x features matrices
#' mat_list <- lapply(muscadet::matLogRatio(muscadet_obj), t)
#'
#' # Run integration & clustering
#' result <- cluster_seurat(mat_list, res_range = c(0.1, 0.3, 0.5))
#'
#' # View results
#' lapply(result$clusters, table)
#' }
#'
cluster_seurat <- function(mat_list,
                           res_range = seq(0.1, 0.5, 0.1),
                           dims_list = rep(list(1:8), length(mat_list)),
                           algorithm = 1,
                           knn_seurat = 20,
                           knn_range_seurat = 200,
                           max_dim = 200,
                           quiet = FALSE
) {
    # Arguments validations
    stopifnot(is.list(mat_list), all(sapply(mat_list, is.matrix)))

    stopifnot(
        "This method is not applicable for more than 2 omics. `mat_list` can contain either 1 or 2 matrices."
        = length(mat_list) <= 2
    )

    if (!is.numeric(res_range) || any(res_range <= 0)) {
        stop("`res_range` must contain only positive numeric values.")
    }

    stopifnot("Length of `dims_list` must match the number of omics in `mat_list`." = length(mat_list) == length(dims_list))

    if (!is.numeric(knn_seurat) || knn_seurat < 2 || knn_seurat != as.integer(knn_seurat)) {
        stop("`knn_seurat` must be an integer >= 2.")
    }
    if (!is.numeric(knn_range_seurat) || knn_range_seurat < 2 || knn_range_seurat != as.integer(knn_range_seurat)) {
        stop("`knn_range_seurat` must be an integer >= 2.")
    }

    # Initialize output
    out <- list()

    # Store clustering parameters
    out[["params"]] <- list(
        method = "seurat",
        res_range = res_range,
        dims_list = dims_list,
        knn_seurat = knn_seurat,
        knn_range_seurat = knn_range_seurat
    )

    # Get common cell barcodes across all omics
    common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames)))

    # Check that number of neighbors does not exceed number of cells
    if (knn_seurat >= length(common_cells)) {
        knn_seurat <- length(common_cells) - 1
        warning(
            paste0(
                "`knn_seurat` exceeds the number of cells (",
                length(common_cells),
                ") and is therefore set to ",
                knn_seurat,
                "."
            )
        )
    }
    if (knn_range_seurat >= length(common_cells)) {
        knn_range_seurat <- length(common_cells) - 1
        warning(
            paste0(
                "`knn_range_seurat` exceeds the number of cells (",
                length(common_cells),
                ") and is therefore set to ",
                knn_range_seurat,
                "."
            )
        )
    }

    # Filter matrices on common cells
    mat_list <- lapply(mat_list, function(mat) {
        mat[common_cells, ]
    })

    # Create assay objects for each modality, use only common cells
    # We transpose the matrix since Seurat assays expect features x cells matrices
    assay_list <- lapply(mat_list, function(mat) {
        SeuratObject::CreateAssay5Object(SeuratObject::as.sparse(t(mat)),
                                         # rows = features x columns = cells
                                         min.cells = 0,
                                         min.features = 0)
    })

    # Join assays into a Seurat object
    seurat <- Seurat::CreateSeuratObject(assay_list[[1]], assay = names(assay_list)[1])
    for (i in seq_along(assay_list)[-1]) {
        seurat[[names(assay_list[i])]] <- assay_list[[i]]
    }

    # Perform PCA individually on log ratio (common cells only)
    if (!quiet) {
        message("Performing PCA...")
    }
    pcs_list <- lapply(mat_list, function(mat) {
        prcomp(x = mat, rank. = max_dim, center = TRUE, scale. = FALSE) # rows = cells x columns = features
    })
    pcs_list_sum <- lapply(pcs_list, summary)

    out[["pcs"]] <- pcs_list_sum

    # Create DimReducObject for Seurat and add to seuratobj
    pca_list <- lapply(setNames(as.list(names(pcs_list)), names(pcs_list)), function(omic) {
        SeuratObject::CreateDimReducObject(
            embeddings = pcs_list[[omic]]$x, # cells × PCs
            loadings = pcs_list[[omic]]$rotation, # features × PCs
            stdev = pcs_list[[omic]]$sdev,
            key = "PC_",
            assay = omic
        )
    })

    # Add reductions to Seurat object
    # Warning: since Seurat v5 the way to add reductions has changed, see:
    # https://github.com/satijalab/seurat/discussions/8726#discussioncomment-10330615
    for (omic in names(pca_list)) {
        seurat@reductions[[paste0("PCA_", omic)]] <- pca_list[[omic]]
    }

    # Find nearest neighbors & run UMAP

    # Case for 2 omics: find multi modal neighbors
    if (length(seurat@reductions) > 1) {
        # Generate SNN graph for UMAP and clustering
        if (!quiet) {
            message("Finding neighbors and constructing graph...")
        }
        suppressWarnings(seurat <- Seurat::FindMultiModalNeighbors(
            seurat,
            reduction.list = as.list(names(seurat@reductions)),
            dims.list = dims_list,
            k.nn = knn_seurat,
            knn.range = knn_range_seurat,  # must be lower than the number of cells
            verbose = FALSE
        ))

        # Run UMAP on the nearest neighbor graph
        if (!quiet) {
            message("Computing UMAP...")
        }
        options(Seurat.warn.umap.uwot = FALSE)
        seurat <- Seurat::RunUMAP(
            seurat,
            nn.name = "weighted.nn",
            reduction.name = "nn.umap",
            reduction.key = "nnUMAP_",
            verbose = FALSE
        )
        # Case for 1 omic: find neighbors
    } else {
        # Generate SNN graph for UMAP and clustering
        if (!quiet) {
            message("Finding neighbors and constructing graph...")
        }
        seurat <- Seurat::FindNeighbors(
            seurat,
            reduction = grep("PCA", names(seurat@reductions), value = TRUE),
            dims = unlist(dims_list),
            k.param = knn_seurat,
            compute.SNN = TRUE,
            graph.name = c("nn", "snn"), # compute graph
            verbose = FALSE
        )
        seurat <- Seurat::FindNeighbors(
            seurat,
            reduction = grep("PCA", names(seurat@reductions), value = TRUE),
            dims = unlist(dims_list),
            k.param = knn_seurat,
            return.neighbor = TRUE, # compute neighbors
            verbose = FALSE
        )

        # Run UMAP on the nearest neighbor graph
        if (!quiet) {
            message("Computing UMAP...")
        }
        seurat <- Seurat::RunUMAP(
            seurat,
            nn.name = names(seurat@neighbors),
            reduction.name = "nn.umap",
            reduction.key = "nnUMAP_",
            verbose = FALSE
        )
    }

    out[["nn"]] <- seurat@neighbors[[names(seurat@neighbors)]]

    # Convert graph into distance
    graph <- seurat@graphs[[grep("snn", names(seurat@graphs))]]
    tmp <- as.matrix(graph)
    diag(tmp) <- 0
    dist <- max(tmp) - as.matrix(graph)

    out[["graph"]] <- graph
    out[["dist"]] <- dist

    dim_reduc <- seurat@reductions[[grep("umap", names(seurat@reductions))]]
    umap <- dim_reduc@cell.embeddings
    colnames(umap) <- c("UMAP_1", "UMAP_2")
    out[["umap"]] <- umap

    # Clustering using the graph with different resolution
    if (!quiet) {
        message("Finding clusters...")
    }
    res_list <- setNames(as.list(res_range), res_range)
    clusters <- lapply(res_list, function(res) {
        seurat <- Seurat::FindClusters(
            seurat,
            algorithm = algorithm,
            resolution = res,
            graph.name = grep("snn", names(seurat@graphs), value = TRUE),
            verbose = FALSE
        )
        cl <- seurat@meta.data$seurat_clusters
        levels(cl) <- seq(1, length(levels(cl)), by = 1)
        cl <- setNames(as.integer(cl), Cells(seurat))
    })

    out[["clusters"]] <- clusters

    return(out)
}

#' Multi Omics Clustering with SNF integration and Hierarchical Clustering
#'
#' Performs the integration of log R ratio matrices (`mat_list`) using
#' Similarity Network Fusion (SNF) followed by hierarchical clustering, on the
#' integrated SNF matrix, to identify clusters of cells for each specified k
#' number of cluster (`k_range`). For more than 1 omic, the integration and
#' clustering are performed only on common cells between omics.
#'
#' @param mat_list A named list of log R ratio matrices (cells x features), one
#'   per omic layer (`list`).
#' @param k_range A numeric vector of integers (≥2) specifying the cluster
#'   numbers (k) to extract from hierarchical clustering (`numeric` vector).
#'   Default is from 2 to 10.
#' @param dist_method A string specifying the distance method for
#'   [Rfast::Dist()] (e.g., `"euclidean"`, `"manhattan"`, `"cosine"`)
#'   (`character` string). Default is `"euclidean"`.
#' @param hclust_method A string specifying the hierarchical clustering linkage
#'   method for [fastcluster::hclust()] (e.g., `"ward.D"`, `"average"`)
#'   (`character` string). Default is `"ward.D"`.
#' @param weights A numeric vector of non-negative values of length equal to the
#'   number of omic (internally normalized to sum to 1) (`numeric` vector). It
#'   specifies the relatives weights of each omic for SNF with
#'   [weightedSNF()]. Omics with a weight of 0 will not contribute to
#'   the clustering. If `NULL` (default), weights are uniform.
#' @param knn_affinity Integer specifying the number of nearest neighbors used
#'   when building affinity matrices with [SNFtool::affinityMatrix()]
#'   (`integer`). Default is `40`.
#' @param var_affinity Numeric value for the variance parameter (Gaussian kernel
#'   width `sigma`) when building affinity matrix with
#'   [SNFtool::affinityMatrix()] (`numeric`). Default is `1`.
#' @param knn_SNF Integer specifying the number of nearest neighbors used during
#'   the Similarity Network Fusion (SNF) with [weightedSNF()]
#'   (`integer`). Default is `40`.
#' @param iter_SNF Integer specifying the number of iterations for SNF with
#'   [weightedSNF()] (`integer`). Default is `50`.
#' @param knn_umap Integer specifying the number of nearest neighbors used for
#'   manifold approximation (UMAP) (see [uwot::umap()]). Default is `20`.
#' @param quiet Logical. If `TRUE`, suppresses informative messages during
#'   execution. Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{params}{List of parameters used for clustering (`list`).}
#'   \item{SNF}{Fused similarity matrix computed with SNF (`matrix`).}
#'   \item{dist}{Distance matrix derived from the SNF similarity (`matrix`).}
#'   \item{hclust}{Hierarchical clustering object from [fastcluster::hclust()] (`hclust`).}
#'   \item{umap}{UMAP coordinates (`matrix`).}
#'   \item{clusters}{A named list of clustering results (vectors of cluster
#'   labels) for each value in `k_range` (`list`).}
#' }
#'
#' @details
#' The function calculates pairwise distances between cells within each omic dataset
#' using the specified `dist_method` (using only common cells between omics).
#'
#' It constructs affinity matrices based on these distances, applies SNF to
#' generate a fused similarity matrix.
#'
#' `weights` can be assigned to each omic dataset to prioritize certain data types
#' over others, allowing users to tailor the analysis based on the
#' characteristics and importance of each dataset.
#'
#' It then performs a hierarchical clustering using the specified
#' `hclust_method` to assign clusters to each common cells.
#'
#' Results are given as cluster assignments for each number of cluster specified
#' by `k_range`.
#'
#' @seealso
#' Similarity Network Fusion: [weightedSNF()].
#'
#' @include objects.R
#'
#' @importFrom Rfast Dist
#' @importFrom SNFtool affinityMatrix
#' @importFrom fastcluster hclust
#' @importFrom stats as.dist
#' @importFrom dendextend cutree
#' @importFrom utils packageVersion
#' @importFrom uwot umap
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Format input
#' # transpose matrices to: cells x features matrices
#' mat_list <- lapply(muscadet::matLogRatio(muscadet_obj), t)
#'
#' # Run integration & clustering
#' result <- cluster_hclust(mat_list, k_range = 2:4)
#'
#' # View results
#' lapply(result$clusters, table)
#' }
#'
cluster_hclust <- function(mat_list,
                           k_range = seq(2, 10, 1),
                           dist_method = "euclidean",
                           hclust_method = "ward.D",
                           weights = rep(1, length(mat_list)),
                           knn_affinity = 40,
                           var_affinity = 1,
                           knn_SNF = 40,
                           iter_SNF = 50,
                           knn_umap = 20,
                           quiet = FALSE) {

    # Arguments validations
    stopifnot("`mat_list` must be a list of matrices." = is.list(mat_list) & all(sapply(mat_list, is.matrix)))
    stopifnot("Length of `weights` must match the number of omics in `mat_list`." = length(weights) == length(mat_list))
    if (!is.numeric(k_range) || any(k_range < 1) || any(k_range != as.integer(k_range))) {
        stop("`k_range` must be a numeric vector of integers greater than or equal to 2.")
    }

    # Initialize output
    out <- list()

    # Store clustering parameters
    out[["params"]] <- list(
        method = "hclust",
        k_range = k_range,
        dist_method = dist_method,
        hclust_method = hclust_method,
        weights = weights,
        knn_affinity = knn_affinity,
        var_affinity = var_affinity,
        knn_SNF = knn_SNF,
        iter_SNF = iter_SNF
    )

    # Get common cell barcodes across all omics
    common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames)))

    # Compute pairwise distances for each omic
    if (!quiet) {
        message("Computing distance matrices...")
    }
    dist_list <- lapply(mat_list, function(mat) {
        mat <- mat[common_cells, ]
        dist <- Rfast::Dist(mat, method = dist_method)
        # for version <=2.1, Rfast returns a similarity matrix for cosine with 0
        # values in diagonal. The result should be converted into distance.
        # See https://github.com/RfastOfficial/Rfast/issues/119 for more info
        if (dist_method == "cosine" &
            utils::packageVersion("Rfast") <= "2.1.0") {
            diag(dist) <- 1
            dist <- 1 - dist
        }
        dimnames(dist) <- list(rownames(mat), rownames(mat))
        return(dist)
    })

    # Compute affinity matrices for each omic
    if (!quiet) {
        message("Computing affinity matrices...")
    }
    aff_list <- lapply(dist_list, function(dist) {
        aff <- SNFtool::affinityMatrix(dist, K = knn_affinity, sigma = var_affinity)
        aff <- aff[sort(rownames(aff)), sort(colnames(aff))]
        return(aff)
    })

    # Perform Similarity Network Fusion (SNF)
    if (!quiet) {
        message("Performing SNF integration...")
    }
    matSNF <- weightedSNF(aff_list,
                          K = knn_SNF,
                          t = iter_SNF,
                          weights = weights)
    out[["SNF"]] <- matSNF

    # Convert the fused affinity matrix to a distance matrix
    tmp <- matSNF
    diag(tmp) <- 0
    dist <- max(tmp) - matSNF
    out[["dist"]] <- dist

    # Compute UMAP
    if (!quiet) {
        message("Computing UMAP...")
    }

    umap <- uwot::umap(
        X = matSNF,
        n_neighbors = knn_umap,
        metric = dist_method,
        verbose = FALSE
    )
    rownames(umap) <- rownames(matSNF)
    colnames(umap) <- c("UMAP_1", "UMAP_2")
    out[["umap"]] <- umap

    # Perform hierarchical clustering
    if (!quiet) {
        message("Performing hierarchical clustering...")
    }
    hc <- fastcluster::hclust(stats::as.dist(dist), hclust_method)
    out[["hclust"]] <- hc

    # Cut the dendrogram to generate clusters for the specified k range
    k_list <- setNames(as.list(k_range), k_range)
    clusters <- lapply(k_list, function(k) {
        dendextend::cutree(hc, k, order_clusters_as_data = FALSE)
    })

    out[["clusters"]] <- clusters

    return(out)
}


#' Weighted Similarity Network Fusion
#'
#' @description
#' Similarity Network Fusion (SNF) takes multiple views of a network and fuses
#' them together to construct a unified similarity matrix. Each affinity matrix
#' can be assigned a weight to control its relative contribution to the final
#' fused matrix. This approach enables combining information from multiple data
#' types or networks, preserving complementary structures across views.
#'
#'
#' @param Wall List of affinity matrices (`list`). Each element of the list is a
#'   square, symmetric matrix that shows affinities of the data points from a
#'   certain view. If only one matrix is provided, the function returns the
#'   unique input matrix.
#' @param K Number of neighbors in K-nearest neighbors (`integer`). Default is
#'   `20`.
#' @param t Number of iterations for the diffusion process (`integer`). Default
#'   is `20`.
#' @param weights Numeric vector of non-negative values of length equal to
#'   `length(Wall)`, specifying the relative weights for each affinity matrix in
#'   the fusion process (internally normalized to sum to 1) (`numeric` vector).
#'   Matrices with a weight of 0 are excluded from the fusion. If only one
#'   matrix has a non-zero weight, the function returns that matrix unchanged.
#'   Defaults to equal weights for all matrices.
#'
#'
#' @return
#' An overall status matrix derived, unified similarity graph of all data
#' types. It contains both complementary information and common
#' structures from all individual network.
#'
#' If only one matrix is provided (`Wall`), the function returns the unique input matrix.
#'
#' @export
#'
#' @source This function is derived from the [SNFtool::SNF()] function, with the
#' addition of the `weights` argument.
#'
#' Wang B, Mezlini A, Demir F, Fiume M, Tu Z, Brudno M, Haibe-Kains B,
#' Goldenberg A (2021). _SNFtool: Similarity Network Fusion_. R package version
#' 2.3.1, [https://CRAN.R-project.org/package=SNFtool](https://CRAN.R-project.org/package=SNFtool).
#'
#' @references
#' Wang B, Mezlini AM, Demir F, Fiume M, Tu Z, Brudno M, Haibe-Kains B,
#' Goldenberg A. Similarity network fusion for aggregating data types on a
#' genomic scale. Nat Methods. 2014 Mar;11(3):333-7. doi: [10.1038/nmeth.2810](http://doi.org/10.1038/nmeth.2810).
#'
#' Concise description of SNF can be found here:
#' [http://compbio.cs.toronto.edu/SNF/SNF/Software.html](http://compbio.cs.toronto.edu/SNF/SNF/Software.html)
#'
#' @examples
#' set.seed(123) # For reproducibility
#'
#' # Number of samples
#' n <- 20
#'
#' # Generate the first affinity matrix: block diagonal structure
#' block1 <- matrix(rnorm(n * n / 4, mean = 1, sd = 0.1), n / 2, n / 2)
#' block2 <- matrix(rnorm(n * n / 4, mean = 1, sd = 0.1), n / 2, n / 2)
#' A1 <- rbind(cbind(block1, matrix(0, n / 2, n / 2)), cbind(matrix(0, n / 2, n / 2), block2))
#' A1 <- A1 + diag(n) # Add self-similarity (diagonal dominance)
#'
#' # Generate the second affinity matrix: random similarity matrix
#' A2 <- matrix(runif(n * n, min = 0, max = 1), n, n)
#' A2 <- (A2 + t(A2)) / 2 # Make it symmetric
#' A2 <- A2 + diag(n) # Add self-similarity (diagonal dominance)
#'
#' # Normalize rows to sum to 1 (transition matrices)
#' A1 <- A1 / rowSums(A1)
#' A2 <- A2 / rowSums(A2)
#'
#' # Create a list of affinity matrices
#' affinity_list <- list(A1, A2)
#'
#' # Visualize the matrices
#' par(mfrow = c(3, 2))
#'
#' matSNF1 <- weightedSNF(affinity_list, K = 10, weights = c(1, 0))
#' matSNF2 <- weightedSNF(affinity_list, K = 10, weights = c(0, 1))
#' matSNF3 <- weightedSNF(affinity_list, K = 10, weights = c(1, 1))
#' matSNF4 <- weightedSNF(affinity_list, K = 10, weights = c(0.8, 0.2))
#'
#' # Set common breaks to get the same color scale for all
#' breaks <- quantile(c(A1, A2, matSNF1, matSNF2, matSNF3, matSNF4), seq(0, 1, 0.1))
#'
#' image(A1, main = "Matrix 1 (Block Structure)", col = heat.colors(10), breaks = breaks)
#' image(A2, main = "Matrix 2 (Random)", col = heat.colors(10), breaks = breaks)
#'
#' image(matSNF1, main = "Weight only Matrix", col = heat.colors(10), breaks = breaks)
#' image(matSNF2, main = "Weight only Matrix", col = heat.colors(10), breaks = breaks)
#'
#' image(matSNF3, main = "Equal weights", col = heat.colors(10), breaks = breaks)
#' image(matSNF4, main = "Respectively 0.8 and 0.2 weights", col = heat.colors(10), breaks = breaks)
#'
weightedSNF <- function(
        Wall,
        K=20,
        t=20,
        weights = rep.int(1, length(Wall)))
{
    # Similarity Network Fusion takes multiple views of a network (Wall) and
    # fuses them together to create a overall affinity matrix.

    # Check list of matrices input
    if (!is.list(Wall)) {
        stop("Wall must be a list of matrices")
    }

    # n: number of matrices
    n <- length(Wall)
    if (n == 0) {
        stop("Wall is empty")
    }
    if (n == 1) {
        # message("Only one matrix provided. Returning the input matrix.")
        return(Wall[[1]])
    }

    check_wall_names <- function(Wall){
        # Checks if dimnames are consistant across all matrices in Wall
        name_match <- function(names_A, names_B){
            return(identical(dimnames(names_A), dimnames(names_B)))
        }
        return(all(unlist(lapply(Wall, FUN=name_match, Wall[[1]]))))
    }

    # Check if Wall names are consistent across all matrices in Wall
    wall.name.check <- check_wall_names(Wall)
    wall.names <- dimnames(Wall[[1]])
    if(!wall.name.check){
        warning("Dim names not consistent across all matrices in Wall.
            Returned matrix will have no dim names.")
    }

    # Check weights length
    if (length(weights) != n) {
        stop("Length of weights must match the number of matrices in Wall.")
    }
    # Ensure no negative weights
    if (any(weights < 0)) {
        stop("Weights must be non-negative.")
    }
    # Normalize weights
    weights <- weights / sum(weights)

    # Filter out matrices with zero weight
    non_zero_w <- which(weights > 0)
    Wall <- Wall[non_zero_w]
    weights <- weights[non_zero_w]

    # If there is only one remaining matrix, return it
    n <- length(Wall)
    if (n == 1) {
        return(Wall[[1]])
    }

    # Check K
    if (K > nrow(Wall[[1]])-1) {
        # k-nearest neighbors does not include itself
        stop("K is too big compared to the number of samples in the matrices.")
    }

    # Normalization method for affinity matrices
    normalize <- function(X){
        row.sum.mdiag <- rowSums(X) - diag(X)
        #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2*(row.sum.mdiag))
        diag(X) <- 0.5
        return(X)
    }

    # Normalize different networks to avoid scale problems.
    for(i in 1:n){
        Wall[[i]] <- normalize(Wall[[i]])
        Wall[[i]] <- (Wall[[i]] + t(Wall[[i]])) / 2 # Symmetrize matrix
    }

    # Internal SNFtool package function: .dominateset()
    .dominateset <- function(xx, KK = 20) {
        # This function outputs the top KK neighbors
        zero <- function(x) {
            s = sort(x, index.return = TRUE)
            x[s$ix[1:(length(x) - KK)]] = 0
            return(x)
        }
        normalize <- function(X)  X / rowSums(X)
        A = matrix(0, nrow(xx), ncol(xx))

        for (i in 1:nrow(xx)) {
            A[i, ] = zero(xx[i, ])
        }
        return(normalize(A))
    }

    # Calculate the local transition matrix using .dominateset() from SNFtool package
    newW <- vector("list", n)
    for(i in 1:n){
        newW[[i]] <- .dominateset(Wall[[i]], K)
    }

    # Perform the diffusion for t iterations
    for (i in 1:t) {

        # List to store the next iteration's matrices
        nextW <- vector("list", n)

        for(j in 1:n){
            sumWJ <- matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
            for(k in 1:n){
                if(k != j) {
                    sumWJ <- sumWJ + (Wall[[k]]*weights[k]) # Weighted sum of the other matrices
                }
            }
            nextW[[j]] <- newW[[j]] %*% (sumWJ/(n-1)) %*% t(newW[[j]]) # Diffusion step
        }

        # Normalize each new obtained networks
        for(j in 1:n){
            Wall[[j]] <- normalize(Wall[[j]])
            Wall[[j]] <- Wall[[j]] * weights[j] # Apply weights to each matrix
            Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2 # Symmetrize after diffusion
        }
    }

    # Construct the combined affinity matrix by summing diffused matrices
    W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:n){
        W <- W + Wall[[i]]
    }

    # Normalize and symmetrize the combined matrix
    W <- W / n
    W <- normalize(W)
    W <- (W + t(W)) / 2

    if(wall.name.check){
        dimnames(W) <- wall.names
    }

    return(W)
}


#' Impute cluster assignments for missing cells by similarity
#'
#' @description This function imputes cluster assignments for cells missing in
#' some omics by leveraging nearest neighbor cells in other omic matrices.
#'
#' @param mat_list A named list of log ratio *cells x features* matrices where
#'   each matrix corresponds to a single omic dataset (`list`). Rows are cells,
#'   and columns are features.
#' @param clusters A named vector of cluster assignments for cells (`numeric` or
#'   `character` vector). The vector names must correspond to the names of the
#'   common cells across omics (matching row names in `mat_list`). The clusters
#'   names can be as integer, numeric or character values.
#' @param knn_imp Integer specifying the number of nearest neighbors cells to
#'   use for imputation (`integer`). Must be a positive integer. Default is
#'   `10`.
#'
#' @return A named vector of combining the original clusters assignments for
#'   common cells across omics (given by the `clusters` argument) and the
#'   imputed cluster assignments for cells missing in at least one omic matrix.
#'
#' @importFrom RANN nn2
#'
#' @details
#' The function operates in the following steps:
#' 1. Identifies cells missing in specific matrices.
#' 2. Finds the k-nearest neighbors for missing cells in matrices where they are
#' present.
#' 3. Imputes cluster assignments for missing cells based on the clusters
#' assigned to their neighbors.
#' 4. Resolves ties (two major clusters found among neighbors) by selecting the
#' one of the first nearest neighbor.
#'
#' The imputation is performed separately for each omic dataset, and the results
#' are aggregated to provide final cluster assignments.
#'
#' @examples
#' # Create matrices with some cells missing in one or the other
#' set.seed(42)
#' mat1 <- matrix(runif(100), nrow = 20)
#' mat2 <- matrix(runif(100), nrow = 20)
#' rownames(mat1) <- paste0("Cell", 1:20)
#' rownames(mat2) <- paste0("Cell", c(1:5, 11:25))
#' mat_list <- list(ATAC = mat1, RNA = mat2)
#'
#' # Create cluster assignments for common cells
#' common_cells <- intersect(rownames(mat1), rownames(mat2))
#' clusters <- setNames(sample(1:4, length(common_cells), replace = TRUE), common_cells)
#'
#' # Check the inputs
#' print(common_cells)
#' print(rownames(mat_list$ATAC))
#' print(rownames(mat_list$RNA))
#'
#' # Impute cluster assignments for missing cells
#' imputed_clusters <- imputeClusters(mat_list, clusters, knn_imp = 3)
#'
#' # View the imputed cluster assignments
#' print(imputed_clusters)
#'
#' @export
#'
imputeClusters <- function(mat_list,
                           clusters,
                           knn_imp = 10) {

    # Argument checks
    if (!is.list(mat_list) || !all(sapply(mat_list, is.matrix))) {
        stop("`mat_list` must be a named list of matrices.")
    }

    stopifnot("`mat_list` must have names corresponding to omic datasets." = !is.null(names(mat_list)))

    stopifnot("`clusters` must be a named vector." = is.vector(clusters) & !is.null(names(clusters)))

    if (!is.numeric(knn_imp) || length(knn_imp) != 1 || knn_imp <= 0 || knn_imp %% 1 != 0) {
        stop("`knn_imp` must be a positive integer.")
    }

    # Get cell barcodes
    all_cells <- sort(Reduce(union, lapply(mat_list, rownames))) # all cells across matrices
    common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames))) # cells common to all matrices
    cells_NA_all <- setdiff(all_cells, common_cells) # cells missing from at least one matrix

    # If no common cells in any matrix
    stopifnot("No common cells found across matrices." = length(common_cells) > 0)

    # If no cells are missing in any matrix, the function exits early
    if (length(setdiff(all_cells, common_cells)) == 0) {
        warning("No missing cells detected; returning original cluster assignments.")
        return(clusters)
    }

    # Check cell names in clusters
    stopifnot("The names of `clusters` must correspond to cell names provided as row names in `mat_list`."
              = all(names(clusters) %in% all_cells))

    # Filter common_cells in clusters if necessary
    diff_cells <- length(setdiff(names(clusters), common_cells))

    if (diff_cells > 0) {
        warning(diff_cells, " cells in `clusters` are not common cells across omics and are removed.")
        clusters <- clusters[names(clusters) %in% common_cells]
    }

    # Check that knn cells are not higher than common cells
    stopifnot("`knn_imp` cannot be greater than the number of common cells." = knn_imp < length(common_cells))

    # ---

    # Initialize the output list for imputed clusters
    out <- list()
    cells_missing <- list()

    # Loop through each omic dataset in the list
    for (omic in names(mat_list)) {

        # Identify cells missing from the current omic dataset
        cells_NA <- setdiff(all_cells, rownames(mat_list[[omic]]))

        # Get all other omic matrices (excluding the current one)
        other_mat <- mat_list[names(mat_list) != omic]

        # Initialize a list to store k-nearest neighbors for each missing cell
        knn_list <- list()

        # Loop through other omic matrices to find neighbors
        for (j in seq_along(other_mat)) {

            # Identify common cells between the current matrix and the other matrix
            common_cells_mat <- intersect(rownames(mat_list[[omic]]), rownames(other_mat[[j]]))

            # Identify missing cells in the current matrix that are present in the other matrix
            cells_NA_other <- cells_NA[which(cells_NA %in% rownames(other_mat[[j]]))]

            # Skip if there are no missing cells to impute
            if (length(cells_NA_other) > 0) {

                # Find k-nearest neighbors using the other matrix
                nearest <- RANN::nn2(
                    data = other_mat[[j]][common_cells_mat, ],
                    query = other_mat[[j]][cells_NA_other, , drop = FALSE],
                    k = knn_imp
                )

                # Extract neighbor indices and map them to cell barcodes
                knn <- nearest[["nn.idx"]]
                rownames(knn) <- cells_NA_other
                knn <- apply(knn, c(1, 2), function(x) common_cells_mat[x])

                # Add NA placeholders for cells not found in this matrix
                missing_knn <- matrix(
                    NA,
                    nrow = length(setdiff(
                        cells_NA_all, cells_NA_other
                    )),
                    ncol = ncol(knn),
                    dimnames = list(
                        setdiff(cells_NA_all, cells_NA_other),
                        colnames(knn)
                    )
                )
                knn <- rbind(knn, missing_knn)
                knn <- knn[sort(rownames(knn)), , drop = FALSE]

                # Assign clusters for each k value using nearest neighbors
                clusters_NA <- matrix(
                    clusters[match(knn, names(clusters))],
                    nrow = nrow(knn),
                    ncol = ncol(knn),
                    dimnames = list(rownames(knn))
                )

                # Store the clusters for this matrix
                knn_list[[names(other_mat[j])]] <- clusters_NA
            }
        }

        # Combine cluster assignments from all other matrices
        clusters_NA <- do.call(cbind, knn_list)

        # Store results
        out[[omic]] <- clusters_NA
        cells_missing[[omic]] <- setdiff(cells_NA_all, cells_NA_other)
    }

    # Reduce the list of cells missing per omic to a vector for all omics
    cells_missing <- Reduce(c, cells_missing)

    # Combine imputed clusters across all omics datasets
    clusters_table <- do.call(cbind, out)
    clusters_table <- clusters_table[cells_missing, , drop = FALSE] # order per omic instead of alphabetically

    # Assign final cluster to each cell based on the majority vote among neighbors
    clusters_imp <- sapply(rownames(clusters_table), function(p) {
        clus_NA <- sort(summary(as.factor(na.omit(
            clusters_table[p, ]
        ))), decreasing = TRUE)

        # Handle ties: pick the cluster of the first matching nearest neighbor
        if (any(duplicated(clus_NA[which(clus_NA == max(clus_NA))]))) {
            # if more than 1 major cluster exists (identical high number of knn cells assigned to two clusters)
            final_cluster <- clusters_table[p, clusters_table[p, ] %in% names(which(clus_NA == max(clus_NA)))][1]
        } else {
            # else take the cluster corresponding to the highest number of knn cells
            final_cluster <- names(clus_NA)[1] # returns character
        }

        # Make sure the format is kept
        if(is.character(clusters)) final_cluster <- as.character(final_cluster)
        if(is.integer(clusters)) final_cluster <- as.integer(final_cluster)
        if(is.numeric(clusters)) final_cluster <- as.numeric(final_cluster)

        return(final_cluster)
    })

    # Combine original clusters with imputed clusters
    clusters_complete <- c(clusters, clusters_imp)

    return(clusters_complete)
}



