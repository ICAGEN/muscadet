#' Perform clustering on a muscadet object
#'
#' This function performs clustering on the cells within a
#' \code{\link{muscadet}} object using the log R ratio matrices. It calculates
#' pairwise distances, constructs affinity matrices, applies Similarity Network
#' Fusion (SNF) to integrate several omic, and assign clusters to cell. Clusters
#' for missing cells are imputed by similarity using k nearest-neighbor method.
#'
#' @param x A `muscadet` object containing omics data and previously calculated
#'   log R ratio matrices (`muscadet`).
#' @param dist_method String specifying the distance method for
#'   [Rfast::Dist()] (`character` string). Default is `"euclidean"`.
#' @param hclust_method String specifying the hierarchical clustering method
#'   for [fastcluster::hclust()] (`character` string). Default is `"ward.D"`.
#' @param knn_affinity Integer specifying the number of nearest neighbors used
#'   in affinity matrix construction with [SNFtool::affinityMatrix()]
#'   (`integer`). Default is `40`.
#' @param var_affinity Numeric value for the variance parameter in affinity
#'   matrix construction with [SNFtool::affinityMatrix()] (`numeric`). Default
#'   is `1`.
#' @param knn_SNF Integer specifying the number of nearest neighbors used for
#'   Similarity Network Fusion (SNF) with [muscadet::weightedSNF()] (`integer`).
#'   Default is `40`.
#' @param iter_SNF Integer specifying the number of iterations for SNF with
#'   [muscadet::weightedSNF()] (`integer`). Default is `50`.
#' @param weights Numeric vector of non-negative values of length equal to the
#'   number of omic (internally normalized to sum to 1) (`numeric` vector). It
#'   specifies the relatives weights of each omic for SNF with
#'   [muscadet::weightedSNF()]. Omics with a weight of 0 will not contribute to
#'   the clustering. Default assigns equal weights to all omics.
#' @param knn_imp Integer specifying the number of nearest neighbors to use for
#'   imputation of missing clusters with [muscadet::imputeClusters()]
#'   (`integer`). Default is `10`.
#' @param k_range Numeric vector specifying the range of clusters (k) to
#'   evaluate (`numeric` vector). Default is from 2 to 10.
#' @param no_aff Logical, whether to not go through a distance-to-affinity
#'   matrix conversion step if only one omic has a non-zero weight. Default is
#'   `FALSE` (distance matrix converted to affinity matrix).
#'
#'
#' @return
#' A \code{\link{muscadet}} object with updated clustering results
#' stored in the `clustering` slot, including:
#' \describe{
#'   \item{`params`}{Clustering parameters used (`list`).}
#'   \item{`SNF`}{Fused similarity matrix, result of Similarity Network Fusion
#'   (SNF) (`matrix`). If the muscadet object contains only one omic, it
#'   corresponds to the affinity matrix from this unique omic.}
#'   \item{`dist`}{Distance matrix, derived from the `SNF` matrix:
#'   `max(SNF) - SNF` (excluding the diagonal values) (`matrix`).}
#'   \item{`hclust`}{Hierarchical clustering object of class \code{\link{hclust}}, result of
#'   the clustering function [fastcluster::hclust()] with the chosen hclust
#'   method used (`hclust`).}
#'   \item{`clusters`}{Cluster assignments as a list for each given `k` from the
#'   `k_range` argument (`list`). Each element of the list is a named vector of
#'   clusters per cell.
#'   Note: element names being integers (`k`) in the character format it should
#'   be called as followed: `muscadet_obj$clustering$clusters[["3"]]` or
#'   `muscadet_obj$clustering$clusters[[as.character(k)]]`}
#'   \item{`silhouette`}{List of silhouette objects (`sil.obj`), average
#'   silhouette widths (`sil.w.avg`), and clusters average silhouette widths
#'   (`sil.w.avg.clusters`), per each `k` partition (`list`).}
#'   \item{`k.opt`}{Optimal `k` number of clusters based on silhouette average
#'   widths (`integer`).}
#' }
#'
#'
#' @details
#' The function calculates pairwise distances for cells within each omic dataset
#' using the specified `dist_method`.
#'
#' It constructs affinity matrices based on these distances, applies SNF to
#' generate a fused similarity matrix.
#'
#' Weights can be assigned to each omic dataset to prioritize certain data types
#' over others, allowing users to tailor the analysis based on the
#' characteristics and importance of each dataset.
#'
#' It then performs a hierarchical clustering using the specified
#' `hclust_method` to assign clusters to each cell with data in all omics.
#'
#' Then, it imputes cluster assignments for missing cells based on their
#' nearest neighbors.
#'
#' Finally, a default optimal number of cluster (`k`) is found based on average
#' silhouette widths.
#'
#' @seealso
#' * \code{\link{muscadet-class}}
#' * Details on Similarity Network Fusion: [muscadet::weightedSNF()]
#' * Plot cluster result as heatmap with [muscadet::heatmapMuscadet()].
#' * Plot silhouette widths with [muscadet::plotSil()].
#' * Plot several cluster validation indexes with [muscadet::plotIndexes()].
#' * After cluster partition validation, assign final cluster assignments with
#' [muscadet::assignClusters()].
#'
#' @importFrom Rfast Dist
#' @importFrom SNFtool affinityMatrix
#' @importFrom stats as.dist
#' @importFrom fastcluster hclust
#' @importFrom dendextend cutree
#' @importFrom cluster silhouette
#' @importFrom methods slot slot<-
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @examples
#' # Load example muscadet object
#' data(muscadet_obj)
#'
#' # Clustering on muscadet object
#' muscadet_obj <- clusterMuscadet(
#'   muscadet_obj,
#'   dist_method = "euclidean",
#'   hclust_method = "ward.D",
#'   k_range = 2:5
#' )
#'
#' # Cluster assignments (for k = 3)
#' muscadet_obj$clustering$clusters[["3"]]
#' table(muscadet_obj$clustering$clusters[["3"]])
#'
#' # Silhouette width per k
#' unlist(muscadet_obj$clustering$silhouette[["sil.w.avg"]])
#'
#' # Optimal k according to silhouette width index
#' muscadet_obj$clustering$k.opt
#'
clusterMuscadet <- function(x, # muscadet object
                            dist_method = "euclidean",
                            hclust_method = "ward.D",
                            knn_affinity = 40,
                            var_affinity = 1,
                            knn_SNF = 40,
                            iter_SNF = 50,
                            weights = rep(1, length(slot(x, "omics"))),
                            knn_imp = 10,
                            k_range = seq(2, 10, 1),
                            no_aff = FALSE) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))

    # Add clustering parameters to the muscadet object
    slot(x, "clustering")[["params"]] <- list(
        dist_method = dist_method,
        hclust_method = hclust_method,
        knn_affinity = knn_affinity,
        var_affinity = var_affinity,
        knn_SNF = knn_SNF,
        iter_SNF = iter_SNF,
        weights = weights,
        knn_imp = knn_imp
    )

    # Prepare a named list of k
    k_list <- as.list(k_range)
    k_list <- setNames(k_list, k_range)

    # Compute distance and affinity matrices -----------------------------------

    # Check that log ratio matrices exists
    stopifnot(
        "At least one log ratio matrix is missing for one omic." =
            !any(unlist(lapply(
                matLogRatio(x), is.null
            )))
    )
    # Extract and transpose log ratio matrices for each omic
    mat_list <- lapply(muscadet::matLogRatio(x), t)

    # Get common cell barcodes across all omics
    common_cells <- sort(Reduce(intersect, lapply(mat_list, rownames)))

    # Compute pairwise distances for each omic
    dist_list <- lapply(
      mat_list,
      function(mat) {
        mat <- mat[common_cells, ]
        dist <- Rfast::Dist(mat, method = dist_method)
        # for version <=2.1, Rfast returns a similarity matrix for cosine with 0
        # values in diagonal. The result should be converted into distance.
        # See https://github.com/RfastOfficial/Rfast/issues/119 for more info
        if(dist_method == "cosine" & utils::packageVersion("Rfast") <= "2.1.0") {
            diag(dist) <- 1
            dist <- 1-dist
        }
        dimnames(dist) <- list(rownames(mat), rownames(mat))
        return(dist)
      }
    )

    omic_index <- which(weights != 0)

    if (no_aff == FALSE) {
        # Compute affinity matrices for each omic
        aff_list <- lapply(dist_list, function(dist) {
            aff <- SNFtool::affinityMatrix(dist, K = knn_affinity, sigma = var_affinity)
            aff <- aff[sort(rownames(aff)), sort(colnames(aff))]
            return(aff)
        })

        # Similarity Network Fusion (SNF) -------------------------------------------

        # Apply Similarity Network Fusion (SNF)
        matSNF <- weightedSNF(aff_list,
                              K = knn_SNF,
                              t = iter_SNF,
                              weights = weights)
        slot(x, "clustering")[["SNF"]] <- matSNF

        # Convert the fused affinity matrix to a distance matrix
        tmp <- matSNF
        diag(tmp) <- 0
        dist <- stats::as.dist(max(tmp) - matSNF)
        slot(x, "clustering")[["dist"]] <- as.matrix(dist)

    } else if (no_aff == TRUE & length(omic_index) == 1) {

        dist <- as.dist(dist_list[[omic_index]])

        slot(x, "clustering")[["SNF"]] <- NULL
        slot(x, "clustering")[["dist"]] <- as.matrix(dist)
    }

    # Clustering ----------------------------------------------------------------

    # Perform hierarchical clustering
    hc <- fastcluster::hclust(dist, hclust_method)
    slot(x, "clustering")[["hclust"]] <- hc

    # Cut the dendrogram to generate clusters for the specified k range
    clusters <- lapply(k_list, function(k) {
        cl <- dendextend::cutree(hc, k, order_clusters_as_data = FALSE)
        return(cl)
    })

    # Clusters imputation ------------------------------------------------------

    # Impute clusters assignments to missing cells based on nearest neighbors (if several omics)
    if (length(mat_list) > 1) {
        clusters_imp <- lapply(clusters, function(clus) {
            imputeClusters(mat_list, clus, knn_imp = knn_imp)
        })
        slot(x, "clustering")[["clusters"]] <- clusters_imp
    } else {
        slot(x, "clustering")[["clusters"]] <- clusters
    }

    # Find optimal k number of clusters  ---------------------------------------
    sil <- lapply(k_list, function(k) {
        cluster::silhouette(as.vector(clusters[[as.character(k)]]), dist)
    })

    sil_widths <- lapply(sil, function(x) summary(x)$avg.width)
    sil_widths_cl <- lapply(sil, function(x) summary(x)$clus.avg.widths)

    sil_widths_vec <- unlist(sil_widths)
    k_opt <- names(sil_widths_vec[sil_widths_vec == max(sil_widths_vec)])

    slot(x, "clustering")[["silhouette"]] <- list(sil.obj = sil,
                                                  sil.w.avg = sil_widths,
                                                  sil.w.avg.clusters = sil_widths_cl)

    slot(x, "clustering")[["k.opt"]] <- k_opt

    # Return the updated muscadet object
    return(x)
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



