
test_that("clusterMuscadet() returns an updated muscadet object", {
    data(muscadet_obj)
    # remove clustering result
    muscadet_obj@clustering <- list()

    obj <- clusterMuscadet(
        muscadet_obj,
        method = "seurat",
        res_range = c(0.6, 0.8),
        knn_seurat = 10,
        knn_range_seurat = 30
    )
    obj2 <- clusterMuscadet(
        muscadet_obj,
        method = "hclust",
        dist_method = "euclidean",
        hclust_method = "ward.D",
        k_range = 2:5
    )
    expect_s4_class(obj, "muscadet")
    expect_s4_class(obj2, "muscadet")
    expect_true(length(obj$clustering) != 0)
    expect_true(length(obj2$clustering) != 0)
})

test_that("imputeClusters() returns a named vectors", {
    # Create matrices with some cells missing in one or the other
    mat1 <- matrix(runif(100), nrow = 20)
    mat2 <- matrix(runif(100), nrow = 20)
    rownames(mat1) <- paste0("Cell", 1:20)
    rownames(mat2) <- paste0("Cell", c(1:5, 11:25))
    mat_list <- list(ATAC = mat1, RNA = mat2)

    # Create cluster assignments for common cells
    clusters <- setNames(sample(1:3, 15, replace = TRUE), intersect(rownames(mat1), rownames(mat2)))

    # Impute cluster assignments for missing cells
    imputed_clusters <- imputeClusters(mat_list, clusters, knn_imp = 5)

    expect_type(imputed_clusters, "double")
    expect_length(imputed_clusters, 25)

    expect_identical(sort(names(imputed_clusters)), sort(union(rownames(mat1), rownames(mat2))))
})





