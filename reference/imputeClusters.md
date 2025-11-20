# Impute cluster assignments for missing cells by similarity

This function imputes cluster assignments for cells missing in some
omics by leveraging nearest neighbor cells in other omic matrices.

## Usage

``` r
imputeClusters(mat_list, clusters, knn_imp = 10)
```

## Arguments

- mat_list:

  A named list of log ratio *cells x features* matrices where each
  matrix corresponds to a single omic dataset (`list`). Rows are cells,
  and columns are features.

- clusters:

  A named vector of cluster assignments for cells (`numeric` or
  `character` vector). The vector names must correspond to the names of
  the common cells across omics (matching row names in `mat_list`). The
  clusters names can be as integer, numeric or character values.

- knn_imp:

  Integer specifying the number of nearest neighbors cells to use for
  imputation (`integer`). Must be a positive integer. Default is `10`.

## Value

A named vector of combining the original clusters assignments for common
cells across omics (given by the `clusters` argument) and the imputed
cluster assignments for cells missing in at least one omic matrix.

## Details

The function operates in the following steps:

1.  Identifies cells missing in specific matrices.

2.  Finds the k-nearest neighbors for missing cells in matrices where
    they are present.

3.  Imputes cluster assignments for missing cells based on the clusters
    assigned to their neighbors.

4.  Resolves ties (two major clusters found among neighbors) by
    selecting the one of the first nearest neighbor.

The imputation is performed separately for each omic dataset, and the
results are aggregated to provide final cluster assignments.

## Examples

``` r
# Create matrices with some cells missing in one or the other
set.seed(42)
mat1 <- matrix(runif(100), nrow = 20)
mat2 <- matrix(runif(100), nrow = 20)
rownames(mat1) <- paste0("Cell", 1:20)
rownames(mat2) <- paste0("Cell", c(1:5, 11:25))
mat_list <- list(ATAC = mat1, RNA = mat2)

# Create cluster assignments for common cells
common_cells <- intersect(rownames(mat1), rownames(mat2))
clusters <- setNames(sample(1:4, length(common_cells), replace = TRUE), common_cells)

# Check the inputs
print(common_cells)
#>  [1] "Cell1"  "Cell2"  "Cell3"  "Cell4"  "Cell5"  "Cell11" "Cell12" "Cell13"
#>  [9] "Cell14" "Cell15" "Cell16" "Cell17" "Cell18" "Cell19" "Cell20"
print(rownames(mat_list$ATAC))
#>  [1] "Cell1"  "Cell2"  "Cell3"  "Cell4"  "Cell5"  "Cell6"  "Cell7"  "Cell8" 
#>  [9] "Cell9"  "Cell10" "Cell11" "Cell12" "Cell13" "Cell14" "Cell15" "Cell16"
#> [17] "Cell17" "Cell18" "Cell19" "Cell20"
print(rownames(mat_list$RNA))
#>  [1] "Cell1"  "Cell2"  "Cell3"  "Cell4"  "Cell5"  "Cell11" "Cell12" "Cell13"
#>  [9] "Cell14" "Cell15" "Cell16" "Cell17" "Cell18" "Cell19" "Cell20" "Cell21"
#> [17] "Cell22" "Cell23" "Cell24" "Cell25"

# Impute cluster assignments for missing cells
imputed_clusters <- imputeClusters(mat_list, clusters, knn_imp = 3)

# View the imputed cluster assignments
print(imputed_clusters)
#>  Cell3 Cell11  Cell8 Cell22  Cell2 Cell12 Cell13 Cell14 Cell16 Cell20 Cell10 
#>      1      1      1      1      2      2      2      2      2      2      2 
#>  Cell6  Cell7  Cell9 Cell23 Cell24  Cell5 Cell17 Cell19 Cell25  Cell1  Cell4 
#>      2      2      2      2      2      3      3      3      3      4      4 
#> Cell15 Cell18 Cell21 
#>      4      4      4 
```
