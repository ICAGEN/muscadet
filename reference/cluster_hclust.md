# Multi Omics Clustering with SNF integration and Hierarchical Clustering

Performs the integration of log R ratio matrices (`mat_list`) using
Similarity Network Fusion (SNF) followed by hierarchical clustering, on
the integrated SNF matrix, to identify clusters of cells for each
specified k number of cluster (`k_range`). For more than 1 omic, the
integration and clustering are performed only on common cells between
omics.

## Usage

``` r
cluster_hclust(
  mat_list,
  k_range = seq(2, 10, 1),
  dist_method = "euclidean",
  hclust_method = "ward.D",
  weights = rep(1, length(mat_list)),
  knn_affinity = 40,
  var_affinity = 1,
  knn_SNF = 40,
  iter_SNF = 50,
  knn_umap = 20,
  quiet = FALSE
)
```

## Arguments

- mat_list:

  A named list of log R ratio matrices (cells x features), one per omic
  layer (`list`).

- k_range:

  A numeric vector of integers (≥2) specifying the cluster numbers (k)
  to extract from hierarchical clustering (`numeric` vector). Default is
  from 2 to 10.

- dist_method:

  A string specifying the distance method for
  [`Rfast::Dist()`](https://rdrr.io/pkg/Rfast/man/Dist.html) (e.g.,
  `"euclidean"`, `"manhattan"`, `"cosine"`) (`character` string).
  Default is `"euclidean"`.

- hclust_method:

  A string specifying the hierarchical clustering linkage method for
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
  (e.g., `"ward.D"`, `"average"`) (`character` string). Default is
  `"ward.D"`.

- weights:

  A numeric vector of non-negative values of length equal to the number
  of omic (internally normalized to sum to 1) (`numeric` vector). It
  specifies the relatives weights of each omic for SNF with
  [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md).
  Omics with a weight of 0 will not contribute to the clustering. If
  `NULL` (default), weights are uniform.

- knn_affinity:

  Integer specifying the number of nearest neighbors used when building
  affinity matrices with
  [`SNFtool::affinityMatrix()`](https://rdrr.io/pkg/SNFtool/man/affinityMatrix.html)
  (`integer`). Default is `40`.

- var_affinity:

  Numeric value for the variance parameter (Gaussian kernel width
  `sigma`) when building affinity matrix with
  [`SNFtool::affinityMatrix()`](https://rdrr.io/pkg/SNFtool/man/affinityMatrix.html)
  (`numeric`). Default is `1`.

- knn_SNF:

  Integer specifying the number of nearest neighbors used during the
  Similarity Network Fusion (SNF) with
  [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md)
  (`integer`). Default is `40`.

- iter_SNF:

  Integer specifying the number of iterations for SNF with
  [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md)
  (`integer`). Default is `50`.

- knn_umap:

  Integer specifying the number of nearest neighbors used for manifold
  approximation (UMAP) (see
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)).
  Default is `20`.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

A list containing:

- params:

  List of parameters used for clustering (`list`).

- SNF:

  Fused similarity matrix computed with SNF (`matrix`).

- dist:

  Distance matrix derived from the SNF similarity (`matrix`).

- hclust:

  Hierarchical clustering object from
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
  (`hclust`).

- umap:

  UMAP coordinates (`matrix`).

- clusters:

  A named list of clustering results (vectors of cluster labels) for
  each value in `k_range` (`list`).

## Details

The function calculates pairwise distances between cells within each
omic dataset using the specified `dist_method` (using only common cells
between omics).

It constructs affinity matrices based on these distances, applies SNF to
generate a fused similarity matrix.

`weights` can be assigned to each omic dataset to prioritize certain
data types over others, allowing users to tailor the analysis based on
the characteristics and importance of each dataset.

It then performs a hierarchical clustering using the specified
`hclust_method` to assign clusters to each common cells.

Results are given as cluster assignments for each number of cluster
specified by `k_range`.

## See also

Similarity Network Fusion:
[`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example muscadet object
# data("muscadet_obj")

# Format input
# transpose matrices to: cells x features matrices
mat_list <- lapply(muscadet::matLogRatio(muscadet_obj), t)

# Run integration & clustering
result <- cluster_hclust(mat_list, k_range = 2:4)

# View results
lapply(result$clusters, table)
} # }
```
