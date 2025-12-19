# Multi Omics Integration and Clustering on a `muscadet` Object

Performs integration of multi omics and clustering of cells based on log
ratio data contained in a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. Two methods are available for integration and clustering of
common cells between omics:

- Method `seurat` uses nearest neighbors for integration followed by
  graph-based clustering.

- Method `hclust` uses Similarity Network Fusion (SNF) for integration
  followed by hierarchical clustering. Then, clusters are imputed for
  cells missing data in at least one omic, by similarity using nearest
  neighbor cells. Finally, silhouette widths are computed on the
  integrated distance matrix to help identify the optimal clustering
  partition.

## Usage

``` r
clusterMuscadet(
  x,
  method = c("seurat", "hclust"),
  omics = NULL,
  knn_imp = 10,
  quiet = FALSE,
  ...
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing omics data and previously computed log R ratio
  matrices (`muscadet`).

- method:

  The clustering method to apply (`character` string). One of
  [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  should be provided, and for method `"hclust"`, arguments for
  [`cluster_hclust()`](https://icagen.github.io/muscadet/reference/cluster_hclust.md)
  should be provided. Note: The `"seurat"` method can only be applied
  for a maximum of 2 omics. Default is `"seurat"`.

- omics:

  Optional character vector specifying omic names to use for clustering.
  Must match names of available omics in the `x` muscadet object. If
  `NULL` (default), all available omics are used.

- knn_imp:

  Number of k nearest neighbors to use for imputing cluster assignments
  of cells missing in one or more omics. Only relevant for more than one
  omic.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

- ...:

  Arguments passed on to
  [`cluster_seurat`](https://icagen.github.io/muscadet/reference/cluster_seurat.md),
  [`cluster_hclust`](https://icagen.github.io/muscadet/reference/cluster_hclust.md)

  `res_range`

  :   A numeric non-negative vector specifying the resolution values to
      use for
      [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
      (`numeric` vector). Default is `c(0.1, 0.2, 0.3, 0.4, 0.5)`.

  `dims_list`

  :   A list of vectors of PC dimensions to use for each omic (`list`).
      Must match the length of `mat_list` (e.g., list(1:8) for 1 omic ;
      list(1:8, 1:8) for 2 omics). Default is the first 8 dimensions for
      each provided omic.

  `algorithm`

  :   Integer specifying the algorithm for modularity optimization by
      [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
      (`1` = original Louvain algorithm; `2` = Louvain algorithm with
      multilevel refinement; `3` = SLM algorithm; `4` = Leiden
      algorithm). Default is `4` for Leiden algorithm as recommended,
      see
      [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
      Details section.

  `knn_seurat`

  :   Integer specifying the number of nearest neighbors used for graph
      construction with
      [`Seurat::Seurat()`](https://satijalab.org/seurat/reference/Seurat-package.html)
      functions
      [`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)
      (`k.param`) or
      [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
      (`k.nn`) (`integer`). Default is `20`.

  `knn_range_seurat`

  :   Integer specifying the approximate number of nearest neighbors to
      compute for
      [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
      (`knn.range`) (`integer`). Default is `200`.

  `k_range`

  :   A numeric vector of integers (≥2) specifying the cluster
      numbers (k) to extract from hierarchical clustering (`numeric`
      vector). Default is from 2 to 10.

  `dist_method`

  :   A string specifying the distance method for
      [`Rfast::Dist()`](https://rdrr.io/pkg/Rfast/man/Dist.html) (e.g.,
      `"euclidean"`, `"manhattan"`, `"cosine"`) (`character` string).
      Default is `"euclidean"`.

  `hclust_method`

  :   A string specifying the hierarchical clustering linkage method for
      [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
      (e.g., `"ward.D"`, `"average"`) (`character` string). Default is
      `"ward.D"`.

  `weights`

  :   A numeric vector of non-negative values of length equal to the
      number of omic (internally normalized to sum to 1) (`numeric`
      vector). It specifies the relatives weights of each omic for SNF
      with
      [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md).
      Omics with a weight of 0 will not contribute to the clustering. If
      `NULL` (default), weights are uniform.

## Value

The input
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object with its `clustering` slot updated. This slot contains:

- params:

  List of parameters used for clustering (`list`).

- ...:

  Output objects depending on the method. See
  [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  or
  [`cluster_hclust()`](https://icagen.github.io/muscadet/reference/cluster_hclust.md).

- clusters:

  A named list of cluster partitions (named vectors of cluster labels)
  for all cells (imputed clusters assignments for non-common cells), for
  each value in `k_range` or `res_range` (`list`).

- silhouette:

  A list of silhouette objects and widths for each cluster partition
  (`list`).

- partition.opt:

  Name of the optimal cluster partition based on maximum average
  silhouette width.

## See also

Methodology and functionality:

- [muscadet](https://icagen.github.io/muscadet/reference/muscadet-class.md)

- [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  for graph-based clustering using
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html).

- [`cluster_hclust()`](https://icagen.github.io/muscadet/reference/cluster_hclust.md)
  for hierarchical clustering of SNF-fused distances.

- [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md)
  for weighted Similarity Network Fusion (SNF).

- [`imputeClusters()`](https://icagen.github.io/muscadet/reference/imputeClusters.md)
  for imputing cluster labels across omics.

Visualization:

- [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  to plot clustering result as heatmap.

- [`plotSil()`](https://icagen.github.io/muscadet/reference/plotSil.md)
  to plot silhouette widths.

- [`plotIndexes()`](https://icagen.github.io/muscadet/reference/plotIndexes.md)
  to plot several normalized cluster validation indexes.

Select clusters to continue with CNA calling:

- [`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md)
  to assign final cluster assignments in the `muscadet` object after
  cluster partition validation.

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")

# Perform clustering with "seurat" method
exdata_muscadet <- clusterMuscadet(
  x = exdata_muscadet,
  method = "seurat",
  res_range = c(0.1, 0.3),
  dims_list = list(1:10, 1:10),
  knn_seurat = 10, # adapted to low number of cells in example data
  knn_range_seurat = 30 # adapted to low number of cells in example data
)
#> Clustering method: 'seurat'
#> Resolutions to compute: 0.1, 0.3
#> Number of selected dimensions: 10, 10
#> Clustering algorithm selected: 4 (Leiden)
#> Performing PCA...
#> Finding neighbors and constructing graph...
#> Computing UMAP...
#> Finding clusters...
#> Imputing clusters...
#> Computing Silhouette scores...
#> Done.

# Perform clustering with "hclust" method
exdata_muscadet2 <- clusterMuscadet(
  x = exdata_muscadet,
  k_range = 2:4,
  method = "hclust",
  dist_method = "euclidean",
  hclust_method = "ward.D",
  weights = c(1, 1)
)
#> Clustering method: 'hclust'
#> Partitions k to compute: 2, 3, 4
#> Distance method selected: euclidean
#> Hierarchical clustering method selected: ward.D
#> Computing distance matrices...
#> Computing affinity matrices...
#> Performing SNF integration...
#> Computing UMAP...
#> Performing hierarchical clustering...
#> Imputing clusters...
#> Computing Silhouette scores...
#> Done.

# Retrieve cluster assignments
clusters <- exdata_muscadet$clustering$clusters
lapply(clusters, table)
#> $`0.1`
#> 
#>  1 
#> 77 
#> 
#> $`0.3`
#> 
#>  1  2 
#> 41 36 
#> 

if (FALSE) { # \dontrun{

# Plot clustree
library(clustree)
partitions <- lapply(exdata_muscadet$clustering$clusters, as.data.frame)
partitions <- do.call(cbind, partitions)
colnames(partitions) <- paste0("res_", names(exdata_muscadet$clustering$clusters))
clustree(partitions, prefix = "res_")
} # }
```
