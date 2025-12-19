# Multi Omics Clustering using Seurat Multi Modal Graph-based Clustering

Performs graph-based clustering of cells using
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html),
based on one or two log R ratio matrices (`mat_list`), including shared
nearest neighbors (SNN) graph construction on selected dimensions from
PCA (`dims_list`), to identify clusters of cells for each specified
resolution (`res_range`).

- For two omics: multimodal integration is performed using
  [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
  (weighted shared nearest neighbors graph). Only common cells between
  omics are used.

- For a single omic:
  [`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  (shared nearest neighbors graph) is used.

## Usage

``` r
cluster_seurat(
  mat_list,
  res_range = seq(0.1, 0.5, 0.1),
  dims_list = rep(list(1:8), length(mat_list)),
  algorithm = 4,
  leiden_method = "igraph",
  knn_seurat = 20,
  knn_range_seurat = 200,
  max_dim = 200,
  random.seed = 1,
  quiet = FALSE
)
```

## Arguments

- mat_list:

  A named list of log R ratio matrices (cells x features), one per omic
  layer (`list`).

- res_range:

  A numeric non-negative vector specifying the resolution values to use
  for
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  (`numeric` vector). Default is `c(0.1, 0.2, 0.3, 0.4, 0.5)`.

- dims_list:

  A list of vectors of PC dimensions to use for each omic (`list`). Must
  match the length of `mat_list` (e.g., list(1:8) for 1 omic ; list(1:8,
  1:8) for 2 omics). Default is the first 8 dimensions for each provided
  omic.

- algorithm:

  Integer specifying the algorithm for modularity optimization by
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  (`1` = original Louvain algorithm; `2` = Louvain algorithm with
  multilevel refinement; `3` = SLM algorithm; `4` = Leiden algorithm).
  Default is `4` for Leiden algorithm as recommended, see
  `cluster_seurat()` Details section.

- leiden_method:

  Character string to choose from which package for running leiden
  algorithm, either leidenbase ("leidenbase") or igraph ("igraph")
  packages (`character`). Default is "igraph".

- knn_seurat:

  Integer specifying the number of nearest neighbors used for graph
  construction with
  [`Seurat::Seurat()`](https://satijalab.org/seurat/reference/Seurat-package.html)
  functions
  [`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  (`k.param`) or
  [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
  (`k.nn`) (`integer`). Default is `20`.

- knn_range_seurat:

  Integer specifying the approximate number of nearest neighbors to
  compute for
  [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
  (`knn.range`) (`integer`). Default is `200`.

- max_dim:

  Integer specifying the maximum number of principal components to be
  used for PCA computation with
  [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) (`integer`).
  Default is `200`.

- random.seed:

  Integer specifying the seed of the random number generator, must be
  greater than 0 for Leiden algorithm. (`integer`). Default is `1`.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

A list containing:

- params:

  List of parameters used for clustering (`list`).

- pcs:

  List of principal components summaries for each omic (`list` of
  [`summary.prcomp`](https://rdrr.io/r/stats/prcomp.html)).

- nn:

  Nearest neighbors object
  ([`Neighbor`](https://satijalab.github.io/seurat-object/reference/Neighbor-class.html)).

- graph:

  Shared nearest neighbors graph
  ([`Graph`](https://satijalab.github.io/seurat-object/reference/Graph-class.html)).

- dist:

  Distance matrix derived from the graph (`matrix`).

- umap:

  UMAP coordinates (`matrix`).

- clusters:

  A named list of clustering results (vectors of cluster labels) for
  each value in `res_range` (`list`).

## Details

The Leiden algorithm (`algorithm = 4`) is recommended based on published
work and best-practice guidelines:

- Traag, V.A., Waltman, L. & van Eck, N.J. *From Louvain to Leiden:
  guaranteeing well-connected communities.* Sci Rep 9, 5233 (2019).
  <https://doi.org/10.1038/s41598-019-41695-z%3E>

- Heumos, L., Schaar, A.C., Lance, C. et al. *Best practices for
  single-cell analysis across modalities.* Nat Rev Genet (2023).
  <https://doi.org/10.1038/s41576-023-00586-w>
  <https://www.sc-best-practices.org/cellular_structure/clustering.html>

Explanation about igraph being the default package for Leiden
implementation (`leiden_method = "igraph"`) here:
https://github.com/satijalab/seurat/issues/9800

## See also

[Weighted Nearest Neighbor Analysis Vignette from
Seurat](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example muscadet object
# data("exdata_muscadet")

# Format input
# transpose matrices to: cells x features matrices
mat_list <- matLogRatio(exdata_muscadet)

# Run integration & clustering
result <- cluster_seurat(
    mat_list,
    res_range = c(0.1, 0.3, 0.5),
    knn_seurat = 10,
    knn_range_seurat = 30
)

# View results
lapply(result$clusters, table)
} # }
```
