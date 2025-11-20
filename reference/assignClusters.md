# Assign a cluster assignment to a `muscadet` object

Add the user selected cluster assignments to cells in a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. This function allows the user to choose the cluster assignments
they consider to fit the data and their requirements, or cluster
assignments based on other data and methods.

## Usage

``` r
assignClusters(
  x,
  partition = NULL,
  clusters = NULL,
  mapping = NULL,
  redo_imputation = TRUE,
  knn_imp = 10
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object (`muscadet`).

- partition:

  Value specifying the clustering partition to choose from the
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object (`numeric` or `character`). It should be either the resolution
  or the k number of cluster (k) used for clustering depending on the
  clustering method (`res_range` or `k_range` with
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).
  Should be provided if `clusters` is `NULL`.

- clusters:

  A custom named vector of cluster assignments (`vector`). The vector
  names must match cell names in the muscadet object `x`, at least
  cluster assignments for all common cells must be provided if
  `redo_imputation` is set to true, otherwise, all cells within the
  muscadet object `x` must be provided. Should be provided if
  `partition` is `NULL`.

- mapping:

  Optional named vector specifying how to remap cluster values
  (`vector`). The names of the vector correspond to the original cluster
  values, and the values are the remapped cluster values. For example,
  `c("1" = 1, "2" = 1, "3" = 2, "4" = 3)` would merge clusters 1 and 2
  into 1, cluster 3 into 2, and cluster 4 into 3.

- redo_imputation:

  Logical. If `TRUE` (default), reruns the imputation process to assign
  clusters to cells with missing data. This ensures that imputed
  clusters are updated if the clustering has changed due to remapping or
  to the use of custom clusters.

- knn_imp:

  Integer specifying the number of nearest neighbors cells to use for
  imputation (`integer`). Must be a positive integer. Default is `10`.

## Value

A
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object updated with the user chosen cluster assignments in
`muscadet_obj$cnacalling$clusters`.

## Details

- The clusters can be taken directly from the
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object clustering results with setting the `partition` argument (e.g.
  `muscadet_obj$clustering$clusters[["0.8"]]` for res=`0.8`).

- A custom vector of cluster assignments can be attributed using the
  `clusters` argument.

- Either way, the clusters assignments can be rearranged using the
  `mapping` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example muscadet object
# data("muscadet_obj")

# Select clustering result for partition = 0.6
muscadet_obj <- assignClusters(muscadet_obj, partition = 0.6)
table(muscadet_obj$cnacalling$clusters)

# Assign custom clusters
set.seed(42)
cell_names <- Reduce(union, SeuratObject::Cells(muscadet_obj))
n1 <- sample(1:length(cell_names), 1)
n2 <- length(cell_names) - n1
custom_clusters <- setNames(c(rep.int(1, n1), rep.int(2, n2)), cell_names)
table(custom_clusters)
muscadet_obj <- assignClusters(muscadet_obj, clusters = custom_clusters)
table(muscadet_obj$cnacalling$clusters)

# Assign clusters with remapping
# example to remap from partition=0.8 with merging of clusters 2 and 3
clusters <- muscadet_obj$clustering$clusters[["0.8"]]
table(clusters) # 3 clusters
mapping <- c("1" = 1, "2" = 2, "3" = 2) # remap to 2 clusters

muscadet_obj <- assignClusters(muscadet_obj, clusters = clusters, mapping = mapping)
table(muscadet_obj$cnacalling$clusters)
# check original and remapped clusters
table(clusters, muscadet_obj$cnacalling$clusters)

muscadet_obj <- assignClusters(muscadet_obj, partition = 0.8, mapping = mapping)
table(muscadet_obj$cnacalling$clusters)
# check original and remapped clusters
table(muscadet_obj$clustering$clusters[["0.8"]],
      muscadet_obj$cnacalling$clusters)

# Visualize clusters on heatmap
heatmapMuscadet(
    muscadet_obj,
    partition = 0.8,
    filename = file.path("heatmap_muscadet_res0.8.png"),
    title = "Example sample | res=0.8"
)
heatmapMuscadet(
    muscadet_obj,
    clusters = muscadet_obj$cnacalling$clusters,
    filename = file.path("heatmap_muscadet_custom_res0.8.png"),
    title = "Example sample | rearranged clusters from res=0.8"
)
} # }
```
