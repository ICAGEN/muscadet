# Plot clustering validation indexes for a `muscadet` object

It generates a plot of clustering validation indexes for a clustering
partition within a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. The index values are computed only using distances between
common cells across omics in the
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Usage

``` r
plotIndexes(x, index = NULL, colors = NULL, title = NULL)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing clustering data (generated using
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).

- index:

  Character vector specifying one or more validation indexes to plot
  among `"silhouette"`, `"dunn2"`, `"daviesbouldin"`, `"pearsongamma"`,
  and `"c"`. If `NULL`, by default all available indexes are included.
  If multiple indexes are selected, the values are normalized for
  comparability.

- colors:

  Vector of colors for each index in the plot (`character` vector).
  Default is `NULL`, which uses predefined colors for the indexes.

- title:

  Character string for the title of the plot (`character` string). If
  `NULL`, a default title is generated.

## Value

A ggplot object visualizing the clustering validation indexes across
different clustering partitions (`res` resolution or `k` number of
clusters depending on the used clustering method).

## Details

The function computes several clustering validation indexes, including:

- **Silhouette**: Measures how similar an object is to its own cluster
  compared to others (see
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)).
  Average of individual silhouette widths.

- **Dunn2**: The ratio of the smallest distance between observations in
  different clusters to the largest within-cluster distance (see
  [`fpc::cluster.stats()`](https://rdrr.io/pkg/fpc/man/cluster.stats.html)
  `$dunn2`). Minimum average dissimilarity between two cluster / maximum
  average within cluster dissimilarity.

- **Davies-Bouldin**: Measures cluster compactness and separation (see
  [`clusterSim::index.DB()`](https://rdrr.io/pkg/clusterSim/man/index.DB.html)).

- **Pearson's Gamma**: Evaluates the goodness of clustering based on
  correlation (see
  [`fpc::cluster.stats()`](https://rdrr.io/pkg/fpc/man/cluster.stats.html)
  `$pearsongamma`). Correlation between distances and a 0-1-vector where
  0 means same cluster, 1 means different clusters. "Normalized gamma"
  in Halkidi et al. (2001).

- **C Index** (Hubert & Levin C index): Measures the internal cluster
  quality compared to random data (see
  [`clusterSim::index.C()`](https://rdrr.io/pkg/clusterSim/man/index.C.html)).

If multiple indexes are selected, the values are normalized to fall
between 0 and 1. For indexes that are better when minimized
("pearsongamma" and "c"), their values are reversed for easier
comparison. The partition for which the mean of indexes is maximal is
highlighted with a dot.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Load example muscadet object
# data("muscadet_obj")

# Plot all indexes
plotIndexes(muscadet_obj)
ggsave("plot_indexes.png", width = 8, height = 4)

# Plot a specific index
plotIndexes(muscadet_obj, index = "silhouette")
ggsave("plot_indexes_sil.png", width = 7, height = 4)
} # }
```
