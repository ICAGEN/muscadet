# Silhouette plot for `muscadet` object

Generate a silhouette plot for a specified clustering partition within a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Usage

``` r
plotSil(x, partition, colors = NULL, title = NULL, annotations = TRUE)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing clustering data (using
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).

- partition:

  Value specifying the clustering partition to plot (`numeric` or
  `character`). It should be either the resolution or the k number of
  cluster (k) used for clustering depending on the clustering method
  (`res_range` or `k_range` with
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).

- colors:

  Vector of colors for the cluster annotation (`character` vector).
  Default is `NULL`, which uses predefined colors.

- title:

  Character string for the title of the plot (`character` string). If
  `NULL`, a default title is generated.

- annotations:

  `TRUE` or `FALSE` (`logical`). Whether to add annotations per
  clusters. By default: `TRUE`.

## Value

A ggplot object representing the silhouette plot.

## Examples

``` r
if (FALSE) { # \dontrun{
library("ggplot2")

# Load example muscadet object
# data("exdata_muscadet")
plotSil(exdata_muscadet, partition = 0.3)

# Loop over partitions
for (p in names(exdata_muscadet$clustering$clusters)) {
    plot <- plotSil(exdata_muscadet, p)
    ggsave(paste0("plot_silhouette_", p, ".png"), plot)
}
} # }
```
