# Plot UMAP Coordinates from a muscadet Object

Visualize UMAP (Uniform Manifold Approximation and Projection)
coordinates stored in a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. The function allows coloring by clusters and optionally adding
cluster labels.

## Usage

``` r
plotUMAP(
  x,
  partition = NULL,
  clusters = NULL,
  colors = NULL,
  title = "",
  lab.x = "UMAP 1",
  lab.y = "UMAP 2",
  add_clusters_labels = FALSE,
  point.size = 0.5,
  legend.point.size = 3,
  ...
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing clustering data (using
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).

- partition:

  (Optional) Value specifying the clustering partition to plot
  (`numeric` or `character`). It should be either the resolution or the
  number of cluster (k) used for clustering depending on the clustering
  method (`res_range` or `k_range` with
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)).
  If both `partition` and `clusters` arguments are `NULL` (default), the
  assigned clusters for CNA calling (`x@cnacalling$cluster`) are used if
  available in `x` (see
  [`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md)).

- clusters:

  (Optional) A custom named vector of cluster assignments (`integer`
  named vector). Names must corresponds to cell names within the
  muscadet object `x`. If it contains less cells than the muscadet
  object `x`, the missing cells are filtered out and not displayed in
  the heatmap. If `show_missing = FALSE` only the provided cells with
  data in all omics will be displayed. If both `partition` and
  `clusters` arguments are `NULL` (default), the assigned clusters for
  CNA calling (`x@cnacalling$cluster`) are used if available in `x` (see
  [`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md)).

- colors:

  Vector of colors for the cluster annotation (`character` vector). If
  `NULL` (default), it uses predefined colors.

- title:

  Character string for the title of the plot (`character` string).
  Default is an empty character string.

- lab.x:

  Label for the x-axis (`character` string). Default is "UMAP 1".

- lab.y:

  Label for the y-axis (`character` string). Default is "UMAP 2".

- add_clusters_labels:

  Logical. If `TRUE`, adds the cluster names as text or boxed labels
  using
  [`add_labels()`](https://icagen.github.io/muscadet/reference/add_labels.md).
  Default is `FALSE`.

- point.size:

  Numeric. Size of the points in the UMAP plot passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default is `0.5`.

- legend.point.size:

  Numeric. Size of the points in the legend of the UMAP plot. Default is
  `3`.

- ...:

  Additional arguments passed to
  [`add_labels()`](https://icagen.github.io/muscadet/reference/add_labels.md)
  providing an underlying geom for label names
  ([`geom_text`](https://ggplot2.tidyverse.org/reference/geom_text.html),
  [`geom_label`](https://ggplot2.tidyverse.org/reference/geom_text.html),
  [`geom_text_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html),
  or
  [`geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)).

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  p <- plotUMAP(muscadet_obj,
                partition = 0.6,
                title = "UMAP copy-number clusters")
  print(p)
} # }
```
