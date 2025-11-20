# Heatmap plot for `muscadet` object

This function generates a heatmap to visualize log R ratio (LRR) data
contained in
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects. One heatmap is generated per omic, rows are cells and columns
are chromosomes, for
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object containing multiple omics, the heatmaps are plotted horizontally
aligned. The cells can be clustered for a specific clustering partition
following the clustering step of
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object, or custom cluster assignments can be used. Additionally, LRR
values from bulk sequencing data can be plotted as an annotation under
the heatmaps.

## Usage

``` r
heatmapMuscadet(
  x,
  filename = NULL,
  partition = NULL,
  clusters = NULL,
  add_bulk_lrr = TRUE,
  show_missing = TRUE,
  averages = FALSE,
  title = "",
  row_annots = NULL,
  white_scale = c(0.3, 0.7),
  colors = NULL,
  png_res = 300,
  quiet = FALSE
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing LRR data for all omics (with
  [`computeLogRatio()`](https://icagen.github.io/muscadet/reference/computeLogRatio.md))
  and clustering data (with
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md))
  (`muscadet`).

- filename:

  (Optional) Character string specifying the file path to save the
  heatmap image in the PNG (if it ends by .png), PDF (if it ends by
  .pdf) or SVG (if it ends by .svg) format (`character` string).

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

- add_bulk_lrr:

  Logical. If `TRUE` (default), adds bulk log R ratio (LRR) data as
  annotation if available in the muscadet object.

- show_missing:

  Logical. If `TRUE` (default), missing cells (i.e., cells with missing
  data in at least one omic) are displayed in the heatmaps.

- averages:

  Logical. If `TRUE`, plots the average log R ratio per cluster. Default
  is `FALSE`.

- title:

  Character string for the title of the plot (`character` string).
  Default is an empty character string.

- row_annots:

  Optional. A list of
  [`HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation-class.html)
  objects from the
  [`ComplexHeatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/ComplexHeatmap-package.html)
  package, specifying row annotations to add on left part of the
  heatmap. Each element in the list must be of class
  [`HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation-class.html),
  must be a row annotation (using
  [`rowAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/rowAnnotation.html)
  or
  [`HeatmapAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
  with `which = 'row'`), and must have a unique name (`name` argument in
  [`rowAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/rowAnnotation.html)
  or
  [`HeatmapAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)).
  If `averages =FALSE`, annotations must concern cells, while if
  `averages = TRUE` they must concern clusters. Default is `NULL`, no
  row annotations is added.

- white_scale:

  Numeric vector of length 2 or a list of numeric vectors (`numeric`
  vector or `list`).

  - If a numeric vector of length 2, the same white color boundaries are
    applied to all omics in the muscadet object. E.g. `c(0.3, 0.7)`
    (default).

  - If a list (named or not with omics name), it must have the same
    length as the number of omics in the muscadet object, where each
    vector element applies the white color boundaries for a specific
    omic. E.g. `list(c(0.3, 0.7), c(0.4, 0.6))` uses 0.3 and 0.7
    quantiles of LRR ref data for the 1st omic heatmap, and 0.4 and 0.6
    quantiles for the second.

    Values of the vectors must be between 0 and 1, specifying the
    quantiles of the LRR reference data that define the boundaries for
    the white color in each heatmap. LRR values falling within this
    range are considered close to the majority of the LRR reference
    data, indicating no significant gain or loss of coverage, and are
    represented as white on the heatmap. Default is `c(0.3, 0.7)`.

- colors:

  Vector of colors for the cluster annotation (`character` vector). If
  `NULL` (default), it uses predefined colors.

- png_res:

  Resolution in ppi for
  [`grDevices::png()`](https://rdrr.io/r/grDevices/png.html) if
  `filename` ends with the .png extension (`numeric`). Default is `300`.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

A list containing:

- `plot`: A [`gTree`](https://rdrr.io/r/grid/grid-defunct.html) object
  created with
  [`grid::grid.grab()`](https://rdrr.io/r/grid/grid.grab.html)
  ([`gTree`](https://rdrr.io/r/grid/grid-defunct.html)).

- `width`: Width of the heatmap plot in mm
  ([`unit`](https://rdrr.io/r/grid/unit.html)).

- `height`: Height of the heatmap plot in mm
  ([`unit`](https://rdrr.io/r/grid/unit.html)).

If the `filename` argument is provided, the heatmap is directly saved as
a PNG image at the provided path.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example muscadet object
# data("muscadet_obj")


# --- Method "seurat" ---

print(muscadet_obj$clustering$params$method)

# Perform clustering if not already done
# muscadet_obj <- clusterMuscadet(
#     x = muscadet_obj,
#     method = "seurat",
#     res_range = c(0.6, 0.8),
#     dims_list = list(1:10, 1:10),
#     knn_seurat = 10, # adapted for low number of cells in example data
#     knn_range_seurat = 30 # adapted for low number of cells in example data
# )

# Plot a single partition
heatmapMuscadet(muscadet_obj,
                filename = file.path("heatmap_res0.6.png"),
                partition = 0.6,
                show_missing = FALSE) # only displaying cells without missing data

ht <- heatmapMuscadet(muscadet_obj, partition = 0.6)
pdf(
    file = file.path("heatmap_res0.6.pdf"),
    width = ht$width * 0.0393701, # convert to inches
    height = ht$height * 0.0393701, # convert to inches
)
grid.draw(ht$plot)
dev.off()


# Loop over partitions
for (p in names(muscadet_obj$clustering$clusters)) {
    filename <- paste0("heatmap_res", p, ".png")
    title <- paste(
        "Example |",
        paste0("method=", muscadet_obj$clustering$params[["method"]]), "|",
        paste0("omics=", paste0(muscadet_obj$clustering$params[["omics"]], collapse = ",")), "|",
        paste0("dims=", "1:10,1:10"), "|",
        paste0("res=", p)
    )
    heatmapMuscadet(muscadet_obj, filename, partition = p, title = title)
}

# --- Plot Averages per Clusters ---

heatmapMuscadet(muscadet_obj,
                filename = file.path("heatmap_res0.6_averages.png"),
                partition = 0.6,
                averages = TRUE,
                add_bulk_lrr = FALSE)

# --- Add Row Annotation ---

library("ComplexHeatmap")
library("grid")

# Define example cell annotation
muscadet_cells <- Reduce(union, SeuratObject::Cells(muscadet_obj))
cells_origin <- setNames(c(
    rep("sample1", ceiling(length(muscadet_cells) / 2)),
    rep("sample2", floor(length(muscadet_cells) / 2))
    ),
    muscadet_cells
)
cells_origin <- cells_origin[sort(names(cells_origin))]
# IMPORTANT: annotation names (cells) must be sorted to match heatmap
# matrices (column names of log ratio matrix are sorted in heatmapMuscadet())

# Create row annotation
ha <- rowAnnotation(
    annot = anno_simple(
        cells_origin[sort(names(cells_origin))],
        col = c("sample1" = "cadetblue3", "sample2" = "orchid3")),
    name = "origin", # unique name
    annotation_label = "origin", # label displayed on heatmap
    annotation_name_gp = gpar(fontsize = 10) # change font size
)

# Plot heatmap with supplementary row annotation
heatmapMuscadet(muscadet_obj,
                filename = file.path("heatmap_res0.6_annot.png"),
                partition = 0.6,
                row_annots = list(ha))

# --- Add Row Annotation for averages ---

library("ComplexHeatmap")
library("grid")

# Define example cluster annotation
clus <- setNames(c("annot1", "annot2"), c(1, 2)) # 2 clusters for partition 0.6
clus <- clus[order(names(clus))]
# IMPORTANT: annotation names (clusters) must be sorted to match heatmap
# matrices (column names of log ratio averages matrix are sorted in heatmapMuscadet())


# Create row annotation
ha2 <- rowAnnotation(
    annot = anno_simple(clus, col = c("annot1" = "tomato", "annot2" = "gold2")),
    name = "annot", # unique name
    annotation_label = "annot", # label displayed on heatmap
    annotation_name_gp = gpar(fontsize = 10) # change font size
)
heatmapMuscadet(muscadet_obj,
                averages = TRUE,
                filename = file.path("heatmap_res0.6_annot_averages.png"),
                partition = 0.6,
                row_annots = list(ha2))


# --- Method "hclust" ---

# Perform clustering if not already done
muscadet_obj2 <- clusterMuscadet(
    x = muscadet_obj,
    method = "hclust",
    k_range = 3:5,
    dist_method = "euclidean",
    hclust_method = "ward.D"
)

print(muscadet_obj2$clustering$params$method)

# Plot a single partition
heatmapMuscadet(muscadet_obj2,
                filename = file.path("heatmap_k3.png"),
                partition = 3,
                show_missing = FALSE)

# Loop over partitions
for (p in names(muscadet_obj2$clustering$clusters)) {

    filename <- paste0("heatmap_k", p, ".png")
    title <- paste(
        "Example |",
        paste0("method=", muscadet_obj2$clustering$params[["method"]]), "|",
        muscadet_obj2$clustering$params[["dist_method"]],
        muscadet_obj2$clustering$params[["hclust_method"]], "|",
        paste0("weights=",
               paste0(muscadet_obj2$clustering$params[["weights"]], collapse = ",")),
        "|",
        paste0("k=", p)
    )

    heatmapMuscadet(muscadet_obj2, filename, partition = p, title = title)
}
} # }
```
