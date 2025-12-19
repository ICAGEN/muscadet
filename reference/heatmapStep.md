# Create heatmap and distribution plots of the different steps of computing log R ratios

This function generates heatmap and distribution plots for tumor and
reference cells for any step of computing log R ratios matrices. The
input object corresponds to the output of
[`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)
or
[`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)
with the argument `all_steps = TRUE`.

## Usage

``` r
heatmapStep(
  obj,
  step,
  filename,
  title = NULL,
  col_quantiles = c(0.1, 0.4, 0.6, 0.9),
  col_breaks = NULL,
  colors = c("#00008E", "white", "white", "#630000")
)
```

## Arguments

- obj:

  A list provided as output by the
  [`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)
  or
  [`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)
  functions with the argument `all_steps = TRUE`. It includes tumor and
  reference matrices at each step of the computing of log R ratio
  matrices.

- step:

  The step within the `obj` list to use for plotting (`character`
  string). It must match one of the names in `obj`.

- filename:

  File path to save the output plot (`character` string). The file
  format is inferred from the extension (".png", ".pdf" or ".svg").

- title:

  Title of the plot (`character` string). If `NULL`, the title is
  automatically generated using the provided step argument and its
  corresponding value name (`obj[[step]]$name`).

- col_quantiles:

  A numeric vector of length 4, specifying the quantiles to use for the
  color breakpoints in the heatmap (`numeric`). Either `col_quantiles`
  or `col_breaks` must be provided, if both are provided `col_breaks` is
  used. Default is `c(0.1, 0.4, 0.6, 0.9)`.

- col_breaks:

  A numeric vector of length 4, specifying custom breakpoints for the
  color scale in the heatmap (`numeric`). Either `col_quantiles` or
  `col_breaks` must be provided, if both are provided `col_breaks` is
  used. Default is `NULL`.

- colors:

  A character vector of 4 colors used for the color scale of the heatmap
  (`character` vector). Default is
  `c("#00008E", "white", "white", "#630000")`.

## Value

The function does not return any value but saves a heatmaps-histograms
plot to the specified file.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create muscomic objects
atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = exdata_mat_counts_atac_tumor,
    features = exdata_peaks
)
rna <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = exdata_mat_counts_rna_tumor,
    features = exdata_genes
)
atac_ref <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = exdata_mat_counts_atac_ref,
    features = exdata_peaks
)
rna_ref <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = exdata_mat_counts_rna_ref,
    features = exdata_genes
)

# Create muscadet objects
muscadet <- CreateMuscadetObject(
    omics = list(atac, rna),
    bulk.lrr = exdata_bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38"
)
muscadet_ref <- CreateMuscadetObject(
    omics = list(atac_ref, rna_ref),
    genome = "hg38"
)

# Compute log R ratios with `all_steps = TRUE`
obj_atac_all <- computeLogRatioATAC(
    matTumor = matCounts(muscadet)$ATAC,
    matRef = matCounts(muscadet_ref)$ATAC,
    peaksCoord = coordFeatures(muscadet)$ATAC,
    genome = muscadet$genome,
    minReads = 1, # low value for example subsampled datasets
    minPeaks = 1, # low value for example subsampled datasets
    all_steps = TRUE
)
names(obj_atac_all)

# Plot heatmap and distribution of values for Step01
heatmapStep(obj = obj_atac_all,
            step = "step01",
            filename = file.path(tempdir(), "step01.png"),
            title = "Example dataset - Step 01")

# Plot heatmap and distribution of values for all steps
for (step in grep("step", names(obj_atac_all), value = TRUE)) {
    heatmapStep(
        obj_atac_all,
        step,
        filename = file.path(tempdir(), paste0("ATAC_", step, ".pdf")),
        title = paste("ATAC -", step)
    )
}
} # }
```
