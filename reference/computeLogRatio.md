# Compute log R ratios

Computes log R ratios from raw count matrices within
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects. The log R ratios values are computed based on read counts from
a sample
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object *versus* read counts a the `reference`
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. In the output sample
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object, the newly computed matrix of log R ratios is added for the
selected `omic`.

## Usage

``` r
computeLogRatio(
  x,
  reference,
  omic,
  method = NULL,
  new.label.features = NULL,
  remove.raw = TRUE,
  quiet = FALSE,
  all_steps = FALSE,
  ...
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing sample data (`muscadet`).

- reference:

  Another
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing reference data (`muscadet`).

- omic:

  Name of the omic to apply this function (`character` string).

- method:

  Method to apply to the selected omic (`character` string). Supported
  methods are "ATAC" and "RNA":

  - "ATAC" method calls for
    [`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)
    function

  - "RNA" method calls for
    [`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)
    function

  If `NULL` and if the omic type is either "ATAC" or "RNA", the
  corresponding method will be applied.

- new.label.features:

  New label for features (`character` string). If `NULL`, the label
  remains unchanged when using the "RNA" method and becomes "windows of
  peaks" when using the "ATAC" method.

- remove.raw:

  `TRUE` or `FALSE` (`logical`). Whether to remove raw count matrices.
  `TRUE` by default to reduce object size. Setting it to `FALSE` will
  keep raw count matrices within the object.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

- all_steps:

  Logical value to indicate whether matrices at each step of the log
  ratio computation are kept and returned (`logical`). IMPORTANT: if
  `TRUE`, the function does not return an updated
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object but a list with matrices at each step of the log ratio
  computation (see Value section). Default is `FALSE`.

- ...:

  Arguments passed on to
  [`computeLogRatioATAC`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md),
  [`computeLogRatioRNA`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)

  `windowSize`

  :   Size of windows in base pairs (`integer` value). By default:
      `10e6` (10 Mbp).

  `slidingSize`

  :   Distance between start positions of sliding windows in base pairs
      (`integer` value). If set to the same value as `windowSize`, the
      windows don't overlap. By default: `2e6` (2 Mbp).

  `minReads`

  :   Minimum read average per window in reference cells (`integer`
      value). By default: `5`.

  `minPeaks`

  :   Minimum number of peaks per window (`integer` value). By default:
      `100`.

  `genesPerWindow`

  :   Number of genes per moving window (`integer` value). By default:
      `101`.

  `refReads`

  :   Minimum of reads in reference cells (`integer` value). By default:
      `100`.

  `refMeanReads`

  :   Minimum of average reads per reference cell (`integer` value). By
      default: `0.01`.

  `thresh_capping`

  :   Threshold to cap the range of log R ratio values (`numeric`
      value). By default: `3`.

## Value

A
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object corresponding to the sample
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object (`x`) containing the computed log R ratio matrix in the
`coverage` slot of the selected `omic`.

If the `all_steps` argument is set to `TRUE`, it returns a list with
intermediate matrices `matTumor` and `matRef` from every step from
`step01` to `step08` (no `step06` for `method` = "ATAC"), `coord` table
of features coordinates and associated data, and `params` list of
parameters used for the function.

## Details

Log R ratios computation steps are described in functions:

- [`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)

- [`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)

## See also

Other computeLogRatio:
[`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md),
[`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)

## Examples

``` r
# Create muscomic objects
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor,
  allele_counts = allele_counts_atac_tumor,
  features = peaks
)
rna <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_tumor,
  allele_counts = allele_counts_rna_tumor,
  features = genes
)
atac_ref <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_ref,
  allele_counts = allele_counts_atac_ref,
  features = peaks
)
rna_ref <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_ref,
  allele_counts = allele_counts_rna_ref,
  features = genes
)

# Create muscadet objects
muscadet <- CreateMuscadetObject(
  omics = list(atac, rna),
  bulk.lrr = bulk_lrr,
  bulk.label = "WGS",
  genome = "hg38"
)
muscadet_ref <- CreateMuscadetObject(
  omics = list(atac_ref, rna_ref),
  genome = "hg38"
)

# compute log R ratios for ATAC
muscadet <- computeLogRatio(
  x = muscadet,
  reference = muscadet_ref,
  omic = "ATAC",
  method = "ATAC",
  minReads = 1, # low value for example subsampled datasets
  minPeaks = 1 # low value for example subsampled datasets
)
#> -- computeLogRatio: Method 'ATAC' using computeLogRatioATAC().
#> Step 01 - Group peaks in windows: window size set at 10 Mb, sliding by 2 Mb
#> Step 02 - Filtering windows: Minimum of 1 peaks per window with a minimum average of 1 read(s)
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - [No step 06 for scATAC-seq]
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
#> Done.

# compute log R ratios for RNA
muscadet <- computeLogRatio(
  x = muscadet,
  reference = muscadet_ref,
  omic = "RNA",
  method = "RNA",
  refReads = 2 # low value for example subsampled datasets
)
#> -- computeLogRatio: Method 'RNA' using computeLogRatioRNA().
#> Step 01 - Match genes in count matrix with coordinates
#> Step 02 - Filtering genes: Minimum of 2 read(s) in reference cells and minimum of 0.01 read(s) in average per reference cell
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - Smoothing values on gene windows: 101 genes per window
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
#> Done.
```
