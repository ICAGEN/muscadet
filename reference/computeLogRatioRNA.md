# Compute log R ratios for scRNA-seq data

Compute log R ratios from raw count matrices with method specifically
adapted for scRNA-seq data.

## Usage

``` r
computeLogRatioRNA(
  matTumor,
  matRef,
  genesCoord,
  genome = "hg38",
  genesPerWindow = 101,
  refReads = 100,
  refMeanReads = 0.01,
  thresh_capping = 3,
  all_steps = FALSE,
  quiet = FALSE
)
```

## Arguments

- matTumor:

  Raw count matrix *features x cells* for tumor/sample cells (`matrix`
  or `dgCMatrix`).

- matRef:

  Raw count matrix *features x cells* for reference cells (`matrix` or
  `dgCMatrix`).

- genesCoord:

  Data frame of gene coordinates with columns `CHROM`, `start`, `end`,
  `id` (`data.frame`).

- genome:

  Reference genome name among: "hg38", "hg19" and "mm10" (`character`).
  By default: "hg38".

- genesPerWindow:

  Number of genes per moving window (`integer` value). By default:
  `101`.

- refReads:

  Minimum of reads in reference cells (`integer` value). By default:
  `100`.

- refMeanReads:

  Minimum of average reads per reference cell (`integer` value). By
  default: `0.01`.

- thresh_capping:

  Threshold to cap the range of log R ratio values (`numeric` value). By
  default: `3`.

- all_steps:

  `TRUE` or `FALSE` (`logical`). Whether to keep intermediate result
  from every step in the final object. By default: `FALSE`.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

If `all_steps` is set to `FALSE`, a list containing:

- `matTumor`:

  Matrix of log R ratio values *cells x features* for the tumor/sample
  cells (`matrix`).

- `matRef`:

  Matrix of log R ratio values *cells x features* for the reference
  cells (`matrix`).

- `params`:

  List of parameters set for the `genesPerWindow`, `refMeans`,
  `refMeanReads` and `thresh_capping` arguments (`list`).

- `coord`:

  Data frame of coordinates for genes and associated data along the
  different steps (`data.frame`). Columns :

  - `CHROM`, `start`, `end`, `id`: coordinates and name of genes.

  - `sumReads.tum/ref`: sum of read counts for all cells in tumor or
    reference cells.

  - `meanReads.tum/ref`: mean of read counts per cells for tumor or
    reference cells.

  - `sdReads.tum/ref`: standard deviation of read counts per cells for
    tumor or reference cells.

  - `keep`: logical, `TRUE` for genes to keep after filtering based on
    reference coverage (depends on `refReads` and `refMeanReads`
    arguments).

  - `meanReads/sdReads.norm.tum/ref`: mean/sd of normalized counts per
    million for tumor/reference cells.

  - `meanLRR/sdReads.raw.tum/ref`: mean/sd of raw log R ratio (LRR) for
    tumor/reference cells.

  - `meanLRR/sdLRR.cap.tum/ref`: mean/sd of capped log R ratio (LRR) for
    tumor/reference cells (depends on `thresh_capping` argument).

  - `meanLRR/sdLRR.smoo.tum/ref`: mean/sd of smoothed log R ratio (LRR)
    for tumor/reference cells (means of moving windows defined by the
    `genesPerWindow` argument).

  - `meanLRR/sdLRR.cent.tum/ref`: mean/sd of centered log R ratio (LRR)
    for tumor/reference cells.

  - `meanLRR/sdLRR.corr.tum/ref`: mean/sd of final log R ratio (LRR)
    corrected by reference variability for tumor/reference cells.

## Details

The raw count matrix is transformed into log R ratios through the
following steps:

- Match genes in count matrix with coordinates

- Filtering on coverage (`refReads` and `refMeanReads` arguments)

- Normalization for sequencing depth

- Log transformation and normalization by reference data: log R ratio

- Capping the range of values (`thresh_capping` argument)

- Smoothing on genes windows

- Centering of cells

- Correcting by reference variability

## See also

Other computeLogRatio:
[`computeLogRatio()`](https://icagen.github.io/muscadet/reference/computeLogRatio.md),
[`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)

## Examples

``` r
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

# Compute log R ratio for RNA
obj_rna <- computeLogRatioRNA(
  matTumor = matCounts(muscadet)$RNA,
  matRef = matCounts(muscadet_ref)$RNA,
  genesCoord = coordFeatures(muscadet)$RNA,
  genome = muscadet$genome,
  refReads = 2 # low value for example subsampled datasets
)
#> Step 01 - Match genes in count matrix with coordinates
#> Step 02 - Filtering genes: Minimum of 2 read(s) in reference cells and minimum of 0.01 read(s) in average per reference cell
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - Smoothing values on gene windows: 101 genes per window
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
table(obj_rna$coord$keep)
#> 
#> FALSE  TRUE 
#>    88   212 

# With results form every step when `all_steps = TRUE`
obj_rna_all <- computeLogRatioRNA(
  matTumor = matCounts(muscadet)$RNA,
  matRef = matCounts(muscadet_ref)$RNA,
  genesCoord = coordFeatures(muscadet)$RNA,
  genome = muscadet$genome,
  refReads = 2, # low value for example subsampled datasets
  all_steps = TRUE
)
#> Step 01 - Match genes in count matrix with coordinates
#> Step 02 - Filtering genes: Minimum of 2 read(s) in reference cells and minimum of 0.01 read(s) in average per reference cell
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - Smoothing values on gene windows: 101 genes per window
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
names(obj_rna_all)
#>  [1] "step01" "step02" "step03" "step04" "step05" "step06" "step07" "step08"
#>  [9] "params" "coord" 
table(obj_rna_all$coord$keep)
#> 
#> FALSE  TRUE 
#>    88   212 
nrow(obj_rna_all$step08$matTumor)
#> [1] 69
```
