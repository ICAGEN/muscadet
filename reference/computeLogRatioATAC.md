# Compute log R ratios for scATAC-seq data

Compute log R ratios from raw count matrices with method specifically
adapted for scATAC-seq data. The counts per peaks are grouped into
counts per windows of peaks, thereby features become windows of peaks
instead of peaks.

## Usage

``` r
computeLogRatioATAC(
  matTumor,
  matRef,
  peaksCoord,
  genome = "hg38",
  windowSize = 1e+07,
  slidingSize = 2e+06,
  minReads = 5,
  minPeaks = 100,
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

- peaksCoord:

  Data frame of peak coordinates with columns `CHROM`, `start`, `end`,
  `id` (`data.frame`).

- genome:

  Reference genome name among: "hg38", "hg19" and "mm10" (`character`
  string). By default: "hg38".

- windowSize:

  Size of windows in base pairs (`integer` value). By default: `10e6`
  (10 Mbp).

- slidingSize:

  Distance between start positions of sliding windows in base pairs
  (`integer` value). If set to the same value as `windowSize`, the
  windows don't overlap. By default: `2e6` (2 Mbp).

- minReads:

  Minimum read average per window in reference cells (`integer` value).
  By default: `5`.

- minPeaks:

  Minimum number of peaks per window (`integer` value). By default:
  `100`.

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

  Matrix of log R ratio values *features x cells* for the tumor/sample
  cells (`matrix`).

- `matRef`:

  Matrix of log R ratio values *features x cells* for the reference
  cells (`matrix`).

- `params`:

  List of parameters set for the `windowSize`, `slidingSize`,
  `minReads`, `minPeaks` and `thresh_capping` arguments (`list`).

- `coord`:

  Data frame of coordinates for windows of peaks and associated data
  along the different steps (`data.frame`). Columns:

  - `CHROM`, `start`, `end`, `width`, `id`: coordinates and unique
    identifier of windows (depends on `windowSize` and `slidingSize`
    arguments).

  - `nPeaks`: number of peaks per window.

  - `sumReads.tum/ref`: sum of read counts for all cells in tumor or
    reference cells.

  - `meanReads.tum/ref`: mean of read counts per cells for tumor or
    reference cells.

  - `sdReads.tum/ref`: standard deviation of read counts per cells for
    tumor or reference cells.

  - `keep`: logical, `TRUE` for windows to keep after filtering based on
    coverage (depends on `minPeaks` and `meanReads` arguments).

  - `meanReads/sdReads.norm.tum/ref`: mean/sd of normalized counts per
    million for tumor/reference cells.

  - `meanLRR/sdReads.raw.tum/ref`: mean/sd of raw log R ratio (LRR) for
    tumor/reference cells.

  - `meanLRR/sdLRR.cap.tum/ref`: mean/sd of capped log R ratio (LRR) for
    tumor/reference cells (depends on `thresh_capping` argument).

  - `meanLRR/sdLRR.cent.tum/ref`: mean/sd of centered log R ratio (LRR)
    for tumor/reference cells.

  - `meanLRR/sdLRR.corr.tum/ref`: mean/sd of final log R ratio (LRR)
    corrected by reference variability for tumor/reference cells.

## Details

The raw count matrix is transformed into log R ratios through the
following steps:

- Group peaks in windows (`windowSize` and `slidingSize` arguments)

- Filtering on coverage (`minPeaks` and `meanReads` arguments)

- Normalization for sequencing depth

- Log transformation and normalization by reference data: log R ratio

- Capping the range of values (`thresh_capping` argument)

- Centering of cells

- Correcting by reference variability

## See also

Other computeLogRatio:
[`computeLogRatio()`](https://icagen.github.io/muscadet/reference/computeLogRatio.md),
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

# Compute log R ratio for ATAC
obj_atac <- computeLogRatioATAC(
  matTumor = matCounts(muscadet)$ATAC,
  matRef = matCounts(muscadet_ref)$ATAC,
  peaksCoord = coordFeatures(muscadet)$ATAC,
  genome = slot(muscadet, "genome"),
  minReads = 1, # low value for example subsampled datasets
  minPeaks = 1 # low value for example subsampled datasets
)
#> Step 01 - Group peaks in windows: window size set at 10 Mb, sliding by 2 Mb
#> Step 02 - Filtering windows: Minimum of 1 peaks per window with a minimum average of 1 read(s)
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - [No step 06 for scATAC-seq]
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
table(obj_atac$coord$keep)
#> 
#> FALSE  TRUE 
#>  1132   133 

# With results form every step when `all_steps = TRUE`
obj_atac_all <- computeLogRatioATAC(
  matTumor = matCounts(muscadet)$ATAC,
  matRef = matCounts(muscadet_ref)$ATAC,
  peaksCoord = coordFeatures(muscadet)$ATAC,
  genome = slot(muscadet, "genome"),
  minReads = 1, # low value for example subsampled datasets
  minPeaks = 1, # low value for example subsampled datasets
  all_steps = TRUE
)
#> Step 01 - Group peaks in windows: window size set at 10 Mb, sliding by 2 Mb
#> Step 02 - Filtering windows: Minimum of 1 peaks per window with a minimum average of 1 read(s)
#> Step 03 - Normalization for sequencing depth: Normalized counts per million
#> Step 04 - Log transformation and normalization by reference data: log R ratio
#> Step 05 - Capping the range of values: threshold = 3
#> Step 06 - [No step 06 for scATAC-seq]
#> Step 07 - Centering of cells
#> Step 08 - Correcting by reference variability
names(obj_atac_all)
#> [1] "step01" "step02" "step03" "step04" "step05" "step07" "step08" "params"
#> [9] "coord" 
table(obj_atac_all$coord$keep)
#> 
#> FALSE  TRUE 
#>  1132   133 
nrow(obj_atac_all$step08$matTumor)
#> [1] 133
```
