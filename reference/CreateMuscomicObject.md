# Create a `muscomic` object

Create a
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
object containing coverage and (optionally) allelic count information
for a single-cell omic dataset.

## Usage

``` r
CreateMuscomicObject(
  type = c("ATAC", "RNA"),
  mat_counts,
  features,
  allele_counts = NULL,
  label.omic = NULL,
  label.features = NULL
)
```

## Arguments

- type:

  Type of single cell omic, either "ATAC" or "RNA" (`character` string).

- mat_counts:

  Matrix of raw counts *cells x features* (`matrix` or `dgCMatrix`).
  Rows are cells and columns are features. (they must correspond to the
  id column of `features`)

- features:

  Data frame of features (peaks, genes...) coordinates on genome
  (`data.frame`). It should contain 4 columns:

  `CHROM`

  :   Chromosome names in character format, e.g. "15", "X"
      (`character`).

  `start`

  :   Start positions (`integer`).

  `end`

  :   End positions (`integer`).

  `id`

  :   Unique identifiers, e.g. gene name "CDH1" or peak identifier
      CHROM_start_end "1_1600338_1600838" (`character`). It should match
      the feature identifiers as row names of `mat_counts`.

- allele_counts:

  Data frame of allele counts at variant positions per cell
  (`data.frame`). Variant positions can be either common single
  nucleotide polymorphisms (SNPs) positions or individual-specific
  heterozygous positions retrieved from bulk sequencing. The data frame
  format is based on the Variant Calling Format (VCF), thereby it must
  contain the following columns : `cell`, `CHROM`, `POS`, `REF`, `ALT`,
  `RD`, `AD`. See
  [`exdata_allele_counts()`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  for details.

- label.omic:

  Label for the single cell omic (`character` string). By default
  "scATAC-seq" is used for "ATAC" type and "scRNA-seq" for "RNA" type.

- label.features:

  Label for features (`character` string). By default "peaks" is used
  for "ATAC" type and "genes" for "RNA" type.

## Value

A
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
object.

## See also

- [`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)

- [muscomic](https://icagen.github.io/muscadet/reference/muscomic-class.md),
  [muscadet](https://icagen.github.io/muscadet/reference/muscadet-class.md)

## Examples

``` r
# Minimal example
mat <- Matrix::sparseMatrix(
    i = c(1, 1, 2),
    j = c(1, 2, 1),
    x = c(5, 3, 2),
    dims = c(2, 2),
    dimnames = list(c("cell1", "cell2"), c("f1", "f2"))
)

features <- data.frame(
    CHROM = c("1", "1"),
    start = c(100, 200),
    end   = c(150, 250),
    id    = c("f1", "f2")
)

allele <- data.frame(
    cell = c("cell1", "cell2"),
    CHROM = c("1", "1"),
    POS = c(100, 200),
    REF = c("A", "G"),
    ALT = c("C", "T"),
    RD = c(10, 5),
    AD = c(3, 1)
)

muscomic <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat,
    features = features,
    allele_counts = allele
)

# On the example dataset
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = exdata_mat_counts_atac_tumor,
  allele_counts = exdata_allele_counts_atac_tumor,
  features = exdata_peaks
)
atac
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: None
#>  variant positions: 681

# without allele counts data (not required for clustering step)
atac2 <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = exdata_mat_counts_atac_tumor,
  features = exdata_peaks
)
atac2
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: None
#>  variant positions: 0
```
