# Create a `muscomic` object

Create a
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
object.

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

  Matrix of raw counts *features x cells* (`matrix` or `dgCMatrix`).
  Rows are features (they must correspond to the id column of
  `features`), and columns are cells.

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
  contain the following columns : `cell`, `id`, `CHROM`, `POS`, `REF`,
  `ALT`, `RD`, `AD`, `DP`, (`GT`). See
  [`allele_counts()`](https://icagen.github.io/muscadet/reference/allele_counts.md)
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
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor,
  allele_counts = allele_counts_atac_tumor,
  features = peaks
)
atac
#> A muscomic object of type ATAC labelled scATAC-seq containing: 
#>  mat.counts coverage data matrix 
#>  112 cells 
#>  1000 features: peaks 
#>  691 variant positions 

rna <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_tumor,
  allele_counts = allele_counts_rna_tumor,
  features = genes
)
rna
#> A muscomic object of type RNA labelled scRNA-seq containing: 
#>  mat.counts coverage data matrix 
#>  119 cells 
#>  500 features: genes 
#>  373 variant positions 

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
rna_ref
#> A muscomic object of type RNA labelled scRNA-seq containing: 
#>  mat.counts coverage data matrix 
#>  97 cells 
#>  500 features: genes 
#>  373 variant positions 

# without allele counts data (not required for clustering step)
atac2 <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor,
  features = peaks
)
atac2
#> A muscomic object of type ATAC labelled scATAC-seq containing: 
#>  mat.counts coverage data matrix 
#>  112 cells 
#>  1000 features: peaks 
#>  0 variant positions 
```
