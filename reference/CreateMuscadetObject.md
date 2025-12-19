# Create a `muscadet` object

Create a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Usage

``` r
CreateMuscadetObject(
  omics,
  bulk.lrr = NULL,
  bulk.label = NULL,
  genome = "hg38"
)
```

## Arguments

- omics:

  A list of
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  objects (`list`). The names of the list will set the names of omics in
  the final object, if the list is unnamed, the type is taken instead.

- bulk.lrr:

  A data frame containing log R ratio per genomic segments from bulk
  sequencing data (`data.frame`). One row per segment and 4 columns
  ordered as followed: chromosome (`character`), start position
  (`integer`), end position (`integer`), and Log R ratio value
  (`numeric`).

- bulk.label:

  Label for bulk data (`character` string).

- genome:

  Reference genome name among: "hg38", "hg19" and "mm10" (`character`
  string). "hg38" by default.

## Value

A
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Note

A
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object can contain several
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
objects of the same type (`ATAC` or `RNA` for slot `type`) but they
can't have identical labels (slot `label.omic`).

## See also

- [`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)

- [muscomic](https://icagen.github.io/muscadet/reference/muscomic-class.md),
  [muscadet](https://icagen.github.io/muscadet/reference/muscadet-class.md)

## Examples

``` r
# Create muscomic objects
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = exdata_mat_counts_atac_tumor,
  allele_counts = exdata_allele_counts_atac_tumor,
  features = exdata_peaks
)
rna <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = exdata_mat_counts_rna_tumor,
  allele_counts = exdata_allele_counts_rna_tumor,
  features = exdata_genes
)
atac_ref <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = exdata_mat_counts_atac_ref,
  allele_counts = exdata_allele_counts_atac_ref,
  features = exdata_peaks
)
rna_ref <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = exdata_mat_counts_rna_ref,
  allele_counts = exdata_allele_counts_rna_ref,
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
```
