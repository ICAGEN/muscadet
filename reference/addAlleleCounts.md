# Add allele counts to a `muscadet` object

This function adds allele counts data to the `omics` of a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. The data frames in the `allele_counts` list are assigned to the
`allelic` slots of `omics`.

## Usage

``` r
addAlleleCounts(x, allele_counts)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object.

- allele_counts:

  A list of data frames where each data frame contains allele counts for
  a specific omic (`list`). The list must have the same length and order
  as the number of omics in the `x` object. Each data frames must
  contain the following columns : `cell`, `id`, `CHROM`, `POS`, `REF`,
  `ALT`, `RD`, `AD`, `DP`, `GT`. See
  [exdata_allele_counts](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  for details.

## Value

A modified
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object with updated allele counts in the `allelic` slot of each
`muscomic` object in the `omics` slot.

## Note

As the allele counts are not used for computing log R ratios and cell
clustering, they are not mandatory at the creation of `muscomic` and
`muscadet` objects. The allele counts data can be added to objects later
with this `addAlleleCounts` function, before using the
[`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)
function.

This function is also useful to add allele counts for
individual-specific variant positions to a common `muscadet` object, for
example for the reference `muscadet` object: a common `muscadet` object
with reference coverage data can be stored as a unique object to use as
reference for computing log R ratios for all samples, and then it can
updated with allele counts at individual-specific variant positions
(e.g. found by bulk sequencing) before Copy Number Alterations (CNAs)
calling.

## See also

[`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md),
[`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")
# data("exdata_muscadet_ref")

# Add allele counts data frames to muscadet objects
exdata_muscadet <- addAlleleCounts(
    exdata_muscadet,
    allele_counts = list(exdata_allele_counts_atac_tumor, exdata_allele_counts_rna_tumor))
exdata_muscadet_ref <- addAlleleCounts(
    exdata_muscadet_ref,
    allele_counts = list(exdata_allele_counts_atac_ref, exdata_allele_counts_rna_ref))
```
