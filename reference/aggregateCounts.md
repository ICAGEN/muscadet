# Aggregate counts for `muscadet` objects

This function combines allelic (counts at variant positions, either
common SNPs or individual-specific heterozygous positions) and coverage
counts (counts on features) from all omics per cluster for both sample
and reference. The resulting merged data is stored in the `cnacalling`
slot of the sample
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Usage

``` r
aggregateCounts(x, reference, nor.het = TRUE, quiet = FALSE)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing sample data (`muscadet`). This object must include
  clustering assignments in the `cnacalling$clusters` slot.

- reference:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing reference data (`muscadet`).

- nor.het:

  A logical value to specify if normal reference allele counts are
  modified to: total normal depth counts divided by 2, to force these
  positions to be heterozygous in the normal reference in allelic data
  (e.g. when reference cells are from a different individual but
  individual-specific heterozygous positions based on matched bulk
  sequencing data are used and thereby assumed to be heterozygous)
  before combining coverage and allelic data. Default is `TRUE`.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

A modified
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object corresponding to the `x`
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object, with updated `cnacalling` slot containing:

- `allelic.counts`: Processed allelic counts on variant positions, for
  all omics.

- `coverage.counts`: Processed coverage counts merged with the
  reference.

- `combined.counts`: Combined data for allelic and coverage counts.

Abbreviations:

- RD = Reference allele read depth

- AD = Alternative allele read depth

- DP = Total read depth

- TUM = tumor sample

- NOR = normal reference

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")
# data("exdata_muscadet_ref")

# Aggregate counts from all omics from both sample and reference
exdata_muscadet <- aggregateCounts(exdata_muscadet, exdata_muscadet_ref)
#> Clusters used: 1 (41 cells), 2 (36 cells)
#> Allelic data processing...
#> Coverage data processing...
#> Combining allelic and coverage data...
```
