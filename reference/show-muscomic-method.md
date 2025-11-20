# muscomic object overview

Overview of a
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
object.

## Usage

``` r
# S4 method for class 'muscomic'
show(object)
```

## Arguments

- object:

  A
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object.

## Value

Prints summary to
[`stdout`](https://rdrr.io/r/base/showConnections.html) and invisibly
returns `NULL`

## See also

[muscomic](https://icagen.github.io/muscadet/reference/muscomic-class.md),
[`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)

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
```
