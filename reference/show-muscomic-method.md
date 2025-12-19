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
```
