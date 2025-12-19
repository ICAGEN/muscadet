# Autocompletion for `$` access on `muscadet` or `muscomic` objects

Enable autocompletion for `$` access for
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
or
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
objects. For `muscadet` objects, it also lists omic datasets contained
inside the `omics` slot.

## Usage

``` r
# S3 method for class 'muscadet'
.DollarNames(x, pattern = "")

# S3 method for class 'muscomic'
.DollarNames(x, pattern = "")
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  or
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object.

- pattern:

  A regular expression. Only matching names are returned.

## Value

Character vector of matching element names.

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")
exdata_muscadet$ATAC
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681

# Load example muscadet object
# data("exdata_muscadet")
exdata_muscadet$ATAC$type
#> [1] "ATAC"
```
