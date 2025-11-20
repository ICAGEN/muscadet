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
# data("muscadet_obj")
muscadet_obj$ATAC
#> A muscomic object of type ATAC labelled scATAC-seq containing: 
#>  log.ratio coverage data matrix 
#>  112 cells 
#>  133 features: windows of peaks 
#>  691 variant positions 

# Load example muscadet object
# data("muscadet_obj")
muscadet_obj$ATAC$type
#> [1] "ATAC"
```
