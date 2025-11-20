# The muscomic class

A `muscomic` object encapsulates coverage and allelic data from one
single-cell omic dataset used as primary input for muscadet analysis.

## Slots

- `type`:

  Type of single-cell omic. Either ATAC or RNA currently supported
  (`character`).

- `label.omic`:

  Label to display for the single-cell omic (`character`).

- `coverage`:

  Coverage data on features (`list`).

- `allelic`:

  Allelic data at variant positions (common SNPs or individual-specific
  heterozygous positions) (`list`).

## See also

Functions related to `muscomic` objects:

- [`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)
  (create `muscomic` objects)

- [muscadet-access](https://icagen.github.io/muscadet/reference/musc-access.md)
  (simplified access and assignment methods)

- [.DollarNames.muscadet](https://icagen.github.io/muscadet/reference/musc-auto.md)
  (autocompletion for `$`)

- [`show,muscomic-method`](https://icagen.github.io/muscadet/reference/show-muscomic-method.md)
  (`show` method)

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
muscadet_obj$RNA
#> A muscomic object of type RNA labelled scRNA-seq containing: 
#>  log.ratio coverage data matrix 
#>  119 cells 
#>  349 features: genes 
#>  373 variant positions 

str(muscadet_obj$ATAC, max.level = 2)
#> Formal class 'muscomic' [package "muscadet"] with 4 slots
#>   ..@ type      : chr "ATAC"
#>   ..@ label.omic: chr "scATAC-seq"
#>   ..@ coverage  :List of 6
#>   ..@ allelic   :List of 1
```
