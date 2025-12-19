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
# data("exdata_muscadet")

exdata_muscadet$ATAC
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681
exdata_muscadet$RNA
#> A muscomic object 
#>  type: RNA 
#>  label: scRNA-seq 
#>  cells: 69 
#>  counts: 69 cells x 300 features (genes)
#>  logratio: 69 cells x 212 features (genes)
#>  variant positions: 359

str(exdata_muscadet$ATAC, max.level = 2)
#> Formal class 'muscomic' [package "muscadet"] with 4 slots
#>   ..@ type      : chr "ATAC"
#>   ..@ label.omic: chr "scATAC-seq"
#>   ..@ coverage  :List of 2
#>   ..@ allelic   :List of 3
```
