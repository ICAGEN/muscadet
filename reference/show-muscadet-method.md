# muscadet object overview

Overview of a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object.

## Usage

``` r
# S4 method for class 'muscadet'
show(object)
```

## Arguments

- object:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object.

## Value

Prints summary to
[`stdout`](https://rdrr.io/r/base/showConnections.html) and invisibly
returns `NULL`

## See also

[muscadet](https://icagen.github.io/muscadet/reference/muscadet-class.md),
[`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)

## Examples

``` r
# Load example muscadet object
# data("muscadet_obj")

# Overview of the muscadet object
show(muscadet_obj)
#> A muscadet object 
#>  2 omics: ATAC, RNA 
#>  types: ATAC, RNA 
#>  labels: scATAC-seq, scRNA-seq 
#>  coverage data matrix: log.ratio, log.ratio 
#>  cells: 112, 119 (common: 84, total: 147) 
#>  features: 133, 349 
#>  feature labels: windows of peaks, genes 
#>  variant positions: 691, 373 
#>  data from paired bulk sequencing: WGS 
#>  clustering: partitions = 0.6, 0.8, 1 ; optimal partition = 1 
#>  CNA calling: 2 clusters ; 47 consensus segments including 1 CNA segments 
#>  genome: hg38 

# Overview of the muscomic objects within
show(slot(muscadet_obj, "omics"))
#> $ATAC
#> A muscomic object of type ATAC labelled scATAC-seq containing: 
#>  log.ratio coverage data matrix 
#>  112 cells 
#>  133 features: windows of peaks 
#>  691 variant positions 
#> 
#> $RNA
#> A muscomic object of type RNA labelled scRNA-seq containing: 
#>  log.ratio coverage data matrix 
#>  119 cells 
#>  349 features: genes 
#>  373 variant positions 
#> 

# Overview of the first muscomic objects within
show(slot(muscadet_obj, "omics")[[1]])
#> A muscomic object of type ATAC labelled scATAC-seq containing: 
#>  log.ratio coverage data matrix 
#>  112 cells 
#>  133 features: windows of peaks 
#>  691 variant positions 
```
