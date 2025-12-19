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
# data("exdata_muscadet")

# Overview of the muscadet object
show(exdata_muscadet)
#> A muscadet object 
#>  2 omics: ATAC, RNA 
#>  types: ATAC, RNA 
#>  labels: scATAC-seq, scRNA-seq 
#>  cells: 71, 69 (common: 63, total: 77) 
#>  counts: 71 cells x 1200 features (peaks), 69 cells x 300 features (genes) 
#>  logratio: 71 cells x 213 features (windows of peaks), 69 cells x 212 features (genes) 
#>  variant positions: 681, 359 
#>  bulk data: WGS 
#>  clustering: partitions = 0.1, 0.3, 0.5 ; optimal partition = 0.5 
#>  CNA calling: 2 clusters ; 3 consensus segments including 0 CNA segments 
#>  genome: hg38 

# Overview of the muscomic objects within
show(slot(exdata_muscadet, "omics"))
#> $ATAC
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681
#> 
#> $RNA
#> A muscomic object 
#>  type: RNA 
#>  label: scRNA-seq 
#>  cells: 69 
#>  counts: 69 cells x 300 features (genes)
#>  logratio: 69 cells x 212 features (genes)
#>  variant positions: 359
#> 

# Overview of the first muscomic objects within
show(slot(exdata_muscadet, "omics")[[1]])
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681
```
