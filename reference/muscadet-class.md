# The muscadet class

A `muscadet` object encapsulates data from different single-cell omics
as
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
objects as well as downstream analysis results after clustering and CNA
calling.

## Slots

- `omics`:

  List of
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  objects, one per single-cell omic (`list`).

- `bulk.data`:

  List of objects containing data from paired bulk sequencing (`list`).

- `clustering`:

  List of objects containing data associated with the clustering of
  cells based on coverage log R ratio values (`list`).

- `cnacalling`:

  List of objects containing data associated with the calling of copy
  number alterations (CNAs) (`list`).

- `genome`:

  Reference genome name among: "hg38", "hg19" and "mm10" (`character`).

## See also

Functions related to `muscadet` objects:

- [`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)
  (create `muscadet` objects)

- [muscadet-access](https://icagen.github.io/muscadet/reference/musc-access.md)
  (simplified access and assignment methods)

- [.DollarNames.muscadet](https://icagen.github.io/muscadet/reference/musc-auto.md)
  (autocompletion for `$`)

- [`show,muscadet-method`](https://icagen.github.io/muscadet/reference/show-muscadet-method.md)
  (`show` method)

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")

exdata_muscadet
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

str(exdata_muscadet, max.level = 2)
#> Formal class 'muscadet' [package "muscadet"] with 5 slots
#>   ..@ omics     :List of 2
#>   ..@ bulk.data :List of 2
#>   ..@ clustering:List of 9
#>   ..@ cnacalling:List of 20
#>   ..@ genome    : chr "hg38"

```
