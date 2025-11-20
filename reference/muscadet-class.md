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
# data("muscadet_obj")

muscadet_obj
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

str(muscadet_obj, max.level = 2)
#> Formal class 'muscadet' [package "muscadet"] with 5 slots
#>   ..@ omics     :List of 2
#>   ..@ bulk.data :List of 2
#>   ..@ clustering:List of 9
#>   ..@ cnacalling:List of 21
#>   ..@ genome    : chr "hg38"

```
