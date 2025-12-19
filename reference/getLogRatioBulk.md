# Retrieve log R ratio from Bulk data on single-cell features (internal)

This internal function assigns the log R ratio (LRR) values from bulk
data segments to single-cell omic features. It matches the bulk data
segments to corresponding genomic features from the single-cell omic
data and returns a table of single-cell features with the corresponding
bulk LRR values.

## Usage

``` r
getLogRatioBulk(x, bulk.lrr)
```

## Arguments

- x:

  A `muscomic` object containing the single-cell omic data, which
  includes the feature coordinates (`muscomic`).

- bulk.lrr:

  A data frame containing log R ratio per genomic segments from bulk
  sequencing data (`data.frame`). One row per segment and 4 columns
  ordered as followed: chromosome (`character`), start position
  (`integer`), end position (`integer`), and Log R ratio value
  (`numeric`).

## Value

A data frame of single-cell omic features, with columns corresponding to
`CHROM`, `start`, `end`, and `bulk.lrr`, where `bulk.lrr` corresponds to
the log R ratio values retrieved from bulk data matching the feature
coordinates.

## See also

- example data:
  [`exdata_bulk_lrr`](https://icagen.github.io/muscadet/reference/exdata_bulk_lrr.md)

- functions:
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md),
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")

# muscomic object inside a muscadet object
exdata_muscadet$ATAC
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681

# Bulk Log R ratio data
head(exdata_muscadet$bulk.data$logratio)
#>   CHROM     start       end        lrr
#> 1     3     11900  64709600  0.1229632
#> 2     3  64710100  65331100 -0.1631371
#> 3     3  65331200  93784800  0.1271516
#> 4     3  93784900 165334700  0.1507467
#> 5     3 165335000 165426322 -0.2682795
#> 6     3 165426671 197495900  0.1484524

features_bulk_lrr <- getLogRatioBulk(
  x = exdata_muscadet$ATAC,
  bulk.lrr = exdata_muscadet$bulk.data$logratio
)
head(features_bulk_lrr)
#>   CHROM    start      end  bulk.lrr
#> 1     3        1 10000000 0.1229632
#> 2     3  2000001 12000000 0.1229632
#> 3     3  4000001 14000000 0.1229632
#> 4     3  6000001 16000000 0.1229632
#> 5     3  8000001 18000000 0.1229632
#> 6     3 10000001 20000000 0.1229632

```
