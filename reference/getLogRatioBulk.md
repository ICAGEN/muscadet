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
  ordered as followed: chromosome (`integer`), start position
  (`integer`), end position (`integer`), and Log R ratio value
  (`numeric`).

## Value

A data frame of single-cell omic features, with columns: `CHROM`,
`start`, `end`, and `bulk.lrr`, where `bulk.lrr` corresponds to the log
R ratio values retrieved from bulk data matching the feature
coordinates.

## See also

example data:
[`bulk_lrr`](https://icagen.github.io/muscadet/reference/bulk_lrr.md)
functions:
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md),
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)

## Examples

``` r
# Create a muscomic object
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor,
  allele_counts = allele_counts_atac_tumor,
  features = peaks
)
# or use a muscomic object inside a muscadet object
atac <- slot(muscadet_obj, "omics")[["ATAC"]]

# Load bulk Log R ratio data frame
data(bulk_lrr)
head(bulk_lrr)
#>   CHROM     start       end         lrr
#> 1     1     16100  78558400 -0.17422897
#> 2     1  78558700  80533900 -0.85006857
#> 3     1  80534279  80603200 -0.17219191
#> 4     1  80603500  99964900 -0.86685361
#> 5     1  99965200 122514900 -0.15990417
#> 6     1 122520400 124588600  0.05119736
# or use the one inside a muscadet object
head(slot(muscadet_obj, "bulk.data")[["log.ratio"]])
#>   CHROM     start       end         lrr
#> 1     1     16100  78558400 -0.17422897
#> 2     1  78558700  80533900 -0.85006857
#> 3     1  80534279  80603200 -0.17219191
#> 4     1  80603500  99964900 -0.86685361
#> 5     1  99965200 122514900 -0.15990417
#> 6     1 122520400 124588600  0.05119736

features_bulk_lrr <- getLogRatioBulk(
  x = atac,
  bulk.lrr = bulk_lrr # or slot(muscadet_obj, "bulk.data")[["log.ratio"]]
)
head(features_bulk_lrr)
#>   CHROM     start       end   bulk.lrr
#> 1     1 102000001 112000000 -0.1599042
#> 2     1 104000001 114000000 -0.1599042
#> 3     1 106000001 116000000 -0.1599042
#> 4     1 108000001 118000000 -0.1599042
#> 5     1 110000001 120000000 -0.1599042
#> 6     1 228000001 238000000 -0.1641702
```
