# Annotate consensus segments across cluster with CNA states

This function annotates consensus segments across clusters using
segment-level data from each cluster. It assigns states (i.e., gain,
loss, cnloh, neutral), computes proportions of cells per cluster, and
defines whether alterations are clonal.

## Usage

``` r
annotateSegments(
  x,
  consensus_segs,
  ncells,
  ploidy = 2,
  minoverlap = 1e+06,
  clonal.thresh = 0.9
)
```

## Arguments

- x:

  A data frame containing cluster segments information with the
  following required columns:

  - `chrom`: Chromosome name (`factor` or `character`).

  - `start`: Start position of the segment (`numeric`).

  - `end`: End position of the segment (`numeric`).

  - `cluster`: Cluster identifier (`numeric` or `integer`).

  - `tcn.em`: Total copy number (EM algorithm) (`numeric`).

  - `lcn.em`: Lower copy number (EM algorithm) (`numeric`).

- consensus_segs:

  A `data.frame` of unique consensus segments.

- ncells:

  A named vector specifying the number of cells per cluster (`numeric`
  vector). The names must match the cluster identifiers in the `cluster`
  column of `x`.

- ploidy:

  An `integer` specifying the expected ploidy (default is `2`). Used to
  determine whether a segment is classified as gain, loss, or neutral.

- minoverlap:

  A non-negative integer specifying the minimum overlap between a
  cluster specific segment and a consensus segment (see
  [`IRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html))
  (`integer`). Default: `1e6`.

- clonal.thresh:

  Minimum cell proportion to label a segment as clonal. Default: `1e6`.

## Value

A `data.frame` with the consensus segments annotated per cluster,
including CNA state, clonal status, and proportion of cells per cluster
and across clusters:

- `ncells`: Number of cells in cluster.

- `prop.cluster`: Proportion of cells in cluster relatively to the total
  number of cells.

- `gnl`: GNL (gain ; neutral ; loss) status as an integer value (1 ; 0 ;
  -1).

- `loh`: Loss of heterozygosity status (logical).

- `state`: State of segment (gain ; loss ; neu ; cnloh).

- `cna`: Segment CNA status (logical).

- `cna_state`: State of CNA segment (gain ; loss ; cnloh).

- `prop.tot`: Proportion of cells for the state across clusters.

- `state_clonal`: Clonal state of segment (gain ; loss ; neu ; cnloh).

- `cna_clonal`: Segment clonal CNA status (logical).

- `cna_clonal_state`: Clonal state of CNA segment (gain ; loss ; cnloh).

## See also

[`getSegConsensus()`](https://icagen.github.io/muscadet/reference/getSegConsensus.md),
[`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md).

## Examples

``` r
#' # Example data frame
segs <- data.frame(
    chrom = c("chr1", "chr1", "chr2", "chr2"),
    start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
    end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
    cluster = c("1", "2", "1", "2"),
    cf.em = c(1, 1, 1, 1),
    tcn.em = c(1, 2, 3, 3),
    lcn.em = c(0, 1, 1, 1)
)

consensus_segs <- getSegConsensus(
    x = segs,
    ncells = c("1" = 50, "2" = 30)
)

# Annotate consensus segments
table <- annotateSegments(
    x = segs,
    consensus_segs = consensus_segs,
    ncells = c("1" = 50, "2" = 30)
)

print(table)
#>   chrom   start     end id cluster cf.em tcn.em lcn.em ncells prop.cluster gnl
#> 1  chr1 1200000 2500000  1       1     1      1      0     50        0.625  -1
#> 2  chr2 3100000 5500000  2       1     1      3      1     50        0.625   1
#> 3  chr1 1200000 2500000  1       2     1      2      1     30        0.375   0
#> 4  chr2 3100000 5500000  2       2     1      3      1     30        0.375   1
#>     loh state   cna cna_state prop.tot state_clonal cna_clonal cna_clonal_state
#> 1  TRUE  loss  TRUE      loss    0.625         <NA>      FALSE             <NA>
#> 2 FALSE  gain  TRUE      gain    1.000         gain       TRUE             gain
#> 3 FALSE   neu FALSE      <NA>    0.375         <NA>      FALSE             <NA>
#> 4 FALSE  gain  TRUE      gain    1.000         gain       TRUE             gain
```
