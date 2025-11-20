# Get consensus segments across clusters

This function processes segments data, in which clusters have different
segment breakpoints, to identify consensus segments across all clusters.
It groups breakpoints within a specified genomic distance
(`dist.breakpoints`) and calculates representative breakpoints for each
group based on cluster size and median coordinates.

## Usage

``` r
getSegConsensus(x, ncells, dist.breakpoints = 1e+06)
```

## Arguments

- x:

  A data frame containing cluster segments information with the
  following required columns:

  - `chrom`: Chromosome name (`factor` or `character`).

  - `start`: Start position of the segment (`numeric`).

  - `end`: End position of the segment (`numeric`).

  - `cluster`: Cluster identifier (`numeric`).

- ncells:

  A named vector specifying the number of cells per cluster (`numeric`
  vector). The names must match the cluster identifiers in the `cluster`
  column of `x`.

- dist.breakpoints:

  A numeric value specifying the minimum genomic distance between
  adjacent breakpoints to be grouped into the same consensus segment
  (`numeric` value). Default: `1e6`.

## Value

A data frame containing consensus segments with the following columns:

- `chrom`: Chromosome name.

- `start`: Start position of the consensus segment.

- `end`: End position of the consensus segment.

## See also

[`annotateSegments()`](https://icagen.github.io/muscadet/reference/annotateSegments.md),
[`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)

## Examples

``` r
# Example data frame
segs <- data.frame(
  chrom = c("chr1", "chr1", "chr2", "chr2"),
  start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
  end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
  cluster = c("1", "2", "1", "2")
)

# Generate consensus segments
consensus_segs <- getSegConsensus(segs,
                                  ncells = c("1" = 50, "2" = 30),
                                  dist.breakpoints = 1e6)
print(consensus_segs)
#>   chrom   start     end
#> 1  chr1 1200000 2500000
#> 2  chr2 3100000 5500000
```
