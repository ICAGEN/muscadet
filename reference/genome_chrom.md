# Genome chromosome sizes (internal data)

[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
objects containing chromosomes sizes for hg38, hg19 and mm10 genome
assemblies.

## Format

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object containing:

- `seqnames`:

  Chromosome name: 1 to 22, X and Y for human ; and 1 to 19, X and Y for
  mouse (`Rle`).

- `ranges`:

  Ranges of chromosomes (`IRanges`).

- `strand`:

  Strand information \* (`Rle`).

## Note

Data obtained from assemblies provided by the `BSgenome` package.

- `BSgenome.Hsapiens.UCSC.hg38` version 1.4.5 - `GRCh38.p14`

- `BSgenome.Hsapiens.UCSC.hg19` version 1.4.3 - `GRCh37.p13`

- `BSgenome.Mmusculus.UCSC.mm10` version 1.4.3 - `GRCm38.p6`

## References

- BSgenome package:

  Pagès H (2024). BSgenome: Software infrastructure for efficient
  representation of full genomes and their SNPs.
  <https://bioconductor.org/packages/BSgenome>.

## Examples

``` r
muscadet:::hg38_chrom
#> GRanges object with 24 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]        1 1-248956422      *
#>    [2]        2 1-242193529      *
#>    [3]        3 1-198295559      *
#>    [4]        4 1-190214555      *
#>    [5]        5 1-181538259      *
#>    ...      ...         ...    ...
#>   [20]       20  1-64444167      *
#>   [21]       21  1-46709983      *
#>   [22]       22  1-50818468      *
#>   [23]        X 1-156040895      *
#>   [24]        Y  1-57227415      *
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
muscadet:::hg19_chrom
#> GRanges object with 24 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]        1 1-249250621      *
#>    [2]        2 1-243199373      *
#>    [3]        3 1-198022430      *
#>    [4]        4 1-191154276      *
#>    [5]        5 1-180915260      *
#>    ...      ...         ...    ...
#>   [20]       20  1-63025520      *
#>   [21]       21  1-48129895      *
#>   [22]       22  1-51304566      *
#>   [23]        X 1-155270560      *
#>   [24]        Y  1-59373566      *
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
muscadet:::mm10_chrom
#> GRanges object with 21 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]        1 1-195471971      *
#>    [2]        2 1-182113224      *
#>    [3]        3 1-160039680      *
#>    [4]        4 1-156508116      *
#>    [5]        5 1-151834684      *
#>    ...      ...         ...    ...
#>   [17]       17  1-94987271      *
#>   [18]       18  1-90702639      *
#>   [19]       19  1-61431566      *
#>   [20]        X 1-171031299      *
#>   [21]        Y  1-91744698      *
#>   -------
#>   seqinfo: 21 sequences from an unspecified genome; no seqlengths
```
