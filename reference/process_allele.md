# Process allele counts results from SCReadCounts format

This function reformats SCReadCounts output into a standardized VCF-like
data frame and optionally merges unphased/phased genotype information
from a VCF-format data frame.

## Usage

``` r
process_allele(data, vcf = NULL, samplename = NULL)
```

## Arguments

- data:

  A data frame of SCReadCounts output table (`data.frame`).

- vcf:

  An optional VCF-format data frame with unphased/phased genotype
  information to add it in output (`data.frame`).

- samplename:

  An optional sample name to prefix cell barcodes (`character`).

## Value

A cleaned data frame with VCF-like standardized columns:

- `cell`:

  Barcodes of cells (`character`).

- `id`:

  Variant unique identifier defined as CHROM_POS_REF_ALT, e.g.
  "1_920949_C_G" (`character`).

- `CHROM`:

  Chromosome in integer format, e.g. 15 (X and Y chromosomes are not
  included) (`integer`).

- `POS`:

  Position of the variant (1-base positions) (`integer`).

- `REF`:

  Reference allele base, "A" "C" "G" or "T" (`character`).

- `ALT`:

  Alternative allele base, "A" "C" "G" or "T" (`character`).

- `RD`:

  Reference allele depth/count (`integer`).

- `AD`:

  Alternative allele depth/count (`integer`).

- `DP`:

  Total depth/count (`integer`).

- `GT`:

  Genotype (if `vcf` provided): "0/1" or "1/0" if unphased; "0\|1" or
  "1\|0" if phased. (`character`).

## Note

- [SCReadCounts website](https://horvathlab.github.io/NGS/SCReadCounts/)

- [Variant Call Format (VCF)
  format](https://en.wikipedia.org/wiki/Variant_Call_Format)

## References

- [SCReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/):

  Liu, H., Bousounis, P., Movassagh, M., Edwards, N., and Horvath, A.
  SCReadCounts: estimation of cell-level SNVs expression from scRNA-seq
  data. BMC Genomics 22, 689 (2021). doi:
  [10.1186/s12864-021-07974-8](https://doi.org/10.1186/s12864-021-07974-8)

## Examples

``` r
# Generate minimal data
sc_data <- data.frame(
  ReadGroup = c("cell1", "cell2"),
  CHROM = c(1, 1),
  POS = c(1001, 1002),
  REF = c("A", "C"),
  ALT = c("G", "T"),
  SNVCount = c(3, 5),
  RefCount = c(7, 5),
  GoodReads = c(10, 10)
)

vcf_data <- data.frame(
  CHROM = c(1, 1),
  POS = c(1001, 1002),
  ID = c(".", "."),
  REF = c("A", "C"),
  ALT = c("G", "T"),
  QUAL = c(".", "."),
  FILTER = c(".", "."),
  INFO = c(".", "."),
  FORMAT = c("GT", "GT"),
  sample1 = c("0|1", "1|0"),
  stringsAsFactors = FALSE
)

data <- process_allele(sc_data, vcf = vcf_data, samplename = "sampleA")
```
