# Example data: Allele counts at variation positions

Data frames of allele counts at variant positions per cell. Variant
positions can be either common single nucleotide polymorphisms (SNPs)
positions or individual-specific heterozygous positions retrieved by
bulk sequencing.

## Usage

``` r
allele_counts_atac_tumor

allele_counts_atac_ref

allele_counts_rna_tumor

allele_counts_rna_ref
```

## Format

A data frame with columns based on the [Variant Call Format
(VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) columns. It
contains the following columns:

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

  Genotype: "0/1" or "1/0" if unphased; "0\|1" or "1\|0" if phased.
  (`character`).

## See also

- [`process_allele()`](https://icagen.github.io/muscadet/reference/process_allele.md)
