# Example data: Allele counts at variant positions

Data frames of allele counts at variant positions per cell. Variant
positions can be either common single nucleotide polymorphisms (SNPs)
positions or individual-specific heterozygous positions retrieved using
bulk sequencing. Only chromosomes 3, 4 and 8 are present in the `exdata`
example dataset.

## Usage

``` r
exdata_allele_counts_atac_tumor

exdata_allele_counts_atac_ref

exdata_allele_counts_rna_tumor

exdata_allele_counts_rna_ref
```

## Format

A data frame with columns based on the [Variant Call Format
(VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) columns. It
contains the following columns:

- `cell`:

  Barcodes of cells.

- `id`:

  Variant unique identifier defined as CHROM_POS_REF_ALT, e.g.
  "3_3126620_T_G".

- `CHROM`:

  Chromosome names, e.g. "3", "X".

- `POS`:

  Position of the variant (1-base positions).

- `REF`:

  Reference allele base, "A" "C" "G" or "T".

- `ALT`:

  Alternative allele base, "A" "C" "G" or "T".

- `RD`:

  Reference allele depth/count.

- `AD`:

  Alternative allele depth/count.

- `DP`:

  Total depth/count.

- `GT`:

  Genotype: "0/1" or "1/0" if unphased; "0\|1" or "1\|0" if phased.

## See also

- [`process_allele()`](https://icagen.github.io/muscadet/reference/process_allele.md)
