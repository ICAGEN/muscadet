# Save a data frame as a VCF File

This function saves a given data frame in [variant calling format
(VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) with necessary
headers.

## Usage

``` r
save_as_vcf(data, file, header = NULL)
```

## Arguments

- data:

  A data frame containing variant data (`data.frame`). The first column
  should represent the chromosome.

- file:

  A character string specifying the output file path with the `.vcf` or
  `.vcf.gz` extension (`character string`). If the file has the
  `.vcf.gz` extension it will be compressed.

- header:

  An optional custom header for the VCF format (`character` vector). By
  default, `NULL` sets a minimal header.

## Value

The function does not return a value but writes a VCF file to the
specified path.

## Details

The function ensures that the first column is named `#CHROM` and writes
standard VCF headers. It then appends the variant data in a
tab-separated format.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf_data <- data.frame(
  CHROM = c("chr1", "chr2"),
  POS = c(12345, 67890),
  ID = c(".", "."),
  REF = c("A", "T"),
  ALT = c("G", "C"),
  QUAL = c(".", "."),
  FILTER = c("PASS", "PASS"),
  INFO = c("AF=0.5", "AF=0.3"),
  FORMAT = c("GT", "GT"),
  SAMPLE = c("0/1", "1/1")
)
save_as_vcf(vcf_data, "output.vcf")
} # }
```
