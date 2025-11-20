# Example data: Matrices of raw counts

Matrices of raw read counts *features x cells* in
[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
format.

## Usage

``` r
mat_counts_atac_tumor

mat_counts_atac_ref

mat_counts_rna_tumor

mat_counts_rna_ref
```

## Format

A [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html) of
numeric values with the following dimensions:

- `rows`:

  Features (peaks, genes).

- `columns`:

  Cell barcodes.
