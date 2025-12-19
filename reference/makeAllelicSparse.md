# Construct Sparse Matrices for Allelic Data

This function converts a long-format allele counts table into two sparse
matrices: one containing reference allele depths (RD) and one containing
alternate allele depths (AD). Each matrix has cells as rows and variant
positions as columns. Variant coordinates and metadata (id, CHROM, POS,
REF, ALT) is also returned in a structured data frame in the list of
outputs.

## Usage

``` r
makeAllelicSparse(allele_counts)
```

## Arguments

- allele_counts:

  A data.frame containing at least the columns: `cell`, `CHROM`, `POS`,
  `REF`, `ALT`, `RD`, `AD`. Each row represents a single variant in a
  specific cell.

## Value

A list with the following elements:

- mat.RD:

  A `Matrix::dgCMatrix` sparse matrix of reference depths.

- mat.AD:

  A `Matrix::dgCMatrix` sparse matrix of alternate depths.

- coord.vars:

  A data.frame of variant metadata with columns: `id`, `CHROM`, `POS`,
  `REF`, `ALT`.

## Examples

``` r
allele_counts <- data.frame(
  cell = c("cell1","cell1","cell2"),
  CHROM = c("1","1","1"),
  POS = c(101,102,101),
  REF = c("A","G","A"),
  ALT = c("C","T","C"),
  RD = c(10, 0, 4),
  AD = c(3, 5, 0)
)
out <- makeAllelicSparse(allele_counts)
str(out)
#> List of 3
#>  $ mat.RD    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:2] 0 1
#>   .. ..@ p       : int [1:3] 0 2 2
#>   .. ..@ Dim     : int [1:2] 2 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:2] "cell1" "cell2"
#>   .. .. ..$ : chr [1:2] "1_101_A_C" "1_102_G_T"
#>   .. ..@ x       : num [1:2] 10 4
#>   .. ..@ factors : list()
#>  $ mat.AD    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:2] 0 0
#>   .. ..@ p       : int [1:3] 0 1 2
#>   .. ..@ Dim     : int [1:2] 2 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:2] "cell1" "cell2"
#>   .. .. ..$ : chr [1:2] "1_101_A_C" "1_102_G_T"
#>   .. ..@ x       : num [1:2] 3 5
#>   .. ..@ factors : list()
#>  $ coord.vars:'data.frame':  2 obs. of  5 variables:
#>   ..$ id   : chr [1:2] "1_101_A_C" "1_102_G_T"
#>   ..$ CHROM: Ord.factor w/ 1 level "1": 1 1
#>   ..$ POS  : int [1:2] 101 102
#>   ..$ REF  : chr [1:2] "A" "G"
#>   ..$ ALT  : chr [1:2] "C" "T"
# access cells names
rownames(out$mat.RD)
#> [1] "cell1" "cell2"
# access variants ids
colnames(out$mat.RD) # or
#> [1] "1_101_A_C" "1_102_G_T"
out$coord.vars$id
#> [1] "1_101_A_C" "1_102_G_T"
```
