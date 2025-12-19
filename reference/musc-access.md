# Access and assignment methods for `muscadet` objects

Simplified access to omic datasets and slots inside
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects.

Assign new data in a `[muscadet::muscadet()]` or
`[muscadet::muscomic()]` object. For `muscadet` objects, the omic
datasets in the `omics` slot can be directly reassigned.

## Usage

``` r
# S3 method for class 'muscadet'
x[i, ...]

# S3 method for class 'muscomic'
x[i, ...]

# S3 method for class 'muscadet'
x$name

# S3 method for class 'muscomic'
x$name

# S3 method for class 'muscadet'
x[i] <- value

# S3 method for class 'muscomic'
x[i] <- value

# S3 method for class 'muscadet'
x$i <- value

# S3 method for class 'muscomic'
x$i <- value
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  or
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object.

- i:

  The name of the slot (or omic).

- ...:

  Other arguments (ignored).

- name:

  The name of the slot (or omic).

- value:

  The new value to assign.

## Value

The selected slot or the omic dataset (`[muscadet::muscomic()]` object)
for `muscadet` objects. The selected slot for `muscomic` objects.

The updated `[muscadet::muscadet()]` or `[muscadet::muscomic()]` object.

## Examples

``` r
# Load example muscadet object
# data("exdata_muscadet")

# Access to muscadet omics or slots
exdata_muscadet["ATAC"]
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681
exdata_muscadet["genome"]
#> [1] "hg38"

# Load example muscadet object
# data("exdata_muscadet")

# Access to muscomic slots
exdata_muscadet["ATAC"]["label.omic"]
#> [1] "scATAC-seq"

# Load example muscadet object
# data("exdata_muscadet")

# Access to muscadet omics or slots
exdata_muscadet$ATAC
#> A muscomic object 
#>  type: ATAC 
#>  label: scATAC-seq 
#>  cells: 71 
#>  counts: 71 cells x 1200 features (peaks)
#>  logratio: 71 cells x 213 features (windows of peaks)
#>  variant positions: 681
exdata_muscadet$genome
#> [1] "hg38"

# Load example muscadet object
# data("exdata_muscadet")

# Access to muscomic slots
exdata_muscadet$ATAC$label.omic
#> [1] "scATAC-seq"

# Load example muscadet object
# data("exdata_muscadet")

# Assign new data in muscadet object
exdata_muscadet["genome"] <- "hg38"

# Load example muscadet object
# data("exdata_muscadet")

# Assign new data in muscomic object
exdata_muscadet["ATAC"]["label.omic"] <- "scATAC-seq"

# Load example muscadet object
# data("exdata_muscadet")

# Assign new data in muscadet object
exdata_muscadet$genome <- "hg38"

# Load example muscadet object
# data("exdata_muscadet")

# Assign new data in muscomic object
exdata_muscadet$ATAC$label.omic <- "scATAC-seq"
```
