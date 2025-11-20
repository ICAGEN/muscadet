# Methods for [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md) and [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md) objects

Methods to facilitate access to data within the
[`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
and
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects.

- `Cells()`: Get cell identifiers (addition of methods for `muscomic`
  and `muscadet` to
  [`SeuratObject::Cells()`](https://satijalab.github.io/seurat-object/reference/Cells.html)).

- `Features()`: Get feature identifiers (addition of methods for
  `muscomic` and `muscadet` to
  [`SeuratObject::Features()`](https://satijalab.github.io/seurat-object/reference/Cells.html)).

- `coordFeatures()`: Get coordinates of features data frames.

- `matCounts()`: Get raw count matrices.

- `matLogRatio()`: Get log R ratio matrices.

## Usage

``` r
coordFeatures(x)

matCounts(x)

matLogRatio(x)

# S3 method for class 'muscomic'
Cells(x, ...)

# S3 method for class 'muscadet'
Cells(x, ...)

# S3 method for class 'muscomic'
Features(x, ...)

# S3 method for class 'muscadet'
Features(x, ...)

# S4 method for class 'muscomic'
coordFeatures(x)

# S4 method for class 'muscadet'
coordFeatures(x)

# S4 method for class 'muscomic'
matCounts(x)

# S4 method for class 'muscadet'
matCounts(x)

# S4 method for class 'muscomic'
matLogRatio(x)

# S4 method for class 'muscadet'
matLogRatio(x)
```

## Arguments

- x:

  A
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  or
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object.

- ...:

  Other arguments (ignored).

## Value

`Cells`:

- if `x` is a
  [`muscomic()`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object: a vector of cell names.

- if `x` is a
  [`muscadet()`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object: a list of cell names vectors, one list element per omic.

`Features`:

- if `x` is a
  [`muscomic()`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object: a vector of feature names.

- if `x` is a
  [`muscadet()`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object: a list of feature names vectors, one list element per omic.

`coordFeatures`:

- if `x` is a
  [`muscomic()`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object: a data frame of feature coordinates.

- if `x` is a
  [`muscadet()`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object: a list of feature coordinates data frames, one list element
  per omic.

`matCounts`:

- if `x` is a
  [`muscomic()`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object: a
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  *features x cells*.

- if `x` is a
  [`muscadet()`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object: a list of
  [`dgCMatrices`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  *features x cells*, one list element per omic.

`matLogRatio`:

- if `x` is a
  [`muscomic()`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  object: a `matrix` *features x cells*.

- if `x` is a
  [`muscadet()`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object: a list of `matrices` *features x cells*, one list element per
  omic.

## See also

[muscomic](https://icagen.github.io/muscadet/reference/muscomic-class.md),
[muscadet](https://icagen.github.io/muscadet/reference/muscadet-class.md)

## Examples

``` r
library("SeuratObject")

# Load example muscadet object
# data("muscadet_obj")

Cells(muscadet_obj) # list of 2 cell names vectors for the 2 omics
#> $ATAC
#>   [1] "samplename_AACAAAGGTCATGCAA-1" "samplename_AACGGTAAGGATAAAC-1"
#>   [3] "samplename_AAGCATGAGGGACGCA-1" "samplename_AAGCTCCCAACCCTCC-1"
#>   [5] "samplename_AATGTCATCATGCAAC-1" "samplename_AATTAGCGTAAAGCGG-1"
#>   [7] "samplename_ACACGGACAGCTCATA-1" "samplename_ACCAATATCCTCATGC-1"
#>   [9] "samplename_ACCTAAATCCAGCACA-1" "samplename_ACTATGTCAGTAGGTG-1"
#>  [11] "samplename_ACTTACAAGTAAGGGC-1" "samplename_AGAAACTAGGTCAAAG-1"
#>  [13] "samplename_AGAAGGTGTAGTAAGA-1" "samplename_AGACCCGGTAACCACA-1"
#>  [15] "samplename_AGCATCCCACCTCACC-1" "samplename_AGCTTAATCGCAAACT-1"
#>  [17] "samplename_AGGATATAGTTAGTGC-1" "samplename_AGTTACATCCGGTTAG-1"
#>  [19] "samplename_ATATGTCCACCTCACC-1" "samplename_ATGACCAGTCACCAAA-1"
#>  [21] "samplename_ATGACTCAGGGCCACT-1" "samplename_ATGGTTATCCAAATCA-1"
#>  [23] "samplename_ATGTTGTCATTGTGCA-1" "samplename_ATGTTTGAGCTTAGTA-1"
#>  [25] "samplename_ATTTGCGCAATACTGT-1" "samplename_CAACTAGGTGTTTGAG-1"
#>  [27] "samplename_CAGATTCAGCATGTCG-1" "samplename_CAGCATGTCGGTTCCT-1"
#>  [29] "samplename_CATGCAAGTTGAGCCG-1" "samplename_CATTTGTTCAATGACC-1"
#>  [31] "samplename_CCAGTTTGTAATCACG-1" "samplename_CCCTCATAGGACAATG-1"
#>  [33] "samplename_CCGCTAAAGGCATGTT-1" "samplename_CCGTTGCGTGTGTCCC-1"
#>  [35] "samplename_CCTGGATCATTCCTCG-1" "samplename_CCTTCAATCCGCCAAA-1"
#>  [37] "samplename_CGCCTCATCGATTATG-1" "samplename_CGCTACTTCCCTCAGT-1"
#>  [39] "samplename_CGCTATGAGTGACCTG-1" "samplename_CGGATAAAGGCATTGT-1"
#>  [41] "samplename_CGGATTAGTCATGAGC-1" "samplename_CGGTGAGAGCTTGCTC-1"
#>  [43] "samplename_CGGTTTCTCAATTACG-1" "samplename_CGTAATGGTCACCTAT-1"
#>  [45] "samplename_CGTGCACAGGTAGCTT-1" "samplename_CGTGTGTCAAGGATTA-1"
#>  [47] "samplename_CTAGCTTGTTAATGCG-1" "samplename_CTAGTGAGTTAACGAT-1"
#>  [49] "samplename_CTATGTTTCTTGGATA-1" "samplename_CTCACACTCCCATAAA-1"
#>  [51] "samplename_CTCACTCAGCGGATAA-1" "samplename_CTCCGGACAATGCGCT-1"
#>  [53] "samplename_CTTCAAGCACTAGGTC-1" "samplename_GAAAGGCTCTAAGTGC-1"
#>  [55] "samplename_GAAGGAACACTCGCTC-1" "samplename_GAGCAAGGTGGTTCTT-1"
#>  [57] "samplename_GAGGTACAGCAAGGTA-1" "samplename_GAGGTTAAGGACCTTG-1"
#>  [59] "samplename_GATCAGGCAACTAGCC-1" "samplename_GATTGCGTCCTCACTA-1"
#>  [61] "samplename_GCAAGTGCAGGACCTT-1" "samplename_GCAATCTAGCTCCTTA-1"
#>  [63] "samplename_GCAGGAAGTGCATCGG-1" "samplename_GCAGGTTGTTATTGCC-1"
#>  [65] "samplename_GCCAGGTTCACAGGAA-1" "samplename_GCTGCAATCCGTAAAC-1"
#>  [67] "samplename_GCTGTAAGTTAAGGCC-1" "samplename_GCTTAGTAGTAATCCA-1"
#>  [69] "samplename_GGAGGTTAGGCTAATC-1" "samplename_GGCGGTAAGAACAAGT-1"
#>  [71] "samplename_GGGCTAACAGCCAGAA-1" "samplename_GGTACTAGTGTGAGGA-1"
#>  [73] "samplename_GGTCAAGCATGAGCAG-1" "samplename_GGTGAGCCAGCTCATA-1"
#>  [75] "samplename_GGTTAATGTATTCGTC-1" "samplename_GGTTAGCGTCAGGCAT-1"
#>  [77] "samplename_GGTTATGGTGTTTCAC-1" "samplename_GGTTGCTCAGGTTATT-1"
#>  [79] "samplename_GTACTTCGTTGTGACA-1" "samplename_GTAGCGCTCCTTTACG-1"
#>  [81] "samplename_GTCTAGCCAGGATGGC-1" "samplename_GTGTGCGGTTACCGGG-1"
#>  [83] "samplename_GTTAAGTGTAGGTGTC-1" "samplename_GTTAGACTCACACAGT-1"
#>  [85] "samplename_GTTGCCCGTTAGGATT-1" "samplename_GTTTAACCACCTGTAA-1"
#>  [87] "samplename_TAATGCATCGATCAGT-1" "samplename_TACATCAAGGCTAATC-1"
#>  [89] "samplename_TACGGTTAGTTTGTCT-1" "samplename_TATATCCTCCGCACAA-1"
#>  [91] "samplename_TATTAGCCAGGCTAAG-1" "samplename_TATTGACCAGCGCTTG-1"
#>  [93] "samplename_TCAGCCTTCACCTGCT-1" "samplename_TCAGCGATCTTTAAGG-1"
#>  [95] "samplename_TCCCTGGTCCCGCAAA-1" "samplename_TCCTGGTTCGAGGAAC-1"
#>  [97] "samplename_TCGATTAAGCCGCTAA-1" "samplename_TCTAACTTCTTTGACT-1"
#>  [99] "samplename_TCTTAGCGTCCAAGAC-1" "samplename_TCTTCAAGTAATGGCC-1"
#> [101] "samplename_TGAGAACCAAAGGTAC-1" "samplename_TGAGAACCAAGGTGCA-1"
#> [103] "samplename_TGCTCCGTCAAACCGT-1" "samplename_TGGTTCCTCTTGCTAT-1"
#> [105] "samplename_TGTGGCCAGGGACTAA-1" "samplename_TGTGGCCAGTAACCAC-1"
#> [107] "samplename_TGTTATGAGTAGCGCC-1" "samplename_TGTTGGCCAGGCTAAG-1"
#> [109] "samplename_TTAGCTGCATAAGCAA-1" "samplename_TTGCAGCCAGAGGGAG-1"
#> [111] "samplename_TTTAACCTCGGTTCCT-1" "samplename_TTTCTTGCAGAGGCTA-1"
#> 
#> $RNA
#>   [1] "samplename_AAAGGTTAGTCACCAG-1" "samplename_AACAAAGGTCATGCAA-1"
#>   [3] "samplename_AACGGTAAGGATAAAC-1" "samplename_AAGAACAGTTTGAGGC-1"
#>   [5] "samplename_AAGCATGAGGGACGCA-1" "samplename_AAGCTCCCAACCCTCC-1"
#>   [7] "samplename_AATAGAGGTTTGCGAA-1" "samplename_AATGTCATCATGCAAC-1"
#>   [9] "samplename_AATTAGCGTAAAGCGG-1" "samplename_AATTGCTCATCGCTTT-1"
#>  [11] "samplename_ACACGGACAGCTCATA-1" "samplename_ACCTAAATCCAGCACA-1"
#>  [13] "samplename_ACTATGTCAGTAGGTG-1" "samplename_ACTATGTCATTAAACC-1"
#>  [15] "samplename_ACTCACCTCATTACAG-1" "samplename_ACTTACAAGTAAGGGC-1"
#>  [17] "samplename_AGAAACTAGGTCAAAG-1" "samplename_AGACCCGGTAACCACA-1"
#>  [19] "samplename_AGCATCCCACCTCACC-1" "samplename_AGCTTAATCGCAAACT-1"
#>  [21] "samplename_AGGATATAGTTAGTGC-1" "samplename_ATACCGGTCTTAGGGT-1"
#>  [23] "samplename_ATATGTCCACCTCACC-1" "samplename_ATGACCAGTCACCAAA-1"
#>  [25] "samplename_ATGACTCAGGGCCACT-1" "samplename_ATGCAGGCAACTAGCC-1"
#>  [27] "samplename_ATGCGATTCAGCAAAG-1" "samplename_ATGGCTTAGTTGTCTT-1"
#>  [29] "samplename_ATGGTTATCCAAATCA-1" "samplename_ATGTTGTCATTGTGCA-1"
#>  [31] "samplename_ATTGCACAGGACCTGC-1" "samplename_CAACTAGGTGTTTGAG-1"
#>  [33] "samplename_CAAGAACCACCTACTT-1" "samplename_CAAGAACCATGAAGTA-1"
#>  [35] "samplename_CAGATTCAGCATGTCG-1" "samplename_CAGCATGTCGGTTCCT-1"
#>  [37] "samplename_CATGCAAGTTGAGCCG-1" "samplename_CATTTGTTCAATGACC-1"
#>  [39] "samplename_CCAGTTTGTAATCACG-1" "samplename_CCATATTTCAAGCGCC-1"
#>  [41] "samplename_CCGCTAAAGGCATGTT-1" "samplename_CCGTTGCGTGTGTCCC-1"
#>  [43] "samplename_CCTGGATCATTCCTCG-1" "samplename_CCTTCAATCCGCCAAA-1"
#>  [45] "samplename_CGAAGTAAGGGCTAAA-1" "samplename_CGCACACAGTTCCTGC-1"
#>  [47] "samplename_CGCCTCATCGATTATG-1" "samplename_CGCTACTTCCCTCAGT-1"
#>  [49] "samplename_CGGATAAAGGCATTGT-1" "samplename_CGGATTAGTCATGAGC-1"
#>  [51] "samplename_CGGTGAGAGCTTGCTC-1" "samplename_CGTAATGGTCACCTAT-1"
#>  [53] "samplename_CGTGCACAGGTAGCTT-1" "samplename_CGTGTGTCACTCGCTC-1"
#>  [55] "samplename_CTACTAAAGGGATTAG-1" "samplename_CTAGCTTGTTAATGCG-1"
#>  [57] "samplename_CTAGTGAGTTAACGAT-1" "samplename_CTATGTTTCTTGGATA-1"
#>  [59] "samplename_CTCACACTCCCATAAA-1" "samplename_CTCACTCAGGCTATGT-1"
#>  [61] "samplename_CTCATTAGTGGATTGC-1" "samplename_CTCCGGACAATGCGCT-1"
#>  [63] "samplename_GAAAGGCTCTAAGTGC-1" "samplename_GAAGGAACACTCGCTC-1"
#>  [65] "samplename_GAGATAAGTTTCCTCC-1" "samplename_GAGGTTAAGGACCTTG-1"
#>  [67] "samplename_GATCAGGCAACTAGCC-1" "samplename_GATTGCGTCCTCACTA-1"
#>  [69] "samplename_GCAAGTGCAGGACCTT-1" "samplename_GCAATCTAGCTCCTTA-1"
#>  [71] "samplename_GCCAGGAAGATGGACA-1" "samplename_GCCAGGTTCACAGGAA-1"
#>  [73] "samplename_GCTGTAAGTTAAGGCC-1" "samplename_GCTGTGATCATCCTCA-1"
#>  [75] "samplename_GCTTAGTAGTAATCCA-1" "samplename_GCTTGACCAGCCAGTT-1"
#>  [77] "samplename_GGAGGTTAGGCTAATC-1" "samplename_GGGCTAACAGCCAGAA-1"
#>  [79] "samplename_GGTACAAAGCAACAAG-1" "samplename_GGTACTAGTGTGAGGA-1"
#>  [81] "samplename_GGTCAAGCATGAGCAG-1" "samplename_GGTGAGCCAGCTCATA-1"
#>  [83] "samplename_GGTTAATGTATTCGTC-1" "samplename_GGTTATGGTGTTTCAC-1"
#>  [85] "samplename_GTACTTCGTTGTGACA-1" "samplename_GTAGCGCTCCTTTACG-1"
#>  [87] "samplename_GTCTAGCCAGGATGGC-1" "samplename_GTCTCACTCTTTAGGA-1"
#>  [89] "samplename_GTGTGCGGTTACCGGG-1" "samplename_GTTAAGTGTAGGTGTC-1"
#>  [91] "samplename_GTTAGACTCACACAGT-1" "samplename_GTTCACCTCATGCGTG-1"
#>  [93] "samplename_GTTGTGAGTTAAGCTG-1" "samplename_TAAAGCCTCTAAGGTC-1"
#>  [95] "samplename_TAATGCATCGATCAGT-1" "samplename_TACATCAAGGCTAATC-1"
#>  [97] "samplename_TACGGTTAGTTTGTCT-1" "samplename_TACTAAGTCCGCATGA-1"
#>  [99] "samplename_TAGTTGTCAAGCCAGA-1" "samplename_TATATCCTCCGCACAA-1"
#> [101] "samplename_TATTACCTCCTGAATA-1" "samplename_TCAGCCTTCACCTGCT-1"
#> [103] "samplename_TCAGCGATCTTTAAGG-1" "samplename_TCCCTGGTCCCGCAAA-1"
#> [105] "samplename_TCCTGGTTCGAGGAAC-1" "samplename_TCGATTAAGCCGCTAA-1"
#> [107] "samplename_TCTAACTTCTTTGACT-1" "samplename_TCTTAGCGTCCAAGAC-1"
#> [109] "samplename_TGAGAACCAAAGGTAC-1" "samplename_TGAGAACCAAGGTGCA-1"
#> [111] "samplename_TGCGCGAGTTGTTGGA-1" "samplename_TGGTTCCTCTTGCTAT-1"
#> [113] "samplename_TGTAACTCAAATTGCT-1" "samplename_TGTGGCCAGGGACTAA-1"
#> [115] "samplename_TGTGGCCAGTAACCAC-1" "samplename_TGTTGGCCAGGCTAAG-1"
#> [117] "samplename_TTAGCTGCATAAGCAA-1" "samplename_TTTAACCTCCTTGCGT-1"
#> [119] "samplename_TTTAACCTCGGTTCCT-1"
#> 
Cells(muscadet_obj)$ATAC # cell names vector from the omic ATAC
#>   [1] "samplename_AACAAAGGTCATGCAA-1" "samplename_AACGGTAAGGATAAAC-1"
#>   [3] "samplename_AAGCATGAGGGACGCA-1" "samplename_AAGCTCCCAACCCTCC-1"
#>   [5] "samplename_AATGTCATCATGCAAC-1" "samplename_AATTAGCGTAAAGCGG-1"
#>   [7] "samplename_ACACGGACAGCTCATA-1" "samplename_ACCAATATCCTCATGC-1"
#>   [9] "samplename_ACCTAAATCCAGCACA-1" "samplename_ACTATGTCAGTAGGTG-1"
#>  [11] "samplename_ACTTACAAGTAAGGGC-1" "samplename_AGAAACTAGGTCAAAG-1"
#>  [13] "samplename_AGAAGGTGTAGTAAGA-1" "samplename_AGACCCGGTAACCACA-1"
#>  [15] "samplename_AGCATCCCACCTCACC-1" "samplename_AGCTTAATCGCAAACT-1"
#>  [17] "samplename_AGGATATAGTTAGTGC-1" "samplename_AGTTACATCCGGTTAG-1"
#>  [19] "samplename_ATATGTCCACCTCACC-1" "samplename_ATGACCAGTCACCAAA-1"
#>  [21] "samplename_ATGACTCAGGGCCACT-1" "samplename_ATGGTTATCCAAATCA-1"
#>  [23] "samplename_ATGTTGTCATTGTGCA-1" "samplename_ATGTTTGAGCTTAGTA-1"
#>  [25] "samplename_ATTTGCGCAATACTGT-1" "samplename_CAACTAGGTGTTTGAG-1"
#>  [27] "samplename_CAGATTCAGCATGTCG-1" "samplename_CAGCATGTCGGTTCCT-1"
#>  [29] "samplename_CATGCAAGTTGAGCCG-1" "samplename_CATTTGTTCAATGACC-1"
#>  [31] "samplename_CCAGTTTGTAATCACG-1" "samplename_CCCTCATAGGACAATG-1"
#>  [33] "samplename_CCGCTAAAGGCATGTT-1" "samplename_CCGTTGCGTGTGTCCC-1"
#>  [35] "samplename_CCTGGATCATTCCTCG-1" "samplename_CCTTCAATCCGCCAAA-1"
#>  [37] "samplename_CGCCTCATCGATTATG-1" "samplename_CGCTACTTCCCTCAGT-1"
#>  [39] "samplename_CGCTATGAGTGACCTG-1" "samplename_CGGATAAAGGCATTGT-1"
#>  [41] "samplename_CGGATTAGTCATGAGC-1" "samplename_CGGTGAGAGCTTGCTC-1"
#>  [43] "samplename_CGGTTTCTCAATTACG-1" "samplename_CGTAATGGTCACCTAT-1"
#>  [45] "samplename_CGTGCACAGGTAGCTT-1" "samplename_CGTGTGTCAAGGATTA-1"
#>  [47] "samplename_CTAGCTTGTTAATGCG-1" "samplename_CTAGTGAGTTAACGAT-1"
#>  [49] "samplename_CTATGTTTCTTGGATA-1" "samplename_CTCACACTCCCATAAA-1"
#>  [51] "samplename_CTCACTCAGCGGATAA-1" "samplename_CTCCGGACAATGCGCT-1"
#>  [53] "samplename_CTTCAAGCACTAGGTC-1" "samplename_GAAAGGCTCTAAGTGC-1"
#>  [55] "samplename_GAAGGAACACTCGCTC-1" "samplename_GAGCAAGGTGGTTCTT-1"
#>  [57] "samplename_GAGGTACAGCAAGGTA-1" "samplename_GAGGTTAAGGACCTTG-1"
#>  [59] "samplename_GATCAGGCAACTAGCC-1" "samplename_GATTGCGTCCTCACTA-1"
#>  [61] "samplename_GCAAGTGCAGGACCTT-1" "samplename_GCAATCTAGCTCCTTA-1"
#>  [63] "samplename_GCAGGAAGTGCATCGG-1" "samplename_GCAGGTTGTTATTGCC-1"
#>  [65] "samplename_GCCAGGTTCACAGGAA-1" "samplename_GCTGCAATCCGTAAAC-1"
#>  [67] "samplename_GCTGTAAGTTAAGGCC-1" "samplename_GCTTAGTAGTAATCCA-1"
#>  [69] "samplename_GGAGGTTAGGCTAATC-1" "samplename_GGCGGTAAGAACAAGT-1"
#>  [71] "samplename_GGGCTAACAGCCAGAA-1" "samplename_GGTACTAGTGTGAGGA-1"
#>  [73] "samplename_GGTCAAGCATGAGCAG-1" "samplename_GGTGAGCCAGCTCATA-1"
#>  [75] "samplename_GGTTAATGTATTCGTC-1" "samplename_GGTTAGCGTCAGGCAT-1"
#>  [77] "samplename_GGTTATGGTGTTTCAC-1" "samplename_GGTTGCTCAGGTTATT-1"
#>  [79] "samplename_GTACTTCGTTGTGACA-1" "samplename_GTAGCGCTCCTTTACG-1"
#>  [81] "samplename_GTCTAGCCAGGATGGC-1" "samplename_GTGTGCGGTTACCGGG-1"
#>  [83] "samplename_GTTAAGTGTAGGTGTC-1" "samplename_GTTAGACTCACACAGT-1"
#>  [85] "samplename_GTTGCCCGTTAGGATT-1" "samplename_GTTTAACCACCTGTAA-1"
#>  [87] "samplename_TAATGCATCGATCAGT-1" "samplename_TACATCAAGGCTAATC-1"
#>  [89] "samplename_TACGGTTAGTTTGTCT-1" "samplename_TATATCCTCCGCACAA-1"
#>  [91] "samplename_TATTAGCCAGGCTAAG-1" "samplename_TATTGACCAGCGCTTG-1"
#>  [93] "samplename_TCAGCCTTCACCTGCT-1" "samplename_TCAGCGATCTTTAAGG-1"
#>  [95] "samplename_TCCCTGGTCCCGCAAA-1" "samplename_TCCTGGTTCGAGGAAC-1"
#>  [97] "samplename_TCGATTAAGCCGCTAA-1" "samplename_TCTAACTTCTTTGACT-1"
#>  [99] "samplename_TCTTAGCGTCCAAGAC-1" "samplename_TCTTCAAGTAATGGCC-1"
#> [101] "samplename_TGAGAACCAAAGGTAC-1" "samplename_TGAGAACCAAGGTGCA-1"
#> [103] "samplename_TGCTCCGTCAAACCGT-1" "samplename_TGGTTCCTCTTGCTAT-1"
#> [105] "samplename_TGTGGCCAGGGACTAA-1" "samplename_TGTGGCCAGTAACCAC-1"
#> [107] "samplename_TGTTATGAGTAGCGCC-1" "samplename_TGTTGGCCAGGCTAAG-1"
#> [109] "samplename_TTAGCTGCATAAGCAA-1" "samplename_TTGCAGCCAGAGGGAG-1"
#> [111] "samplename_TTTAACCTCGGTTCCT-1" "samplename_TTTCTTGCAGAGGCTA-1"
Cells(muscadet_obj$ATAC) # cell names vector from the ATAC muscomic object
#>   [1] "samplename_AACAAAGGTCATGCAA-1" "samplename_AACGGTAAGGATAAAC-1"
#>   [3] "samplename_AAGCATGAGGGACGCA-1" "samplename_AAGCTCCCAACCCTCC-1"
#>   [5] "samplename_AATGTCATCATGCAAC-1" "samplename_AATTAGCGTAAAGCGG-1"
#>   [7] "samplename_ACACGGACAGCTCATA-1" "samplename_ACCAATATCCTCATGC-1"
#>   [9] "samplename_ACCTAAATCCAGCACA-1" "samplename_ACTATGTCAGTAGGTG-1"
#>  [11] "samplename_ACTTACAAGTAAGGGC-1" "samplename_AGAAACTAGGTCAAAG-1"
#>  [13] "samplename_AGAAGGTGTAGTAAGA-1" "samplename_AGACCCGGTAACCACA-1"
#>  [15] "samplename_AGCATCCCACCTCACC-1" "samplename_AGCTTAATCGCAAACT-1"
#>  [17] "samplename_AGGATATAGTTAGTGC-1" "samplename_AGTTACATCCGGTTAG-1"
#>  [19] "samplename_ATATGTCCACCTCACC-1" "samplename_ATGACCAGTCACCAAA-1"
#>  [21] "samplename_ATGACTCAGGGCCACT-1" "samplename_ATGGTTATCCAAATCA-1"
#>  [23] "samplename_ATGTTGTCATTGTGCA-1" "samplename_ATGTTTGAGCTTAGTA-1"
#>  [25] "samplename_ATTTGCGCAATACTGT-1" "samplename_CAACTAGGTGTTTGAG-1"
#>  [27] "samplename_CAGATTCAGCATGTCG-1" "samplename_CAGCATGTCGGTTCCT-1"
#>  [29] "samplename_CATGCAAGTTGAGCCG-1" "samplename_CATTTGTTCAATGACC-1"
#>  [31] "samplename_CCAGTTTGTAATCACG-1" "samplename_CCCTCATAGGACAATG-1"
#>  [33] "samplename_CCGCTAAAGGCATGTT-1" "samplename_CCGTTGCGTGTGTCCC-1"
#>  [35] "samplename_CCTGGATCATTCCTCG-1" "samplename_CCTTCAATCCGCCAAA-1"
#>  [37] "samplename_CGCCTCATCGATTATG-1" "samplename_CGCTACTTCCCTCAGT-1"
#>  [39] "samplename_CGCTATGAGTGACCTG-1" "samplename_CGGATAAAGGCATTGT-1"
#>  [41] "samplename_CGGATTAGTCATGAGC-1" "samplename_CGGTGAGAGCTTGCTC-1"
#>  [43] "samplename_CGGTTTCTCAATTACG-1" "samplename_CGTAATGGTCACCTAT-1"
#>  [45] "samplename_CGTGCACAGGTAGCTT-1" "samplename_CGTGTGTCAAGGATTA-1"
#>  [47] "samplename_CTAGCTTGTTAATGCG-1" "samplename_CTAGTGAGTTAACGAT-1"
#>  [49] "samplename_CTATGTTTCTTGGATA-1" "samplename_CTCACACTCCCATAAA-1"
#>  [51] "samplename_CTCACTCAGCGGATAA-1" "samplename_CTCCGGACAATGCGCT-1"
#>  [53] "samplename_CTTCAAGCACTAGGTC-1" "samplename_GAAAGGCTCTAAGTGC-1"
#>  [55] "samplename_GAAGGAACACTCGCTC-1" "samplename_GAGCAAGGTGGTTCTT-1"
#>  [57] "samplename_GAGGTACAGCAAGGTA-1" "samplename_GAGGTTAAGGACCTTG-1"
#>  [59] "samplename_GATCAGGCAACTAGCC-1" "samplename_GATTGCGTCCTCACTA-1"
#>  [61] "samplename_GCAAGTGCAGGACCTT-1" "samplename_GCAATCTAGCTCCTTA-1"
#>  [63] "samplename_GCAGGAAGTGCATCGG-1" "samplename_GCAGGTTGTTATTGCC-1"
#>  [65] "samplename_GCCAGGTTCACAGGAA-1" "samplename_GCTGCAATCCGTAAAC-1"
#>  [67] "samplename_GCTGTAAGTTAAGGCC-1" "samplename_GCTTAGTAGTAATCCA-1"
#>  [69] "samplename_GGAGGTTAGGCTAATC-1" "samplename_GGCGGTAAGAACAAGT-1"
#>  [71] "samplename_GGGCTAACAGCCAGAA-1" "samplename_GGTACTAGTGTGAGGA-1"
#>  [73] "samplename_GGTCAAGCATGAGCAG-1" "samplename_GGTGAGCCAGCTCATA-1"
#>  [75] "samplename_GGTTAATGTATTCGTC-1" "samplename_GGTTAGCGTCAGGCAT-1"
#>  [77] "samplename_GGTTATGGTGTTTCAC-1" "samplename_GGTTGCTCAGGTTATT-1"
#>  [79] "samplename_GTACTTCGTTGTGACA-1" "samplename_GTAGCGCTCCTTTACG-1"
#>  [81] "samplename_GTCTAGCCAGGATGGC-1" "samplename_GTGTGCGGTTACCGGG-1"
#>  [83] "samplename_GTTAAGTGTAGGTGTC-1" "samplename_GTTAGACTCACACAGT-1"
#>  [85] "samplename_GTTGCCCGTTAGGATT-1" "samplename_GTTTAACCACCTGTAA-1"
#>  [87] "samplename_TAATGCATCGATCAGT-1" "samplename_TACATCAAGGCTAATC-1"
#>  [89] "samplename_TACGGTTAGTTTGTCT-1" "samplename_TATATCCTCCGCACAA-1"
#>  [91] "samplename_TATTAGCCAGGCTAAG-1" "samplename_TATTGACCAGCGCTTG-1"
#>  [93] "samplename_TCAGCCTTCACCTGCT-1" "samplename_TCAGCGATCTTTAAGG-1"
#>  [95] "samplename_TCCCTGGTCCCGCAAA-1" "samplename_TCCTGGTTCGAGGAAC-1"
#>  [97] "samplename_TCGATTAAGCCGCTAA-1" "samplename_TCTAACTTCTTTGACT-1"
#>  [99] "samplename_TCTTAGCGTCCAAGAC-1" "samplename_TCTTCAAGTAATGGCC-1"
#> [101] "samplename_TGAGAACCAAAGGTAC-1" "samplename_TGAGAACCAAGGTGCA-1"
#> [103] "samplename_TGCTCCGTCAAACCGT-1" "samplename_TGGTTCCTCTTGCTAT-1"
#> [105] "samplename_TGTGGCCAGGGACTAA-1" "samplename_TGTGGCCAGTAACCAC-1"
#> [107] "samplename_TGTTATGAGTAGCGCC-1" "samplename_TGTTGGCCAGGCTAAG-1"
#> [109] "samplename_TTAGCTGCATAAGCAA-1" "samplename_TTGCAGCCAGAGGGAG-1"
#> [111] "samplename_TTTAACCTCGGTTCCT-1" "samplename_TTTCTTGCAGAGGCTA-1"
```
