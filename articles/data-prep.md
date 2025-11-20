# Preparation of input data

## 1 Count matrices

The `mat_counts` input is required for muscadet analysis and must be
provided with
[`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md).

This input should be a matrix (sparse `dgCMatrix` or not) of raw read
counts with features as rows and cells as columns
([Table 1](#tbl-mat-counts)).

See
[`?mat_counts`](https://icagen.github.io/muscadet/reference/mat_counts.md)
for full documentation.

``` r
data("mat_counts_atac_tumor")
kable(mat_counts_atac_tumor[698:700,34:36])
```

|                        | samplename_CCGTTGCGTGTGTCCC-1 | samplename_CCTGGATCATTCCTCG-1 | samplename_CCTTCAATCCGCCAAA-1 |
|:-----------------------|------------------------------:|------------------------------:|------------------------------:|
| 12_118016401_118016901 |                             1 |                             0 |                             0 |
| 12_124564171_124564671 |                             0 |                             0 |                             0 |
| 12_132710447_132710947 |                             2 |                             0 |                             2 |

Table 1: Allele counts table example

- From Seurat/Signac analysis:

``` r
library(Seurat) # Seurat v5
library(Signac)

data("atac_small") # Example Seurat class object from Signac

# dgCMatrix matrix from Assay named "RNA"
mat_counts_RNA <- atac_small[["RNA"]]$counts 
# dgCMatrix matrix from ChromatinAssay named "peaks"
mat_counts_ATAC <- atac_small[["peaks"]]$counts 
```

- From ArchR/SummarizedExperiment analysis:

``` r
library(SummarizedExperiment)
# RangedSummarizedExperiment class object
se <- readRDS("<SummarizedExperiment_path>") 
# dgCMatrix matrix
mat_counts_ATAC <- assay(se) 
```

## 2 Allele counts tables

The `allele_counts` input is optional for the clustering step (allelic
information not used for clustering) but required for the CNA calling
step: it can be provided in the initial object creation step with
[`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)
or added to `muscadet` objects later on with
[`addAlleleCounts()`](https://icagen.github.io/muscadet/reference/addAlleleCounts.md).

This input should be in a data frame format containing allelic read
counts at specific variant positions across cells
([Table 2](#tbl-allele-counts)). The format follows a Variant Calling
Format (VCF)-like structure and should include:

- `cell`: unique cell barcode
- `id`: variant identifier (e.g. CHROM_POS_REF_ALT)
- `CHROM`: chromosome
- `POS`: position
- `REF` / `ALT`: reference and alternate alleles
- `RD` / `AD`: counts for reference and alternate alleles
- `DP`: total depth
- `GT`: genotype *(optional)*

See
[`?allele_counts`](https://icagen.github.io/muscadet/reference/allele_counts.md)
for full documentation.

``` r
data("allele_counts_atac_tumor")
kable(head(allele_counts_atac_tumor))
```

|       | cell                          | id             | CHROM |      POS | REF | ALT |  RD |  AD |  DP | GT   |
|:------|:------------------------------|:---------------|------:|---------:|:----|:----|----:|----:|----:|:-----|
| 4694  | samplename_GAAAGGCTCTAAGTGC-1 | 1_2323263_T_G  |     1 |  2323263 | T   | G   |   1 |   0 |   1 | 1\|0 |
| 13631 | samplename_ACTTACAAGTAAGGGC-1 | 1_4007524_T_C  |     1 |  4007524 | T   | C   |   1 |   0 |   1 | 1\|0 |
| 24816 | samplename_GTAGCGCTCCTTTACG-1 | 1_6599445_C_A  |     1 |  6599445 | C   | A   |   0 |   1 |   1 | 1\|0 |
| 31038 | samplename_GTACTTCGTTGTGACA-1 | 1_8201062_C_T  |     1 |  8201062 | C   | T   |   0 |   1 |   1 | 0\|1 |
| 34534 | samplename_AAGCTCCCAACCCTCC-1 | 1_9491207_C_A  |     1 |  9491207 | C   | A   |   0 |   1 |   1 | 0\|1 |
| 49929 | samplename_ATGGTTATCCAAATCA-1 | 1_15417791_C_A |     1 | 15417791 | C   | A   |   0 |   1 |   1 | 0\|1 |

Table 2: Allele counts table example

### 2.1 Variant positions

Variant positions can be derived from:

- Matched bulk sequencing to identify individual-specific heterozygous
  SNPs, or
- Reference panels of common SNPs (e.g. gnomAD, 1000G).

#### 2.1.1 Individual-specific heterozygous positions from bulk data

Positions can be retrieved by running
[FACETS](https://github.com/mskcc/facets)[¹](#fn1) on matched WGS/WES
normal samples.

Use `snp-pileup` to extract allele counts from bulk BAM files at known
sites (VCF of common SNPs from the [NCBI
database](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz)).
The output of `snp-pileup` can be filtered on depth and allele frequency
to keep only heterozygous positions well covered.

See [FACETS snp-pileup
documentation](https://github.com/mskcc/facets/blob/master/inst/extcode/README.txt).

#### 2.1.2 Panels of common SNPs

- [gnomAD](https://gnomad.broadinstitute.org/) database (4,099 SNPs):
  [data
  here](https://cloud.google.com/life-sciences/docs/resources/public-datasets/gnomad?hl=en)
- [1000G](https://www.internationalgenome.org/) database (2,548 SNPs):
  [data here](https://sourceforge.net/projects/cellsnp/files/SNPlist)

### 2.2 Count reads from single cells

[SCReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/)[²](#fn2)
is used to get read counts per allele from single cell BAM files.

> The program `scReadCounts` manages the sequential execution of
> programs `readCounts` and `readCountsMatrix`, collects the necessary
> arguments for successful execution, and avoids unnecessary execution
> of the expensive `readCounts` tool if possible. `readCounts` requires
> three input files: a pooled single cell alignment, a list of genomic
> positions of interest, and the barcodes file (barcodes.tsv) .
> Optionally, `readCounts` can be user-configured for *read filtering*
> and *cell-barcode handling*, including restriction to barcodes of
> interest (achieved through specifying barcodes file - barcodes.tsv).
> `readCounts` utilizes the barcode information from the pooled single
> cell alignments and **outputs the variant and reference read counts
> (nvar and nref, respectively), for each barcode (cell), restricted to
> those present in the barcodes file, in a tab separated text file**.
> This file is then used as an input for the second program -
> `readCountsMatrix` - which, upon providing an output prefix, generates
> two outputs: (1) a cell-position matrix with absolute nvar and nref
> counts, and (2) a cell-position matrix with the expressed VAFRNA.
> VAFRNA is estimated at a user-defined threshold of minimum required
> sequencing reads (minR); default minR = 5. `readCountsMatrix` is
> time-efficient and can be re-run multiple times at various minR
> thresholds. Together, these tools facilitate single-cell level
> assessment of read counts. *From
> <https://horvathlab.github.io/NGS/SCReadCounts/>*

As the desired output is the tab separated text file with variant and
reference read counts for each barcode, the command
[`readCounts`](https://horvathlab.github.io/NGS/ReadCounts/) is
preferred to avoid unnecessary large matrix outputs from
`readCountsMatrix`.

- [Command
  usage](https://horvathlab.github.io/NGS/ReadCounts/docs/Usage.html)
  for`readCounts`
- Details on [Read
  grouping](https://horvathlab.github.io/NGS/ReadCounts/docs/Grouping.html)
- Details on [Alignment
  filter](https://horvathlab.github.io/NGS/ReadCounts/docs/Filtering.html)
- Example of `readCounts` command, with `STARsolo_CB` read grouping and
  `MPileup` alignment filter:

``` bash
readCounts \
    -s <vcf_file.vcf> \
    -r <bam_file.bam> \
    -b None \
    -G STARsolo_CB \
    -f MPileup \
    -o <output_file.tsv> \
    -t <threads> \
    &> <log_file>
```

Then, the output of `readCounts` (table of variant and reference read
counts) can be formatted to fit `muscadet` input requirement with
[`process_allele()`](https://icagen.github.io/muscadet/reference/process_allele.md)
([Table 3](#tbl-process-allele1), [Table 4](#tbl-process-allele2)).

``` r
# Example data frame of readCounts results
readcounts <- data.frame(
  CHROM = c("1", "1", "2"),
  POS = c(10101, 20202, 30303),
  REF = c("A", "G", "T"),
  ALT = c("G", "A", "C"),
  ReadGroup = c("cell1", "cell1", "cell1"),
  SNVCountForward = c(5, 10, 3),
  SNVCountReverse = c(4, 6, 2),
  RefCountForward = c(20, 15, 10),
  RefCountReverse = c(18, 12, 8)
)
readcounts$SNVCount <- readcounts$SNVCountForward + readcounts$SNVCountReverse
readcounts$RefCount <- readcounts$RefCountForward +readcounts$RefCountReverse
readcounts$GoodReads <- readcounts$SNVCount + readcounts$RefCount
readcounts[, "%BadRead"] <- c(0, 0, 0)
readcounts$VAF <- round(readcounts$SNVCount / readcounts$GoodReads, 4)
```

``` r
kable(readcounts)
```

| CHROM |   POS | REF | ALT | ReadGroup | SNVCountForward | SNVCountReverse | RefCountForward | RefCountReverse | SNVCount | RefCount | GoodReads | %BadRead |    VAF |
|:------|------:|:----|:----|:----------|----------------:|----------------:|----------------:|----------------:|---------:|---------:|----------:|---------:|-------:|
| 1     | 10101 | A   | G   | cell1     |               5 |               4 |              20 |              18 |        9 |       38 |        47 |        0 | 0.1915 |
| 1     | 20202 | G   | A   | cell1     |              10 |               6 |              15 |              12 |       16 |       27 |        43 |        0 | 0.3721 |
| 2     | 30303 | T   | C   | cell1     |               3 |               2 |              10 |               8 |        5 |       18 |        23 |        0 | 0.2174 |

Table 3: Example of readCounts output

``` r
kable(process_allele(readcounts))
```

| cell  | id          | CHROM |   POS | REF | ALT |  RD |  AD |  DP |
|:------|:------------|:------|------:|:----|:----|----:|----:|----:|
| cell1 | 1_10101_A_G | 1     | 10101 | A   | G   |  38 |   9 |  47 |
| cell1 | 1_20202_G_A | 1     | 20202 | G   | A   |  27 |  16 |  43 |
| cell1 | 2_30303_T_C | 2     | 30303 | T   | C   |  18 |   5 |  23 |

Table 4: Formatted table of allele counts

## 3 Feature coordinates

A table of features coordinates with 4 columns (`CHROM`, `start`, `end`,
`id`) matching features of count matrices (section
[Section 1](#sec-mat-counts)) is required ([Table 5](#tbl-features)).

See
[`?features`](https://icagen.github.io/muscadet/reference/features.md)
for full documentation.

``` r
data("peaks")
kable(head(peaks))
```

| CHROM |    start |      end | id                  |
|:------|---------:|---------:|:--------------------|
| 1     |   941542 |   942042 | 1_941542_942042     |
| 1     |  3976396 |  3976896 | 1_3976396_3976896   |
| 1     |  6465944 |  6466444 | 1_6465944_6466444   |
| 1     |  6901466 |  6901966 | 1_6901466_6901966   |
| 1     |  7402483 |  7402983 | 1_7402483_7402983   |
| 1     | 10794813 | 10795313 | 1_10794813_10795313 |

Table 5: Features coordinates table example

### 3.1 Peaks

- From Seurat/Signac

``` r
library(Signac)
# Example Seurat class object from Signac
data("atac_small") 
# GRanges class object
peaks_coord <- atac_small[["peaks"]]$ranges 

peaks_coord <- as.data.frame(peaks_coord) 
peaks_coord$id <- paste(peaks_coord$seqnames,
                        peaks_coord$start,
                        peaks_coord$end,
                        sep = "-")
peaks_coord <- peaks_coord[, c("seqnames", "start", "end", "id")]
colnames(peaks_coord) <- c("CHROM", "start", "end", "id")
```

- From ArchR/SummarizedExperiment:

``` r
library(SummarizedExperiment)
# RangedSummarizedExperiment class object
se <- readRDS("<SummarizedExperiment_path>") 
# GRanges class object
peaks_coord <- rowRanges(se) 

peaks_coord <- as.data.frame(peaks_coord)
peaks_coord$id <- paste(peaks_coord$seqnames, peaks_coord$start, peaks_coord$end, sep = "-")
peaks_coord <- peaks_coord[, c("seqnames", "start", "end", "id")]
colnames(peaks_coord) <- c("CHROM", "start", "end", "id")
```

### 3.2 Genes

``` r
# Use your own annotation
genes_coord <- read.delim( "<genes_gtf_file>")
genes_coord <- genes_coord[, c("seqnames", "start", "end", "gene_name")]
colnames(genes_coord) <- c("CHROM", "start", "end", "id")

# Or get coordinates 
library(EnsDb.Hsapiens.v86) # for human
library(AnnotationDbi)

genes <- rownames(mat_counts_RNA)
genes_coord <- genes(EnsDb.Hsapiens.v86, filter = GeneNameFilter(genes))
genes_coord <- as.data.frame(genes_coord)
genes_coord <- genes_coord[genes_coord$seqnames %in% c(1:22, "X", "Y"), ]
genes_coord <- genes_coord[, c("seqnames", "start", "end", "gene_name")]
colnames(genes_coord) <- c("CHROM", "start", "end", "id")
rownames(genes_coord) <- NULL
```

### 3.3 Match cells and features

> **Important**
>
> Cell names (columns) and features names (rows) must match between
> assays and features coordinates table ids.

``` r
# Make sure row names of mat_counts match features ids
identical(rownames(mat_counts_ATAC), peaks_coord[, "id"])
rownames(mat_counts_ATAC) <- peaks_coord[, "id"]

table(rownames(mat_counts_RNA) %in% genes_coord[, "id"])

# Cell names format must match between different assays
intersect(colnames(mat_counts_ATAC), colnames(mat_counts_RNA))
```

## 4 Bulk log-ratios

A table of log ratios computed from matched bulk sequencing (e.g. Whole
Genome Sequencing) can optionally be provided to
[`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)
and will be displayed at the bottom of log ratio heatmaps for validation
purposes ([Table 6](#tbl-bulk-lrr)).

Log ratios for bulk sequencing can be obtained through
[FACETS](https://github.com/mskcc/facets)[³](#fn3) analysis.

See
[`?bulk_lrr`](https://icagen.github.io/muscadet/reference/bulk_lrr.md)
for full documentation.

``` r
data("bulk_lrr")
kable(head(bulk_lrr))
```

| CHROM |     start |       end |        lrr |
|------:|----------:|----------:|-----------:|
|     1 |     16100 |  78558400 | -0.1742290 |
|     1 |  78558700 |  80533900 | -0.8500686 |
|     1 |  80534279 |  80603200 | -0.1721919 |
|     1 |  80603500 |  99964900 | -0.8668536 |
|     1 |  99965200 | 122514900 | -0.1599042 |
|     1 | 122520400 | 124588600 |  0.0511974 |

Table 6: Bulk log ratios table example

------------------------------------------------------------------------

1.  Shen R, Seshan VE. **FACETS: allele-specific copy number and clonal
    heterogeneity analysis tool for high-throughput DNA sequencing.**
    *Nucleic Acids Res* (2016). <https://www.doi.org/10.1093/nar/gkw520>

2.  Prashant, N.M., Alomran, N., Chen, Y. et al. **SCReadCounts:
    estimation of cell-level SNVs expression from scRNA-seq data.** *BMC
    Genomics* (2021). <https://doi.org/10.1186/s12864-021-07974-8>

3.  Shen R, Seshan VE. **FACETS: allele-specific copy number and clonal
    heterogeneity analysis tool for high-throughput DNA sequencing.**
    *Nucleic Acids Res* (2016). <https://www.doi.org/10.1093/nar/gkw520>
