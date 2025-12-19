# Preparation of input data

Documentation on example data inputs:

- Count matrices:
  [`?exdata_mat_counts`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.html)
- Allele counts data frames:
  [`?exdata_allele_counts`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.html)
- Features coordinate data frames:
  [`?exdata_features`](https://icagen.github.io/muscadet/reference/exdata_features.html)
- Matched bulk coverage (log R ratio):
  [`?exdata_bulk_lrr`](https://icagen.github.io/muscadet/reference/exdata_bulk_lrr.html)

> **Important note on example data**
>
> The example dataset included in `muscadet` is a toy dataset designed
> for demonstration purposes only. It is deliberately minimal and
> contains a reduced genomic representation (three chromosomes). Its
> sole purpose is to illustrate how to run the main functions and
> explore the package features. Because of this strong simplification,
> the results obtained from this dataset should not be interpreted as
> biologically meaningful or methodologically representative.

## 1 Count matrices

The `mat_counts` input is required for muscadet analysis and must be
provided with
[`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md).

This input should be a matrix (sparse `dgCMatrix` or not) of raw read
counts with features as rows and cells as columns
([Table 1](#tbl-mat-counts)).

See
[`?exdata_mat_counts`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.md)
for full documentation.

``` r
data("exdata_mat_counts_atac_tumor")
kable(exdata_mat_counts_atac_tumor[1:8, 30:34])
```

|                               | 3_14380209_14380709 | 3_14432174_14432674 | 3_14946687_14947187 | 3_15050333_15050833 | 3_15241296_15241796 |
|:------------------------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
| samplename_AACTTAGTCCTCATCA-1 |                   0 |                   0 |                   0 |                   0 |                   0 |
| samplename_AAGAATCAGGTGTTAC-1 |                   0 |                   2 |                   0 |                   0 |                   0 |
| samplename_AAGCCACGTTAGTACG-1 |                   0 |                   0 |                   0 |                   0 |                   0 |
| samplename_AATAACCGTAGTTGGC-1 |                   0 |                   0 |                   0 |                   0 |                   0 |
| samplename_AATCCTAAGGTCCTAG-1 |                   0 |                   1 |                   0 |                   0 |                   0 |
| samplename_AATTTCCTCGAAGTGA-1 |                   0 |                   0 |                   0 |                   0 |                   0 |
| samplename_ACATCATCAGGCTTCG-1 |                   0 |                   3 |                   0 |                   0 |                   0 |
| samplename_ACCAAACTCTAACCTT-1 |                   0 |                   0 |                   0 |                   0 |                   0 |

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
- `REF`: reference allele
- `ALT`: alternate allele
- `RD`: count for reference allele
- `AD`: count for alternate allele
- `DP`: total depth
- `GT`: genotype *(optional)*

See
[`?exdata_allele_counts`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
for full documentation.

``` r
data("exdata_allele_counts_atac_tumor")
rownames(exdata_allele_counts_atac_tumor) <- NULL
kable(head(exdata_allele_counts_atac_tumor))
```

| cell                          | id            | CHROM |     POS | REF | ALT |  RD |  AD |  DP | GT   |
|:------------------------------|:--------------|------:|--------:|:----|:----|----:|----:|----:|:-----|
| samplename_CCGGTTAAGGAGCAAC-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   1 |   0 |   1 | 0\|1 |
| samplename_CTTGCTCAGTTAGCCG-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   0 |   1 |   1 | 0\|1 |
| samplename_CTTTAGTTCTAGCTAA-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   1 |   0 |   1 | 0\|1 |
| samplename_GCGCCTTGTTCCGGGA-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   0 |   1 |   1 | 0\|1 |
| samplename_GGTCTTGAGCAAGGGT-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   0 |   1 |   1 | 0\|1 |
| samplename_GTCAAACTCTAGCGTG-1 | 3_3126620_T_G |     3 | 3126620 | T   | G   |   1 |   0 |   1 | 0\|1 |

Table 2: Allele counts table example

### 2.1 Variant positions

Single variant positions can be derived from:

- Matched bulk sequencing to identify individual-specific heterozygous
  SNPs, or
- Reference panels of common SNPs (e.g. gnomAD, 1000G).

#### 2.1.1 Individual-specific heterozygous positions from bulk data

Single heterozygous positions can be retrieved by running
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

- [gnomAD](https://gnomad.broadinstitute.org/) database: [data
  here](https://gnomad.broadinstitute.org/downloads)
- [1000G](https://www.internationalgenome.org/) database: [data
  here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)

### 2.2 Count reads from single cells

[SCReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/)[²](#fn2)
is used to get read counts per allele from single cell BAM files.

> *From [SCReadCounts
> documentation](https://horvathlab.github.io/NGS/SCReadCounts/)*: The
> program `scReadCounts` manages the sequential execution of programs
> `readCounts` and `readCountsMatrix`, collects the necessary arguments
> for successful execution, and avoids unnecessary execution of the
> expensive `readCounts` tool if possible. `readCounts` requires three
> input files: a pooled single cell alignment, a list of genomic
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
> assessment of read counts.

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
[`?exdata_features`](https://icagen.github.io/muscadet/reference/exdata_features.md)
for full documentation.

``` r
data("exdata_peaks")
kable(head(exdata_peaks))
```

| CHROM |   start |     end | id                |
|:------|--------:|--------:|:------------------|
| 3     | 1535287 | 1535787 | 3_1535287_1535787 |
| 3     | 1594942 | 1595442 | 3_1594942_1595442 |
| 3     | 1632524 | 1633024 | 3_1632524_1633024 |
| 3     | 1964186 | 1964686 | 3_1964186_1964686 |
| 3     | 2793516 | 2794016 | 3_2793516_2794016 |
| 3     | 2970389 | 2970889 | 3_2970389_2970889 |

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

genes <- rownames(exdata_mat_counts_RNA)
genes_coord <- genes(EnsDb.Hsapiens.v86, filter = GeneNameFilter(genes))
genes_coord <- as.data.frame(genes_coord)
genes_coord <- genes_coord[genes_coord$seqnames %in% c(1:22, "X", "Y"), ]
genes_coord <- genes_coord[, c("seqnames", "start", "end", "gene_name")]
colnames(genes_coord) <- c("CHROM", "start", "end", "id")
rownames(genes_coord) <- NULL
```

> **Important note: Match cells and features**
>
> Cell names of matrices (rows) must match between assays. Features
> names of matrices (columns) must match with the ids in the features
> coordinates data frames.
>
> ``` r
> # Cell names format must match between different assays
> intersect(rownames(mat_counts_ATAC), rownames(mat_counts_RNA))
>
> # Make sure columns names of mat_counts match features ids
> identical(colnames(mat_counts_ATAC), peaks_coord[, "id"])
> colnames(mat_counts_ATAC) <- peaks_coord[, "id"]
>
> table(colnames(mat_counts_RNA) %in% genes_coord[, "id"])
> ```

## 4 Bulk log-ratios

A table of log ratios computed from matched bulk sequencing (e.g. Whole
Genome Sequencing) can optionally be provided to
[`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)
and will be displayed at the bottom of log ratio heatmaps for validation
purposes ([Table 6](#tbl-bulk-lrr)).

Log ratios for bulk sequencing can be obtained through
[FACETS](https://github.com/mskcc/facets)[³](#fn3) analysis.

See
[`?exdata_bulk_lrr`](https://icagen.github.io/muscadet/reference/exdata_bulk_lrr.md)
for full documentation.

``` r
data("exdata_bulk_lrr")
kable(head(exdata_bulk_lrr))
```

| CHROM |     start |       end |        lrr |
|------:|----------:|----------:|-----------:|
|     3 |     11900 |  64709600 |  0.1229632 |
|     3 |  64710100 |  65331100 | -0.1631371 |
|     3 |  65331200 |  93784800 |  0.1271516 |
|     3 |  93784900 | 165334700 |  0.1507467 |
|     3 | 165335000 | 165426322 | -0.2682795 |
|     3 | 165426671 | 197495900 |  0.1484524 |

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
