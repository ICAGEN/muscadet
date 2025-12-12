# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature coordinates ----------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Example data: Feature coordinates
#'
#' Data frames of features (peaks, genes...) coordinates on genome. Only
#' chromosomes 3, 4 and 8 are present in the `exdata` example dataset.
#'
#' @name exdata_features
#' @rdname exdata_features
#'
#' @format
#' A data frame with the following columns:
#' \describe{
#'   \item{`CHROM`}{Chromosome names, e.g. "3", "X".}
#'   \item{`start`}{Start positions.}
#'   \item{`end`}{End positions.}
#'   \item{`id`}{Unique identifiers, e.g. gene name "CDH1" or peak identifier
#'   CHROM_start_end "1_1600338_1600838".}
#' }
#'
"exdata_genes"

#' @name exdata_features
#' @rdname exdata_features
#' @format NULL
#'
"exdata_peaks"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Matrices of counts -----------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Example data: Matrices of raw counts
#'
#' Matrices of raw read counts *features x cells* in
#' [`dgCMatrix`][Matrix::dgCMatrix-class] format.
#'
#' @name exdata_mat_counts
#' @rdname exdata_mat_counts
#'
#' @format
#' A [`dgCMatrix`][Matrix::dgCMatrix-class] of numeric values with the following
#' dimensions:
#' \describe{
#'   \item{`rows`}{Features (peaks, genes, ...).}
#'   \item{`columns`}{Cell barcodes.}
#' }
#'
"exdata_mat_counts_atac_tumor"

#' @name exdata_mat_counts
#' @rdname exdata_mat_counts
#' @format NULL
#'
"exdata_mat_counts_atac_ref"

#' @name exdata_mat_counts
#' @rdname exdata_mat_counts
#' @format NULL
#'
"exdata_mat_counts_rna_tumor"

#' @name exdata_mat_counts
#' @rdname exdata_mat_counts
#' @format NULL
#'
"exdata_mat_counts_rna_ref"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Allele counts ----------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Example data: Allele counts at variant positions
#'
#' Data frames of allele counts at variant positions per cell. Variant positions
#' can be either common single nucleotide polymorphisms (SNPs) positions or
#' individual-specific heterozygous positions retrieved using bulk sequencing.
#' Only chromosomes 3, 4 and 8 are present in the `exdata` example dataset.
#'
#' @name exdata_allele_counts
#' @rdname exdata_allele_counts
#'
#' @seealso
#' - [process_allele()]
#'
#' @format
#' A data frame with columns based on the [Variant Call Format
#' (VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) columns. It
#' contains the following columns:
#' \describe{
#'   \item{`cell`}{Barcodes of cells.}
#'   \item{`id`}{Variant unique identifier defined as
#'   CHROM_POS_REF_ALT, e.g. "3_3126620_T_G".}
#'   \item{`CHROM`}{Chromosome names, e.g. "3", "X".}
#'   \item{`POS`}{Position of the variant (1-base positions).}
#'   \item{`REF`}{Reference allele base, "A" "C" "G" or "T".}
#'   \item{`ALT`}{Alternative allele base, "A" "C" "G" or "T".}
#'   \item{`RD`}{Reference allele depth/count.}
#'   \item{`AD`}{Alternative allele depth/count.}
#'   \item{`DP`}{Total depth/count.}
#'   \item{`GT`}{Genotype: "0/1" or "1/0" if unphased; "0|1" or "1|0" if phased.}
#' }
#'
"exdata_allele_counts_atac_tumor"

#' @name exdata_allele_counts
#' @rdname exdata_allele_counts
#' @format NULL
#'
"exdata_allele_counts_atac_ref"

#' @name exdata_allele_counts
#' @rdname exdata_allele_counts
#' @format NULL
#'
"exdata_allele_counts_rna_tumor"

#' @name exdata_allele_counts
#' @rdname exdata_allele_counts
#' @format NULL
#'
"exdata_allele_counts_rna_ref"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Log ratio from bulk data -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Example data: Log R ratio from bulk sequencing data
#'
#' Data frame containing log R ratio values per genomic segments from bulk
#' sequencing data. Only chromosomes 3, 4 and 8 are present in the `exdata`
#' example dataset.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{`CHROM`}{Chromosome names, e.g. "3", "X".}
#'   \item{`start`}{Start position of the segment.}
#'   \item{`end`}{End position of the segment.}
#'   \item{`lrr`}{Log R ratio of the segment ("cnlr.median" column from
#'   [facets::fitcncf()] `$cncf` data frame).}
#' }
#'
#' @note Data obtained from whole genome sequencing (WGS) after using
#'   [facets::fitcncf()] from [`facets`][facets::facets-package] "Cellular
#'   Fraction and Copy Numbers from Tumor Sequencing" version `0.6.2`: `$cncf`
#'   data frame columns `chrom`, `start`, `end`, and `cnlr.median`.
#'
#' @references
#' \describe{
#'   \item{[`facets`][facets::facets-package] package}{Shen R, Seshan VE. FACETS:
#'   allele-specific copy number and clonal heterogeneity analysis tool for
#'   high-throughput DNA sequencing. Nucleic Acids Res. 2016 Sep 19;44(16):e131.
#'   doi: [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
#'   27270079; PMCID: PMC5027494.}
#' }
#'
"exdata_bulk_lrr"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genome chromosome sizes -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Genome chromosome sizes (internal data)
#'
#' @description [`GRanges`][GenomicRanges::GRanges()] objects containing
#'   chromosomes sizes for hg38, hg19 and mm10 genome assemblies.
#'
#' @name genome_chrom
#' @aliases hg38_chrom hg19_chrom mm10_chrom
#' @rdname genome_chrom
#'
#' @keywords internal
#'
#' @format A [`GRanges`][GenomicRanges::GRanges()] object containing:
#' \describe{
#'   \item{`seqnames`}{Chromosome name: 1 to 22, X and Y for human ; and 1 to
#'   19, X and Y for mouse (`Rle`).}
#'   \item{`ranges`}{Ranges of chromosomes (`IRanges`).}
#'   \item{`strand`}{Strand information * (`Rle`).}
#' }
#'
#' @note Data obtained from assemblies provided by the `BSgenome` package.
#' - `BSgenome.Hsapiens.UCSC.hg38` version 1.4.5 - `GRCh38.p14`
#' - `BSgenome.Hsapiens.UCSC.hg19` version 1.4.3 - `GRCh37.p13`
#' - `BSgenome.Mmusculus.UCSC.mm10` version 1.4.3 - `GRCm38.p6`
#'
#'
#' @references
#' \describe{
#'   \item{BSgenome package}{Pagès H (2024). BSgenome: Software infrastructure
#'   for efficient representation of full genomes and their SNPs.
#'   [https://bioconductor.org/packages/BSgenome](https://bioconductor.org/packages/BSgenome).}
#' }
#'
#' @examples
#' muscadet:::hg38_chrom
#' muscadet:::hg19_chrom
#' muscadet:::mm10_chrom
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# muscadet objects -------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Example data: muscadet objects
#'
#' [`muscadet`][muscadet-class] objects, containing two
#'   single-cell omic datasets: scATAC-seq and scRNA-seq.
#' - `exdata_muscadet` with tumor cells data: sample cells
#' - `exdata_muscadet_ref` with normal cells data: reference cells
#'
#' @name exdata_muscadet
#' @rdname exdata_muscadet
#'
#' @format [`muscadet`][muscadet-class] objects with the following slots:
#' \describe{
#'   \item{`omics`}{List of [`muscomic`][muscomic-class] objects, one per
#'   single-cell omic (`list`).}
#'   \item{`bulk.data`}{List of data from paired bulk sequencing (`list`).}
#'   \item{`clustering`}{List of data associated with the clustering of cells
#'   based on coverage log R ratio values (`list`).}
#'   \item{`cnacalling`}{List of data associated with the calling of copy number
#'   alterations (CNAs) (`list`).}
#'   \item{`genome`}{Reference genome name among: "hg38", "hg19" and "mm10" (`character`).}
#' }
#'
"exdata_muscadet"

#' @name exdata_muscadet_ref
#' @rdname exdata_muscadet
#' @format NULL
#'
"exdata_muscadet_ref"


