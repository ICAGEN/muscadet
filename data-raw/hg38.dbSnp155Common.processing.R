# Processing script for UCSC dbSNP 155 Common SNP data (hg38)
#
# Input:
# - hg38.dbSnp155Common.gz — downloaded from UCSC Table Browser
#
# Output:
# - hg38.dbSnp155Common.vcf.gz - all positions, VCF-like format
# - hg38.dbSnp155Common.filtered.vcf.gz - MAF-filtered positions, VCF-like format
#
# The hg38.dbSnp155Common.filtered.vcf.gz file is provided as GitHub release asset.

library(dplyr)
library(data.table)
library(stringi)

# --- Parameters ---------------------------------------------------------------

# Path to the raw downloaded file (not versioned, must be present locally)
input_file <- "data-raw/hg38.dbSnp155Common.gz"

# Output directory for processed files
output_dir <- "data-raw"

# MAF threshold for filtered output
MAF_THRESHOLD <- 0.2

# MAF projects used: 1000 Genomes (col 1), TOPMED (col 3), gnomAD (col 13)
MAF_COL_IDX <- c(1, 3, 13)


# --- Load raw data ------------------------------------------------------------

# UCSC column names:
# V1: chrom, V2: chromEnd (1-based POS), V3: name (rsID), V4: ref, V5: alts,
# V6: minorAlleleFreq, V7: majorAllele, V8: minorAllele
# Values in V6, V7, V8 are comma-separated across 31 population frequency
# projects in a fixed order.

commonSNP <- read.table(
    input_file,
    sep = "\t",
    comment.char = "#",
    stringsAsFactors = FALSE,
    header = FALSE
)

colnames(commonSNP) <- c("CHROM",
                         "POS",
                         "ID",
                         "REF",
                         "ALT",
                         "MAF",
                         "majorAllele",
                         "minorAllele")

# --- Filter chromosomes and sort ---------------------------------------------

chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
commonSNP <- commonSNP %>%
    filter(CHROM %in% chroms) %>%
    mutate(CHROM = factor(CHROM, levels = chroms)) %>%
    arrange(CHROM, POS)

setDT(commonSNP)


# --- Extract max MAF and corresponding minor allele --------------------------
# For each position, retain the highest MAF across the selected projects
# and the corresponding minor allele.

maf_mat <- stri_split_fixed(commonSNP$MAF, ",", simplify = TRUE)[, MAF_COL_IDX]
minor_mat <- stri_split_fixed(commonSNP$minorAllele, ",", simplify = TRUE)[, MAF_COL_IDX]

# Convert to numeric (-inf strings become -Inf, handled by is.finite)
storage.mode(maf_mat) <- "numeric"

# Mask non-finite MAF or empty allele entries
valid_mask <- is.finite(maf_mat) & minor_mat != ""
maf_mat[!valid_mask] <- NA_real_

# Vectorized rowwise argmax — replace NA with -Inf for comparison only
maf_for_maxcol <- maf_mat
maf_for_maxcol[is.na(maf_for_maxcol)] <- -Inf
best_col <- max.col(maf_for_maxcol, ties.method = "first")

# Rows where all projects have no data
all_na <- rowSums(!is.na(maf_mat)) == 0L
row_idx <- seq_len(nrow(maf_mat))

commonSNP[, max_MAF := maf_mat[cbind(row_idx, best_col)]]
commonSNP[, max_ALT := minor_mat[cbind(row_idx, best_col)]]
commonSNP[all_na, `:=`(max_MAF = NA_real_, max_ALT = NA_character_)]


# --- Build clean tables -------------------------------------------------------

commonSNP_clean <- commonSNP[!is.na(max_MAF), .(CHROM, POS, ID, REF, ALT = max_ALT, MAF = max_MAF)]
commonSNP_filtered <- commonSNP_clean[MAF >= MAF_THRESHOLD]


# --- Build VCF-formatted tables -----------------------------------------------

grch38_header <- c(
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    "##source=UCSC_dbSnp155Common_Export",
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chr2,length=242193529>",
    "##contig=<ID=chr3,length=198295559>",
    "##contig=<ID=chr4,length=190214555>",
    "##contig=<ID=chr5,length=181538259>",
    "##contig=<ID=chr6,length=170805979>",
    "##contig=<ID=chr7,length=159345973>",
    "##contig=<ID=chr8,length=145138636>",
    "##contig=<ID=chr9,length=138394717>",
    "##contig=<ID=chr10,length=133797422>",
    "##contig=<ID=chr11,length=135086622>",
    "##contig=<ID=chr12,length=133275309>",
    "##contig=<ID=chr13,length=114364328>",
    "##contig=<ID=chr14,length=107043718>",
    "##contig=<ID=chr15,length=101991189>",
    "##contig=<ID=chr16,length=90338345>",
    "##contig=<ID=chr17,length=83257441>",
    "##contig=<ID=chr18,length=80373285>",
    "##contig=<ID=chr19,length=58617616>",
    "##contig=<ID=chr20,length=64444167>",
    "##contig=<ID=chr21,length=46709983>",
    "##contig=<ID=chr22,length=50818468>",
    "##contig=<ID=chrX,length=156040895>",
    "##contig=<ID=chrY,length=57227415>"
)

add_vcf_cols <- function(dt) {
    dt %>%
        dplyr::mutate(QUAL = ".", FILTER = "PASS", INFO = ".") %>%
        dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
}

vcf_clean <- add_vcf_cols(commonSNP_clean)
vcf_filtered <- add_vcf_cols(commonSNP_filtered)


# --- Save outputs -------------------------------------------------------------

muscadet::save_as_vcf(
    vcf_clean,
    file.path(output_dir, "hg38.dbSnp155Common.vcf.gz"),
    header = grch38_header
)

muscadet::save_as_vcf(
    vcf_filtered,
    file.path(output_dir, "hg38.dbSnp155Common.filtered.vcf.gz"),
    header = grch38_header
)

message("Done. Output files written to: ", output_dir)
