# Preparation of data for muscadet

# directory of files: dir

library(usethis)
library(Matrix)
library(SeuratObject)
library(facets)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(GenomicRanges)

set.seed(123)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input count matrices ---------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Access dgCMatrix from RangedSummarizedExperiment object with SummarizedExperiment::assay(se)

mat_counts_atac_tumor <- readRDS(file.path(dir, "sample.atac.peak-matrix.Rds"))
rownames(mat_counts_atac_tumor) <- 1:nrow(mat_counts_atac_tumor)

mat_counts_atac_ref <- readRDS(file.path(dir, "reference.atac.peak-matrix.Rds"))
rownames(mat_counts_atac_ref) <- 1:nrow(mat_counts_atac_ref)

seurat_tumor <- readRDS(file.path(dir, "sample.rna.filtered.seurat-object.Rds"))
mat_counts_rna_tumor <- seurat_tumor[["RNA"]]$counts
rownames(mat_counts_rna_tumor) <- Features(seurat_tumor)
colnames(mat_counts_rna_tumor) <- paste("samplename", colnames(mat_counts_rna_tumor), sep = "_")

seurat_ref <- readRDS(file.path(dir, "reference.rna.filtered.seurat-object.Rds"))
mat_counts_rna_ref <- seurat_ref[["RNA"]]$counts
rownames(mat_counts_rna_ref) <- Features(seurat_ref)
colnames(mat_counts_rna_ref) <- gsub("refname.", "refname", colnames(mat_counts_rna_ref))
mat_counts_rna_ref <- mat_counts_rna_ref[, !duplicated(colnames(mat_counts_rna_ref))]



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling cells and features -----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## lists of barcodes -----------------------------------------------------------

# tumor
barcodes_atac_tumor <- scan(file.path(dir, "sample.atac.barcodes.txt"), what = "character")
barcodes_rna_tumor <- scan(file.path(dir, "sample.rna.barcodes.txt"), what = "character")

barcodes_atac_tumor <- barcodes_atac_tumor[which(barcodes_atac_tumor %in% colnames(mat_counts_atac_tumor))]
barcodes_rna_tumor <- barcodes_rna_tumor[which(barcodes_rna_tumor %in% colnames(mat_counts_rna_tumor))]

common_barcodes_tumor <- intersect(barcodes_atac_tumor, barcodes_rna_tumor)
common_barcodes_tumor_sub <- sort(sample(common_barcodes_tumor, 84))

barcodes_atac_tumor <- sort(c(common_barcodes_tumor_sub, sample(
  setdiff(barcodes_atac_tumor, common_barcodes_tumor), 28
)))
barcodes_rna_tumor <- sort(c(common_barcodes_tumor_sub, sample(
  setdiff(barcodes_rna_tumor, common_barcodes_tumor), 35
)))


# reference
barcodes_atac_ref <- scan(file.path(dir, "reference.atac.barcodes.txt"), what = "character")
barcodes_rna_ref <- scan(file.path(dir, "reference.rna.barcodes.txt"), what = "character")

barcodes_atac_ref <- barcodes_atac_ref[which(barcodes_atac_ref %in% colnames(mat_counts_atac_ref))]
barcodes_rna_ref <- barcodes_rna_ref[which(barcodes_rna_ref %in% colnames(mat_counts_rna_ref))]

common_barcodes_ref <- intersect(barcodes_atac_ref, barcodes_rna_ref)
common_barcodes_ref_sub <- sort(sample(common_barcodes_ref, 78))

barcodes_atac_ref <- sort(c(common_barcodes_ref_sub, sample(
  setdiff(barcodes_atac_ref, common_barcodes_ref), 21
)))
barcodes_rna_ref <- sort(c(common_barcodes_ref_sub, sample(
  setdiff(barcodes_rna_ref, common_barcodes_ref), 19
)))




## features coordinates --------------------------------------------------------

# from peak calling
peaks <- read.delim(
  file.path(dir, "peaks_coordinates.tsv"),
  header = F,
  col.names = c("CHROM", "start", "end", "id")
)
peaks <- peaks[sort(sample(as.integer(rownames(peaks)), 1000)), ]
rownames(peaks) <- NULL

# from gtf file corresponding to genome
genes <- read.delim(
  file.path(dir, "genes_coordinates.tsv"),
  header = F,
  col.names = c("CHROM", "start", "end", "id")
)
genes <- genes[which(genes[, "id"] %in% intersect(
  rownames(mat_counts_rna_tumor),
  rownames(mat_counts_rna_ref)
)), ]
rownames(genes) <- NULL
genes <- genes[sort(sample(as.integer(rownames(genes)), 500)), ]

usethis::use_data(peaks, overwrite = TRUE)
usethis::use_data(genes, overwrite = TRUE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling raw count matrices -----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_counts_atac_tumor <- mat_counts_atac_tumor[rownames(peaks), barcodes_atac_tumor]
mat_counts_atac_ref <- mat_counts_atac_ref[rownames(peaks), barcodes_atac_ref]
rownames(mat_counts_atac_tumor) <- peaks[, "id"]
rownames(mat_counts_atac_ref) <- peaks[, "id"]

mat_counts_rna_tumor <- mat_counts_rna_tumor[genes[, "id"], barcodes_rna_tumor]
mat_counts_rna_ref <- mat_counts_rna_ref[genes[, "id"], barcodes_rna_ref]

usethis::use_data(mat_counts_atac_tumor, overwrite = TRUE)
usethis::use_data(mat_counts_atac_ref, overwrite = TRUE)

usethis::use_data(mat_counts_rna_tumor, overwrite = TRUE)
usethis::use_data(mat_counts_rna_ref, overwrite = TRUE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling allele count tables ----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The variants REF and ALT alleles have been previously randomized for example data
allele_counts_atac_tumor <- read.delim(file.path(dir, "sample.atac.allele_counts.tsv"))
allele_counts_atac_ref <- read.delim(file.path(dir, "reference.atac.allele_counts.tsv"))
allele_counts_rna_tumor <- read.delim(file.path(dir, "sample.rna.allele_counts.tsv"))
allele_counts_rna_ref <- read.delim(file.path(dir, "reference.rna.allele_counts.tsv"))

colnames(allele_counts_atac_tumor)[which(colnames(allele_counts_atac_tumor) == "snp_id")] <- "id"
colnames(allele_counts_atac_ref)[which(colnames(allele_counts_atac_ref) == "snp_id")] <- "id"
colnames(allele_counts_rna_tumor)[which(colnames(allele_counts_rna_tumor) == "snp_id")] <- "id"
colnames(allele_counts_rna_ref)[which(colnames(allele_counts_rna_ref) == "snp_id")] <- "id"

# subsampling by cell barcodes
allele_counts_atac_tumor <- allele_counts_atac_tumor[which(allele_counts_atac_tumor[, "cell"] %in% barcodes_atac_tumor), ]
allele_counts_atac_ref <- allele_counts_atac_ref[which(allele_counts_atac_ref[, "cell"] %in% barcodes_atac_ref), ]

allele_counts_rna_tumor <- allele_counts_rna_tumor[which(allele_counts_rna_tumor[, "cell"] %in% barcodes_rna_tumor), ]
allele_counts_rna_ref <- allele_counts_rna_ref[which(allele_counts_rna_ref[, "cell"] %in% barcodes_rna_ref), ]


# subsampling by common variant identifiers

var_atac <- intersect(unique(allele_counts_atac_tumor[, "id"]),
                      unique(allele_counts_atac_ref[, "id"]))

var_atac_tumor <- c(var_atac, sample(
    setdiff(unique(allele_counts_atac_tumor[, "id"]), var_atac), 100))

var_atac_ref <- c(var_atac, sample(
    setdiff(unique(allele_counts_atac_ref[, "id"]), var_atac), 100))

var_rna <- intersect(unique(allele_counts_rna_tumor[, "id"]),
                      unique(allele_counts_rna_ref[, "id"]))

var_rna_tumor <- c(var_rna, sample(
    setdiff(unique(allele_counts_rna_tumor[, "id"]), var_rna), 100))

var_rna_ref <- c(var_rna, sample(
    setdiff(unique(allele_counts_rna_ref[, "id"]), var_rna), 100))

allele_counts_atac_tumor <- allele_counts_atac_tumor[which(allele_counts_atac_tumor[, "id"] %in% var_atac_tumor), ]
allele_counts_atac_tumor <- allele_counts_atac_tumor[order(allele_counts_atac_tumor[, "POS"]), ]
allele_counts_atac_tumor <- allele_counts_atac_tumor[order(allele_counts_atac_tumor[, "CHROM"]), ]

allele_counts_rna_tumor <- allele_counts_rna_tumor[which(allele_counts_rna_tumor[, "id"] %in% var_rna_tumor), ]
allele_counts_rna_tumor <- allele_counts_rna_tumor[order(allele_counts_rna_tumor[, "POS"]), ]
allele_counts_rna_tumor <- allele_counts_rna_tumor[order(allele_counts_atac_tumor[, "CHROM"]), ]

allele_counts_atac_ref <- allele_counts_atac_ref[which(allele_counts_atac_ref[, "id"] %in% var_atac_ref), ]
allele_counts_atac_ref <- allele_counts_atac_ref[order(allele_counts_atac_ref[, "POS"]), ]
allele_counts_atac_ref <- allele_counts_atac_ref[order(allele_counts_atac_ref[, "CHROM"]), ]

allele_counts_rna_ref <- allele_counts_rna_ref[which(allele_counts_rna_ref[, "id"] %in% var_rna_ref), ]
allele_counts_rna_ref <- allele_counts_rna_ref[order(allele_counts_rna_ref[, "POS"]), ]
allele_counts_rna_ref <- allele_counts_rna_ref[order(allele_counts_atac_ref[, "CHROM"]), ]

usethis::use_data(allele_counts_atac_tumor, overwrite = TRUE)
usethis::use_data(allele_counts_atac_ref, overwrite = TRUE)

usethis::use_data(allele_counts_rna_tumor, overwrite = TRUE)
usethis::use_data(allele_counts_rna_ref, overwrite = TRUE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Bulk log ratio result --------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# result after fitcncf() from FACETS 0.6.2 "Cellular Fraction and Copy Numbers from Tumor Sequencing"
wgs <- readRDS(file.path(dir, "sample.wgs.facets_result.Rds"))
bulk_lrr <- wgs$cncf[, c("chrom", "start", "end", "cnlr.median")]
colnames(bulk_lrr) <- c("CHROM", "start", "end", "lrr")

usethis::use_data(bulk_lrr, overwrite = TRUE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genome chromosome sizes ------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genome_list <- list()
for (genome in c("hg38", "hg19", "mm10")) {

    bs <- BSgenome::getBSgenome(genome)
    genome_chrom <- GenomicRanges::GRanges(names(seqlengths(bs)), IRanges(1, seqlengths(bs)))
    genome_chrom <- GenomeInfoDb::keepStandardChromosomes(genome_chrom, pruning.mode = "coarse")
    seqlevels(genome_chrom) <- stringr::str_remove(seqlevels(genome_chrom), "chr")
    genome_chrom <- GenomeInfoDb::dropSeqlevels(genome_chrom, "M", pruning.mode = "coarse")

    genome_list[[genome]] <- genome_chrom
}

hg38_chrom <- genome_list[["hg38"]]
hg19_chrom <- genome_list[["hg19"]]
mm10_chrom <- genome_list[["mm10"]]

usethis::use_data(hg38_chrom, hg19_chrom, mm10_chrom, overwrite = TRUE, internal = TRUE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# muscadet object --------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library("muscadet")

# muscomic objects
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor,
  allele_counts = allele_counts_atac_tumor,
  features = peaks
)
rna <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_tumor,
  allele_counts = allele_counts_rna_tumor,
  features = genes
)
atac_ref <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_ref,
  allele_counts = allele_counts_atac_ref,
  features = peaks
)
rna_ref <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_ref,
  allele_counts = allele_counts_rna_ref,
  features = genes
)

# raw muscadet objects
muscadet <- CreateMuscadetObject(
    omics = list(atac, rna),
    bulk.lrr = bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38"
)
muscadet_ref <- CreateMuscadetObject(
    omics = list(atac_ref, rna_ref),
    genome = "hg38"
)

# compute log R ratios for ATAC
muscadet <- computeLogRatio(
    x = muscadet,
    reference = muscadet_ref,
    omic = "ATAC",
    method = "ATAC",
    minReads = 1,
    minPeaks = 1,
    remove.raw = F
)
# compute log R ratios for RNA
muscadet <- computeLogRatio(
    x = muscadet,
    reference = muscadet_ref,
    omic = "RNA",
    method = "RNA",
    refReads = 2,
    remove.raw = F
)

# clustering
set.seed(123)
muscadet <- clusterMuscadet(
    muscadet,
    method = "seurat",
    res_range = c(0.6, 0.8, 1),
    dims_list = list(1:10, 1:10),
    knn_seurat = 10,
    knn_range_seurat = 30
)

# assign clusters for partition res=0.6
muscadet <- assignClusters(muscadet, partition = 0.6)

# merge counts
muscadet <- mergeCounts(muscadet, muscadet_ref)

# CNA calling
muscadet <- cnaCalling(
    muscadet,
    omics.coverage = "ATAC",
    depthmin.a.clusters = 3,
    depthmin.c.clusters = 5,
    depthmin.a.allcells = 3,
    depthmin.c.allcells = 5,
    depthmin.c.nor = 0
)


muscadet_obj <- muscadet
usethis::use_data(muscadet_obj, compress = "xz", overwrite = TRUE)
muscadet_obj_ref <- muscadet_ref
usethis::use_data(muscadet_obj_ref, compress = "xz", overwrite = TRUE)
