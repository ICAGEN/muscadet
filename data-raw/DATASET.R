
# Preparation of data for muscadet package - example dataset exdata

library(usethis)
library(muscadet)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(GenomicRanges)
library(clustree)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genome chromosome sizes ------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genome_list <- list()
for (genome in c("hg38", "hg19", "mm10")) {

    bs <- BSgenome::getBSgenome(genome)
    genome_chrom <- GenomicRanges::GRanges(names(seqlengths(bs)), IRanges(1, seqlengths(bs)))
    genome_chrom <- GenomeInfoDb::keepStandardChromosomes(genome_chrom, pruning.mode = "coarse")
    seqlevels(genome_chrom) <- gsub("chr", "", seqlevels(genome_chrom))
    genome_chrom <- GenomeInfoDb::dropSeqlevels(genome_chrom, "M", pruning.mode = "coarse")

    genome_list[[genome]] <- genome_chrom
}

hg38_chrom <- genome_list[["hg38"]]
hg19_chrom <- genome_list[["hg19"]]
mm10_chrom <- genome_list[["mm10"]]

usethis::use_data(hg38_chrom, hg19_chrom, mm10_chrom, overwrite = TRUE, internal = TRUE)




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# exdata -----------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Set directory of files: dir <- "path/path"

set.seed(123)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Inputs -----------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Matrices of raw counts -------------------------------------------------------

# Access dgCMatrix from RangedSummarizedExperiment object with SummarizedExperiment::assay(se)
mat_counts_atac_tumor <- t(readRDS(file.path(dir, "sample.atac.peak-matrix.Rds")))
colnames(mat_counts_atac_tumor) <- 1:ncol(mat_counts_atac_tumor)

mat_counts_atac_ref <- t(readRDS(file.path(dir, "reference.atac.peak-matrix.Rds")))
colnames(mat_counts_atac_ref) <- 1:ncol(mat_counts_atac_ref)

seurat_tumor <- readRDS(file.path(dir, "sample.rna.filtered.seurat-object.Rds"))
mat_counts_rna_tumor <- t(seurat_tumor[["RNA"]]$counts)
colnames(mat_counts_rna_tumor) <- Features(seurat_tumor)
rownames(mat_counts_rna_tumor) <- paste("samplename", rownames(mat_counts_rna_tumor), sep = "_")

seurat_ref <- readRDS(file.path(dir, "reference.rna.filtered.seurat-object.Rds"))
mat_counts_rna_ref <- t(seurat_ref[["RNA"]]$counts)
colnames(mat_counts_rna_ref) <- Features(seurat_ref)
rownames(mat_counts_rna_ref) <- gsub("refname.", "refname", rownames(mat_counts_rna_ref))
mat_counts_rna_ref <- mat_counts_rna_ref[, !duplicated(rownames(mat_counts_rna_ref))]

rm(seurat_tumor, seurat_ref)

# Allele counts ----------------------------------------------------------------

# (variants positions have been previously randomized for example data)
allele_counts_atac_tumor <- read.delim(file.path(dir, "sample.atac.allele_counts.tsv"))
allele_counts_atac_ref <- read.delim(file.path(dir, "reference.atac.allele_counts.tsv"))
allele_counts_rna_tumor <- read.delim(file.path(dir, "sample.rna.allele_counts.tsv"))
allele_counts_rna_ref <- read.delim(file.path(dir, "reference.rna.allele_counts.tsv"))

## Features coordinates --------------------------------------------------------

# from peak calling
peaks <- read.delim(
    file.path(dir, "peaks_coordinates.tsv"),
    header = F,
    col.names = c("CHROM", "start", "end", "id")
)

# from gtf file corresponding to genome
genes <- read.delim(
    file.path(dir, "genes_coordinates.tsv"),
    header = F,
    col.names = c("CHROM", "start", "end", "id")
)
# remove genes not in matrices
genes <- genes[genes[, "id"] %in% intersect(colnames(mat_counts_rna_tumor),
                                            colnames(mat_counts_rna_ref)), ]

## Barcodes --------------------------------------------------------------------

b_atac_tumor <- scan(file.path(dir, "sample.atac.barcodes.txt"), what = "character")
b_rna_tumor <- scan(file.path(dir, "sample.rna.barcodes.txt"), what = "character")
b_atac_ref <- scan(file.path(dir, "reference.atac.barcodes.txt"), what = "character")
b_rna_ref <- scan(file.path(dir, "reference.rna.barcodes.txt"), what = "character")

# filter barcodes to match cells in matrices
b_atac_tumor <- b_atac_tumor[b_atac_tumor %in% rownames(mat_counts_atac_tumor)]
b_rna_tumor <- b_rna_tumor[b_rna_tumor %in% rownames(mat_counts_rna_tumor)]
b_atac_ref <- b_atac_ref[b_atac_ref %in% rownames(mat_counts_atac_ref)]
b_rna_ref <- b_rna_ref[b_rna_ref %in% rownames(mat_counts_rna_ref)]


# Bulk log ratio result --------------------------------------------------------

wgs <- readRDS(file.path(dir, "sample.wgs.facets_result.Rds"))
bulk_lrr <- wgs$cncf[, c("chrom", "start", "end", "cnlr.median")]
colnames(bulk_lrr) <- c("CHROM", "start", "end", "lrr")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling cells and features -----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Select features from chromosomes 3, 4 and 8 for the example
chroms <- c(3, 4, 8)

allele_counts_atac_tumor <- allele_counts_atac_tumor[allele_counts_atac_tumor$cell %in% b_atac_tumor, ]
allele_counts_rna_tumor <- allele_counts_rna_tumor[allele_counts_rna_tumor$cell %in% b_rna_tumor, ]
allele_counts_atac_ref <- allele_counts_atac_ref[allele_counts_atac_ref$cell %in% b_atac_ref, ]
allele_counts_rna_ref <- allele_counts_rna_ref[allele_counts_rna_ref$cell %in% b_rna_ref, ]

# Filter cells with allele counts data in selected chromosomes
allele_counts_tum <- lapply(list(
    allele_counts_atac_tumor,
    allele_counts_rna_tumor
), function(x) {
    x[x$CHROM %in% chroms, ]
})

allele_counts_ref <- lapply(list(
    allele_counts_atac_ref,
    allele_counts_rna_ref
), function(x) {
    x[x$CHROM %in% chroms, ]
})

tum_cells = Reduce(intersect, lapply(allele_counts_tum, `[[`, "cell"))
length(tum_cells)
ref_cells = Reduce(intersect, lapply(allele_counts_ref, `[[`, "cell"))
length(ref_cells)

# tumor
# random sample of cells
tum_cells_sub <- sort(sample(tum_cells, 63))
# add non-common cells
barcodes_atac_tumor <- sort(c(tum_cells_sub, sample(
    setdiff(b_atac_tumor, tum_cells), 8
)))
length(barcodes_atac_tumor)
barcodes_rna_tumor <- sort(c(tum_cells_sub, sample(
    setdiff(b_rna_tumor, tum_cells), 6
)))
length(barcodes_rna_tumor)

# reference
# random sample of cells
ref_cells_sub <- sort(sample(ref_cells, 85))
# add non-common cells
barcodes_atac_ref <- sort(c(ref_cells_sub, sample(
    setdiff(b_atac_ref, ref_cells), 10
)))
length(barcodes_atac_ref)
barcodes_rna_ref <- sort(c(ref_cells_sub, sample(
    setdiff(b_rna_ref, ref_cells), 8
)))
length(barcodes_rna_ref)

# Filter features
exdata_peaks <- rbind(
    peaks[peaks$id %in% sample(peaks[peaks$CHROM %in% chroms[1], "id"], 400), ],
    peaks[peaks$id %in% sample(peaks[peaks$CHROM %in% chroms[2], "id"], 400), ],
    peaks[peaks$id %in% sample(peaks[peaks$CHROM %in% chroms[3], "id"], 400), ]
)

exdata_genes <- rbind(
    genes[genes$id %in% sample(genes[genes$CHROM %in% chroms[1], "id"], 100), ],
    genes[genes$id %in% sample(genes[genes$CHROM %in% chroms[2], "id"], 100), ],
    genes[genes$id %in% sample(genes[genes$CHROM %in% chroms[3], "id"], 100), ]
)

# Filter bulk
exdata_bulk_lrr <- bulk_lrr[bulk_lrr$CHROM %in% chroms, ]
rownames(exdata_bulk_lrr) <- NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling raw count matrices -----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exdata_mat_counts_atac_tumor <- mat_counts_atac_tumor[barcodes_atac_tumor, rownames(exdata_peaks)]
exdata_mat_counts_atac_ref <- mat_counts_atac_ref[barcodes_atac_ref, rownames(exdata_peaks)]
colnames(exdata_mat_counts_atac_tumor) <- exdata_peaks$id
colnames(exdata_mat_counts_atac_ref) <- exdata_peaks$id

exdata_mat_counts_rna_tumor <- mat_counts_rna_tumor[barcodes_rna_tumor, exdata_genes$id]
exdata_mat_counts_rna_ref <- mat_counts_rna_ref[barcodes_rna_ref, exdata_genes$id]

rownames(exdata_peaks) <- NULL
rownames(exdata_genes) <- NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampling allele count tables ----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exdata_allele_counts_atac_tumor <- allele_counts_tum[[1]][allele_counts_tum[[1]][, "cell"] %in% barcodes_atac_tumor, ]
exdata_allele_counts_rna_tumor <- allele_counts_tum[[2]][allele_counts_tum[[2]][, "cell"] %in% barcodes_rna_tumor, ]

exdata_allele_counts_atac_ref <- allele_counts_ref[[1]][allele_counts_ref[[1]][, "cell"] %in% barcodes_atac_ref, ]
exdata_allele_counts_rna_ref <- allele_counts_ref[[2]][allele_counts_ref[[2]][, "cell"] %in% barcodes_rna_ref, ]


# subsampling by common variant identifiers

var_atac <- intersect(unique(exdata_allele_counts_atac_tumor[, "id"]), unique(exdata_allele_counts_atac_ref[, "id"]))
var_atac_tumor <- c(var_atac, sample(setdiff(unique(exdata_allele_counts_atac_tumor[, "id"]), var_atac), 10))
var_atac_ref <- c(var_atac, sample(setdiff(unique(exdata_allele_counts_atac_ref[, "id"]), var_atac), 10))

var_rna <- intersect(unique(exdata_allele_counts_rna_tumor[, "id"]), unique(exdata_allele_counts_rna_ref[, "id"]))
var_rna_tumor <- c(var_rna, sample(setdiff(unique(exdata_allele_counts_rna_tumor[, "id"]), var_rna), 10))
var_rna_ref <- c(var_rna, sample(setdiff(unique(exdata_allele_counts_rna_ref[, "id"]), var_rna), 10))

exdata_allele_counts_atac_tumor <- exdata_allele_counts_atac_tumor[exdata_allele_counts_atac_tumor[, "id"] %in% var_atac_tumor, ]
exdata_allele_counts_atac_tumor <- exdata_allele_counts_atac_tumor[order(exdata_allele_counts_atac_tumor$CHROM, exdata_allele_counts_atac_tumor$POS), ]

exdata_allele_counts_rna_tumor <- exdata_allele_counts_rna_tumor[exdata_allele_counts_rna_tumor[, "id"] %in% var_rna_tumor, ]
exdata_allele_counts_rna_tumor <- exdata_allele_counts_rna_tumor[order(exdata_allele_counts_rna_tumor$CHROM, exdata_allele_counts_rna_tumor$POS), ]

exdata_allele_counts_atac_ref <- exdata_allele_counts_atac_ref[exdata_allele_counts_atac_ref[, "id"] %in% var_atac_ref, ]
exdata_allele_counts_atac_ref <- exdata_allele_counts_atac_ref[order(exdata_allele_counts_atac_ref$CHROM, exdata_allele_counts_atac_ref$POS), ]

exdata_allele_counts_rna_ref <- exdata_allele_counts_rna_ref[exdata_allele_counts_rna_ref[, "id"] %in% var_rna_ref, ]
exdata_allele_counts_rna_ref <- exdata_allele_counts_rna_ref[order(exdata_allele_counts_rna_ref$CHROM, exdata_allele_counts_rna_ref$POS), ]

rownames(exdata_allele_counts_atac_tumor) <- NULL
rownames(exdata_allele_counts_rna_tumor) <- NULL
rownames(exdata_allele_counts_atac_ref) <- NULL
rownames(exdata_allele_counts_rna_ref) <- NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# muscadet object --------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# muscomic objects
exdata_atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = exdata_mat_counts_atac_tumor,
    allele_counts = exdata_allele_counts_atac_tumor,
    features = exdata_peaks
)
exdata_rna <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = exdata_mat_counts_rna_tumor,
    allele_counts = exdata_allele_counts_rna_tumor,
    features = exdata_genes
)
exdata_atac_ref <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = exdata_mat_counts_atac_ref,
    allele_counts = exdata_allele_counts_atac_ref,
    features = exdata_peaks
)
exdata_rna_ref <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = exdata_mat_counts_rna_ref,
    allele_counts = exdata_allele_counts_rna_ref,
    features = exdata_genes
)

# raw muscadet objects
exdata_muscadet <- CreateMuscadetObject(
    omics = list(exdata_atac, exdata_rna),
    bulk.lrr = exdata_bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38"
)
exdata_muscadet_ref <- CreateMuscadetObject(
    omics = list(exdata_atac_ref, exdata_rna_ref),
    genome = "hg38"
)

# compute log R ratios for ATAC
exdata_muscadet <- computeLogRatio(
    x = exdata_muscadet,
    reference = exdata_muscadet_ref,
    omic = "ATAC",
    method = "ATAC",
    minReads = 0.5,
    minPeaks = 1
)
# compute log R ratios for RNA
exdata_muscadet <- computeLogRatio(
    x = exdata_muscadet,
    reference = exdata_muscadet_ref,
    omic = "RNA",
    method = "RNA",
    refReads = 2
)

# clustering
set.seed(123)
exdata_muscadet <- clusterMuscadet(
    exdata_muscadet,
    method = "seurat",
    res_range = c(0.1, 0.3, 0.5),
    dims_list = list(1:10, 1:10),
    knn_seurat = 10,
    knn_range_seurat = 30
)

# clustree
partitions <- lapply(exdata_muscadet$clustering$clusters, as.data.frame)
partitions <- do.call(cbind, partitions)
colnames(partitions) <- paste0("res_", names(exdata_muscadet$clustering$clusters))
clustree(partitions, prefix = "res_")
# heatmap
heatmapMuscadet(exdata_muscadet, "exdata_heatmap_res0.3.png", partition = "0.3")

# assign clusters for partition
exdata_muscadet <- assignClusters(exdata_muscadet, partition = 0.3)

# merge counts
exdata_muscadet <- aggregateCounts(exdata_muscadet, exdata_muscadet_ref)

# CNA calling
exdata_muscadet <- cnaCalling(
    exdata_muscadet,
    depthmin.a.clusters = 3,
    depthmin.c.clusters = 5,
    depthmin.a.allcells = 3,
    depthmin.c.allcells = 5,
    depthmin.c.nor = 1
)

# plots
plotProfile(exdata_muscadet, data = 1, title = "Cluster 1 profile", point.cex = 0.8)
plotProfile(exdata_muscadet, data = 2, title = "Cluster 2 profile", point.cex = 0.8)
plotProfile(exdata_muscadet, data = "allcells", title = "All cells profile", point.cex = 0.8)
plotCNA(exdata_muscadet)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save objects -----------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# mat_counts
usethis::use_data(exdata_mat_counts_atac_tumor, overwrite = TRUE)
usethis::use_data(exdata_mat_counts_atac_ref, overwrite = TRUE)
usethis::use_data(exdata_mat_counts_rna_tumor, overwrite = TRUE)
usethis::use_data(exdata_mat_counts_rna_ref, overwrite = TRUE)

# allele_counts
usethis::use_data(exdata_allele_counts_atac_tumor, overwrite = TRUE)
usethis::use_data(exdata_allele_counts_atac_ref, overwrite = TRUE)
usethis::use_data(exdata_allele_counts_rna_tumor, overwrite = TRUE)
usethis::use_data(exdata_allele_counts_rna_ref, overwrite = TRUE)

# features
usethis::use_data(exdata_peaks, overwrite = TRUE)
usethis::use_data(exdata_genes, overwrite = TRUE)

# bulk data
usethis::use_data(exdata_bulk_lrr, overwrite = TRUE)

# muscadet objects
usethis::use_data(exdata_muscadet, compress = "xz", overwrite = TRUE)
usethis::use_data(exdata_muscadet_ref, compress = "xz", overwrite = TRUE)




