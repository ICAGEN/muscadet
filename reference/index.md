# Package index

## Object creation and manipulation

### Objects, classes and methods

Classes and functions to create muscomic/muscadet objects and data
access

- [`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)
  :

  Create a `muscomic` object

- [`muscomic-class`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  : The muscomic class

- [`CreateMuscadetObject()`](https://icagen.github.io/muscadet/reference/CreateMuscadetObject.md)
  :

  Create a `muscadet` object

- [`muscadet-class`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  : The muscadet class

- [`coordFeatures()`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`matCounts()`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`matLogRatio()`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`Cells(`*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`Cells(`*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`Features(`*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  [`Features(`*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/muscadet-methods.md)
  :

  Methods for `muscomic` and `muscadet` objects

- [`` `[`( ``*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `[`( ``*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `$`( ``*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `$`( ``*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `[<-`( ``*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `[<-`( ``*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `$<-`( ``*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  [`` `$<-`( ``*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/musc-access.md)
  :

  Access and assignment methods for `muscadet` objects

- [`.DollarNames(`*`<muscadet>`*`)`](https://icagen.github.io/muscadet/reference/musc-auto.md)
  [`.DollarNames(`*`<muscomic>`*`)`](https://icagen.github.io/muscadet/reference/musc-auto.md)
  :

  Autocompletion for `$` access on `muscadet` or `muscomic` objects

- [`makeAllelicSparse()`](https://icagen.github.io/muscadet/reference/makeAllelicSparse.md)
  : Construct Sparse Matrices for Allelic Data

### Object incremental additions

Functions for incremental additions to muscadet objects during analysis

- [`addAlleleCounts()`](https://icagen.github.io/muscadet/reference/addAlleleCounts.md)
  :

  Add allele counts to a `muscadet` object

- [`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md)
  :

  Assign a cluster assignment to a `muscadet` object

## Log-ratio computation

Functions to compute log R ratios from the count matrices (coverage
information)

- [`computeLogRatio()`](https://icagen.github.io/muscadet/reference/computeLogRatio.md)
  : Compute log R ratios
- [`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md)
  : Compute log R ratios for scATAC-seq data
- [`computeLogRatioRNA()`](https://icagen.github.io/muscadet/reference/computeLogRatioRNA.md)
  : Compute log R ratios for scRNA-seq data
- [`getLogRatioBulk()`](https://icagen.github.io/muscadet/reference/getLogRatioBulk.md)
  : Retrieve log R ratio from Bulk data on single-cell features
  (internal)

## Integration and clustering

Functions to integrate log-ratios of omics and cluster cells

- [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)
  :

  Multi Omics Integration and Clustering on a `muscadet` Object

- [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  : Multi Omics Clustering using Seurat Multi Modal Graph-based
  Clustering

- [`cluster_hclust()`](https://icagen.github.io/muscadet/reference/cluster_hclust.md)
  : Multi Omics Clustering with SNF integration and Hierarchical
  Clustering

- [`weightedSNF()`](https://icagen.github.io/muscadet/reference/weightedSNF.md)
  : Weighted Similarity Network Fusion

- [`imputeClusters()`](https://icagen.github.io/muscadet/reference/imputeClusters.md)
  : Impute cluster assignments for missing cells by similarity

## Visualization - log-ratios and clusters

Functions to visualize log-ratios data on chromosomes, clusters of
cells, and usefull for clustering validation

- [`heatmapStep()`](https://icagen.github.io/muscadet/reference/heatmapStep.md)
  : Create heatmap and distribution plots of the different steps of
  computing log R ratios

- [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  :

  Heatmap plot for `muscadet` object

- [`plotSil()`](https://icagen.github.io/muscadet/reference/plotSil.md)
  :

  Silhouette plot for `muscadet` object

- [`plotIndexes()`](https://icagen.github.io/muscadet/reference/plotIndexes.md)
  :

  Plot clustering validation indexes for a `muscadet` object

- [`plotUMAP()`](https://icagen.github.io/muscadet/reference/plotUMAP.md)
  : Plot UMAP Coordinates from a muscadet Object

- [`add_labels()`](https://icagen.github.io/muscadet/reference/add_labels.md)
  : Add Labels to a ggplot Object

## Calling of Copy Number Alterations

Functions to identify Copy Number Alterations

- [`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)
  :

  Aggregate counts for `muscadet` objects

- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  : Copy Number Alteration (CNA) Calling from muscadet object

- [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)
  : Process read count matrix and segmentation

- [`getSegConsensus()`](https://icagen.github.io/muscadet/reference/getSegConsensus.md)
  : Get consensus segments across clusters

- [`annotateSegments()`](https://icagen.github.io/muscadet/reference/annotateSegments.md)
  : Annotate consensus segments across cluster with CNA states

## Visualization - Copy Number Alterations profiles

Functions to visualize Copy Number Alterations profiles at clone-level
and sample-level

- [`plotProfile()`](https://icagen.github.io/muscadet/reference/plotProfile.md)
  : Plot CNA profiles from muscadet object
- [`plotCNA()`](https://icagen.github.io/muscadet/reference/plotCNA.md)
  : Plot CNA segments across clusters from a muscadet object

## Data preparation

Functions facilitating allelic data preparation

- [`process_allele()`](https://icagen.github.io/muscadet/reference/process_allele.md)
  : Process allele counts results from SCReadCounts format
- [`save_as_vcf()`](https://icagen.github.io/muscadet/reference/save_as_vcf.md)
  : Save a data frame as a VCF File

## Example data

Example data used for examples and practicing functions

- [`exdata_muscadet`](https://icagen.github.io/muscadet/reference/exdata_muscadet.md)
  [`exdata_muscadet_ref`](https://icagen.github.io/muscadet/reference/exdata_muscadet.md)
  : Example data: muscadet objects
- [`exdata_mat_counts_atac_tumor`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.md)
  [`exdata_mat_counts_atac_ref`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.md)
  [`exdata_mat_counts_rna_tumor`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.md)
  [`exdata_mat_counts_rna_ref`](https://icagen.github.io/muscadet/reference/exdata_mat_counts.md)
  : Example data: Matrices of raw counts
- [`exdata_allele_counts_atac_tumor`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  [`exdata_allele_counts_atac_ref`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  [`exdata_allele_counts_rna_tumor`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  [`exdata_allele_counts_rna_ref`](https://icagen.github.io/muscadet/reference/exdata_allele_counts.md)
  : Example data: Allele counts at variant positions
- [`exdata_genes`](https://icagen.github.io/muscadet/reference/exdata_features.md)
  [`exdata_peaks`](https://icagen.github.io/muscadet/reference/exdata_features.md)
  : Example data: Feature coordinates
- [`exdata_bulk_lrr`](https://icagen.github.io/muscadet/reference/exdata_bulk_lrr.md)
  : Example data: Log R ratio from bulk sequencing data

## Misc

- [`genome_chrom`](https://icagen.github.io/muscadet/reference/genome_chrom.md)
  [`hg38_chrom`](https://icagen.github.io/muscadet/reference/genome_chrom.md)
  [`hg19_chrom`](https://icagen.github.io/muscadet/reference/genome_chrom.md)
  [`mm10_chrom`](https://icagen.github.io/muscadet/reference/genome_chrom.md)
  : Genome chromosome sizes (internal data)
