<br>
<img src="https://img.shields.io/badge/version-dev-blue" alt="dev-version" align="left"/>
<br>

# muscadet <img src="man/figures/logo.png" align="right" height="139" alt="" />

### *multiomics single-cell copy number alterations detection*

**R package for the identification of copy number alterations (CNAs) in cancer cells from single-cell multiomics data.**

The package `muscadet` is designed for:

 * **Integration**: Integrate information from multiple omics (e.g. Multiome with matched scATAC-seq and scRNA-seq data).
 * **Clustering**: Cluster cells based on genome-wide coverage profiles.
 * **Imputation**: Impute clusters for cells missing data in one of the omic by nearest neighbor similarity.
 * **Detection**: Detect and call CNA segments using both coverage (log ratio of read counts) and allelic (read counts per allele) data.
 * **Visualization**: Visualize genome-wide coverage, clusters of cells, UMAP, CNA profiles, etc.
 
 

### Installation

You can install the development version of `muscadet` from GitHub:

```r
library(devtools)
devtools::install_github("ICAGEN/muscadet")
```

To get started, read the [Get Started vignette](https://icagen.github.io/muscadet/articles/muscadet.html) 
to learn about the `muscadet` workflow and function usage on demo example data.


<br>

***

### Detection of Somatic Copy Number Alterations from Single-Cell Multiomics Data

> Marie Denoulet<sup>1</sup>, Mia Cherkaoui<sup>1</sup>, Nils
Giordano<sup>1</sup>, Robin Lanée<sup>1</sup>, Elise Douillard<sup>1,2</sup>,
Magali Devic<sup>1,2</sup>, Florence Magrangeas<sup>1,2</sup>, Stéphane
Minvielle<sup>1,2</sup>, Céline Vallot<sup>3,4</sup>, Eric Letouzé<sup>1,2</sup>

> <sup>1</sup>Nantes Université, INSERM, CNRS, Université d'Angers, CRCI2NA,
Nantes, France. 
<sup>2</sup>University Hospital Hôtel-Dieu, Nantes, France.
<sup>3</sup>CNRS UMR3244, Institut Curie, PSL University, Paris, France.
<sup>4</sup>Translational Research Department, Institut Curie, PSL University,
Paris, France


The identification of somatic copy number alterations (CNAs) in cancer cells is
crucial for understanding tumor evolution, including clonal dynamics causing
relapse, and identifying potential therapeutic targets. While existing tools
provide valuable insights into subclonal CNAs, they are typically limited to
analyzing one type of omics data. In response to the growing use of cutting-edge
technologies enabling simultaneous sequencing of multiple omics from individual
cells, there emerges a need for new approaches that leverage multiomics data
integration to improve the detection of CNAs.
Addressing this need, we developed an R package, *muscadet*, that integrates
multiple single-cell datasets across different omics modalities to enhance the
accuracy and resolution of CNA detection within tumoral subclones. We
demonstrated the potency of our approach through the analysis of single-cell
Multiome data, integrating both single-cell RNA-seq and single-cell ATAC-seq
datasets from a common pool of cells, across several multiple myeloma samples.
*muscadet* outperformed existing copy number analysis tools for both scRNA-seq
and scATAC-seq data, revealing accurate CNA profiles and subclones, validated by
matched whole genome sequencing data.
By providing a unified CNA analysis framework applicable to any combination of
single-cell omics data, *muscadet* empowers researchers to unravel the clonal
structure of tumor samples and uncover complex genomic alterations driving
cancer progression.

***
