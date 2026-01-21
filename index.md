# muscadet

### *multiomics single-cell copy number alterations detection*

`muscadet` is an R package for identifying copy number alterations
(CNAs) in cancer cells from single-cell multiomics data.

#### Key features of `muscadet`:

- **Data**: Copy number analysis from a single omic to multiple omics at
  once on the same cells.
- **Integration**: Integrate information at low-level from the matched
  omics.
- **Clustering**: Cluster cells based on genome-wide coverage profiles.
- **Imputation**: Impute cluster for cells missing data in one of the
  omic by nearest neighbor similarity.
- **Detection**: Detect and call CNA segments using both coverage (log
  ratio of read counts) and allelic imbalance (read counts per allele).
- **Visualization**: Explore genome-wide coverage, clusters of cells,
  UMAP embeddings, CNA profiles, and more.

### Installation

You can install the latest version of `muscadet` from GitHub:

``` r
library(devtools)
devtools::install_github("ICAGEN/muscadet")
```

- To get started, read the [Get Started
  vignette](https://icagen.github.io/muscadet/articles/muscadet.html) to
  learn about the `muscadet` workflow and function usage on demo example
  data.
- To prepare inputs, read the [Preparation of input data
  vignette](https://icagen.github.io/muscadet/articles/data-prep.html)
  to get details on the inputs for `muscadet` and how to prepare them.

### Citation

To cite the `muscadet` package in publications, use:

> Denoulet M, Giordano N, Minvielle S, Vallot C, Letouzé E (2025).
> *muscadet: Multiomics Single-Cell Copy Number Alterations Detection*.
> R package version 0.2.0, <https://icagen.github.io/muscadet/>,
> <https://github.com/ICAGEN/muscadet>.

### Communications

- \[Poster\] **Detection of Somatic Copy Number Alterations from
  Single-cell Multiomics Data with the R Package muscadet** Marie
  Denoulet, Nils Giordano, Mia Cherkaoui, Elise Douillard, Magali Devic,
  Florence Magrangeas, Stéphane Minvielle, Céline Vallot, Eric Letouzé.
  *Journées Ouvertes en Biologie, Informatique et Mathématiques (JOBIM
  2025)*, Jul 2025, Bordeaux, France.
  [⟨hal-05150254⟩](https://nantes-universite.hal.science/hal-05150254v1)
- \[Poster\] **Detection of Somatic Copy Number Alterations from
  Single-cell Multiomics Data** Marie Denoulet, Nils Giordano, Mia
  Cherkaoui, Elise Douillard, Magali Devic, Florence Magrangeas,
  Stéphane Minvielle, Céline Vallot, Eric Letouzé. *23rd European
  Conference on Computational Biology (ECCB2024)*, Sep 2024, Turku,
  Finland.
  [⟨hal-05150028⟩](https://nantes-universite.hal.science/hal-05150254v1)

  

------------------------------------------------------------------------

### Detection of Somatic Copy Number Alterations from Single-Cell Multiomics Data with the R Package `muscadet`

> Marie Denoulet¹, Mia Cherkaoui¹, Nils Giordano¹, Robin Lanée¹, Elise
> Douillard^(1,2), Magali Devic^(1,2), Florence Magrangeas^(1,2),
> Stéphane Minvielle^(1,2), Céline Vallot^(3,4), Eric Letouzé^(1,2)

> ¹Nantes Université, INSERM, CNRS, Université d’Angers, CRCI2NA,
> Nantes, France. ²University Hospital Hôtel-Dieu, Nantes, France. ³CNRS
> UMR3244, Institut Curie, PSL University, Paris, France. ⁴Translational
> Research Department, Institut Curie, PSL University, Paris, France

The identification of somatic copy number alterations (CNAs) in cancer
cells is crucial for understanding tumor evolution, including clonal
dynamics causing relapse, and identifying potential therapeutic targets.
While existing tools provide valuable insights into subclonal CNAs, they
are typically limited to analyzing one type of omics data. In response
to the growing use of cutting-edge technologies enabling simultaneous
sequencing of multiple omics from individual cells, there emerges a need
for new approaches that leverage multiomics data integration to improve
the detection of CNAs.

Addressing this need, we developed an R package, *muscadet*, that
integrates multiple single-cell datasets across different omics
modalities to enhance the accuracy and resolution of CNA detection
within tumoral subclones. The *muscadet* framework enables integration
of joint multiomics data, such as single-cell transcriptomics
(scRNA-seq) and chromatin accessibility (scATAC-seq), captured from the
same cells. Cells are initially clustered by relative coverage (log
ratio of read counts) integrated across multiple omics layers. CNA
calling is then performed on aggregated coverage and allelic (read
counts per allele) data at the cluster level, leveraging combined
signals from omics and cells. Lastly, *muscadet* includes visualization
tools for genome-wide coverage assessment, subclone identification, and
comprehensive CNA profiling.

By providing a unified CNA analysis framework applicable to any
combination of single-cell omics data, *muscadet* empowers researchers
to unravel the clonal structure of tumor samples and uncover complex
genomic alterations driving cancer progression.

------------------------------------------------------------------------
