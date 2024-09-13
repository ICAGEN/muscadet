> [!NOTE]
> This R package is currently in development

<br>

# Detection of Somatic Copy Number Alterations from Single-Cell Multiomics Data

<img align="left" src="https://github.com/ICAGEN/muscadet/blob/main/logo.png" width="150" height="150" border="10"/> 

### *muscadet* R package
### *multiomics single-cell copy number alterations detection*

<br clear="left"/>

Marie Denoulet<sup>1,2,3</sup>, Mia Cherkaoui<sup>1</sup>, Nils Giordano<sup>1</sup>, Robin Lanée<sup>1</sup>, Elise Douillard<sup>1,4</sup>, Magali Devic<sup>1,4</sup>, Florence Magrangeas<sup>1,4</sup>, Stéphane Minvielle<sup>1,4</sup>, Céline Vallot<sup>5,6</sup>, Eric Letouzé<sup>1,4</sup>

<sup>1</sup>Nantes Université, INSERM, CNRS, Université d'Angers, CRCI2NA, Nantes, France. <sup>2</sup>SIRIC ILIAD, Angers, Nantes, France. <sup>3</sup>SIRIC Curie, Institut Curie, Paris, France. <sup>4</sup>University Hospital Hôtel-Dieu, Nantes, France. <sup>5</sup>CNRS UMR3244, Institut Curie, PSL University, Paris, France. <sup>6</sup>Translational Research Department, Institut Curie, PSL University, Paris, France

***

The identification of somatic copy number alterations (CNAs) in cancer cells is crucial for understanding tumor evolution, including clonal dynamics causing relapse, and identifying potential therapeutic targets. While existing tools provide valuable insights into subclonal CNAs, they are typically limited to analyzing one type of omics data. In response to the growing use of cutting-edge technologies enabling simultaneous sequencing of multiple omics from individual cells, there emerges a need for new approaches that leverage multiomics data integration to improve the detection of CNAs. Addressing this need, we developed a tool, *muscadet*, that integrates multiple single-cell datasets across different omics modalities to enhance the accuracy and resolution of CNA detection within tumoral subclones. We demonstrated the potency of our approach through the analysis of single-cell Multiome data, integrating both single-cell RNA-seq and single-cell ATAC-seq datasets from a common pool of cells, across 10 multiple myeloma samples. *muscadet* outperformed existing copy number analysis tools for both scRNA-seq and scATAC-seq data, revealing accurate CNA profiles and subclones, validated by matched whole genome sequencing data. By providing a unified CNA analysis framework applicable to any combination of single-cell omics data, *muscadet* empowers researchers to unravel the clonal structure of tumor samples and uncover complex genomic alterations driving cancer progression.

***
