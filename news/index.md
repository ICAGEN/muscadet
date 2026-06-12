# Changelog

## muscadet 0.2.2

#### Additions

- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md):
  add `filter.homozygous` argument (default: `TRUE`) to filter
  homozygous allelic positions in non-deletion regions prior to CNA
  calling, to reduce noise in the allelic imbalance signal. Designed for
  use cases where allele counts are derived from a common SNPs database
  rather than individual-specific heterozygous positions called from
  bulk sequencing data. Supporting arguments `vaf.thresh`,
  `vafmean.win`, and `vafmean.thresh` control the filtering behaviour.
- [`plotCNA()`](https://icagen.github.io/muscadet/reference/plotCNA.md):
  add `labels` argument to control cluster label display (`"auto"`,
  `"clusters"`, `"cells"`, or a custom character vector).
- [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  and
  [`heatmapStep()`](https://icagen.github.io/muscadet/reference/heatmapStep.md):
  add `dim_scale` argument to scale output dimensions and
  `raster_quality` argument to control heatmap tile rasterization
  quality. Both arguments allow finer control over output file size.
  Defaults preserve previous behaviour (`dim_scale = 1`,
  `raster_quality = 3`).

#### Bug fixes

- [`plotCNA()`](https://icagen.github.io/muscadet/reference/plotCNA.md):
  fix incorrect behaviour when clusters have identical proportions,
  caused by applying [`unique()`](https://rdrr.io/r/base/unique.html) to
  proportion values and an incorrect variable name.
- [`plotCNA()`](https://icagen.github.io/muscadet/reference/plotCNA.md):
  improve `cna.colors` argument validation, named vectors with missing
  or unknown state names now produce informative warnings and are filled
  with defaults rather than failing silently.
- [`plotCNA()`](https://icagen.github.io/muscadet/reference/plotCNA.md):
  set `cf.gradient` default to `FALSE`.
- [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md):
  fix incorrect number of cells displayed in the cluster annotation when
  a single cluster is present or when `averages = TRUE`.
- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md):
  fix unnamed cluster vector in `$cnacalling$clusters` slot when some
  clusters are skipped during CNA calling (due to low cell number),
  which compromised downstream usage.
- `createMuscadetObject()`: prevent silent introduction of `NA`
  chromosome names from bulk LRR data. Chromosome names in `bulk.lrr`
  are now validated against omics coordinate features; unmatched
  chromosomes are removed with a warning.
- [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md):
  fix missing `leiden_method` argument inheritance from
  [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md).

#### Documentation

- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md):
  update documentation with new `filter.homozygous` arguments and a
  dedicated details section describing the filtering logic.
- Update *get started* and *tutorial* vignettes to add usage notes for
  `filter.homozygous`depending on the allele counts input type.
- Update *preparation of input data* vignette to add details about the
  common SNP database, including MAF filtering recommendations, and
  download link to hg38 common SNP data (UCSC dbSNP 155 Common SNP
  filtered MAF \> 0.2) provided along the release.
- Add the script for processing and MAF filtering of hg38 common SNP
  data at `data-raw/hg38.dbSnp155Common.processing.R`.

## muscadet 0.2.1

#### Bug Fixes, Minor Updates

- Correct the [`show()`](https://rdrr.io/r/methods/show.html) method for
  muscadet objects.
- Update and fix
  [`plotProfile()`](https://icagen.github.io/muscadet/reference/plotProfile.md):
  `allele.type` argument default setting is now “vaf” ; mirror segment
  medians for “vaf” allelic type are fixed ; `lor.colors` argument
  becomes `var.colors` argument ; minor corrections of y axis labels
  rotation and x axis margins.

#### Documentation

- Add a new article to the documentation website, titled **“Multiomic
  copy number analysis tutorial”**, presenting a complete step-by-step
  workflow on real data. Note that the input data are not included with
  the package; all results in the tutorial are pre-computed.

## muscadet 0.2.0

#### Additions

- Changes in `muscomic` objects to optimize storage and computation, by
  avoiding long-format tables in favors of sparse matrices and by
  reworking the object structure. **muscadet and muscomic objects from
  previous versions are not supported by this version.**
  - [`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)
    now expects input `mat_counts` as **cells x features** (instead of
    features x cells).
  - `muscomic` object new structure:
    - **coverage slot** now contains `counts` and `logratio` layers,
      each with:
      - `mat`: matrix of counts or log ratio *cells × features* (sparse
        or dense depending on data type).
      - `coord.features`: coordinates of features matching the matrix.
      - `label.features`: feature type labels matching counts or log
        ratio.
    - **allelic slot** now contains 3 elements:
      - `mat.RD`: sparse matrix of reference allele counts *cells x
        features*.
      - `mat.AD`: sparse matrix of alternate allele counts *cells x
        features*.
      - `coord.vars`: coordinates of variant positions matching the
        matrices. Implemented new function `makeAllellicSparse()` to
        transform of long-format table (input) into sparse matrices.
        Updated
        [`addAlleleCounts()`](https://icagen.github.io/muscadet/reference/addAlleleCounts.md)
        [`CreateMuscomicObject()`](https://icagen.github.io/muscadet/reference/CreateMuscomicObject.md)
        for sparse format.
  - All package functions have been updated consequently to support the
    new format.
  - Implemented
    [`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)
    to replace `mergeCounts()` and support the new format.
  - Updated function examples, tests, and documentation consequently.
- [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)
  with Seurat-based clustering now supports `leiden_method` argument:
  - Default `"igraph"` implementation; `"leidenbase"` also available.
  - Seurat dependency set to **v5.3.1**.
  - Added `random.seed` argument for reproducible clustering (default
    set to 1).
- Updated example dataset objects `exdata`: subset of cells and subset
  of features/variant on 3 chromosomes for object size optimization.
  Updated function examples and vignettes consequently.

#### Deprecated and Defunct

- `mergeCounts()` removed in favor of
  [`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md).

#### Bug Fixes, Minor Updates

- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  refactored to handle separate aggregation and filtering of allelic and
  coverage counts.
- Fixed handling of clusters in character format for
  [`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)
  outputs.
- [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  `ndepth` argument default set to 1 to remove zero-depth allelic
  positions.
- `quiet` argument is now correctly propagated in
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md).
- Fix, move (from Imports to Suggests) and remove some dependencies.

## muscadet 0.1.3

#### Bug Fixes

- Fixed issue on methods applied to one-omic muscadet objects to get
  back to methods returning a list of objects instead of returning the
  object itself in the one-element list case.
- Fixed the table of features coordinates for ATAC in
  [`computeLogRatioATAC()`](https://icagen.github.io/muscadet/reference/computeLogRatioATAC.md),
  now the windows not covered by peaks are dropped.
- Fixed issue about nX number of chromosomes in
  [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  to take the maximum number of chromosomes depending on genome and not
  on the number of chromosomes in the object (case with missing
  chromosomes).

## muscadet 0.1.2

#### Bug Fixes

- Seurat dependency set to `5.3.0` at minimum to avoid the installation
  of `leidenalg` python package when selecting the Leiden method for
  clustering.
- The `clusters` and `ncells` elements of the `@cnacalling` slot in
  muscadet objects are now updated during the
  [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  function right after some clusters need to be skipped due to not
  enough data. This resolves discordance between `clusters` and `ncells`
  elements with other `@cnacalling` elements, and avoids issues in
  plotting functions.
- The CNA calling functions
  ([`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md),
  `mergeCounts()` and
  [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)
  now all in `cnaCalling.R`) now include a reset of the `@cnacalling`
  slot of the muscadet object, to avoid discordance between results if
  clusters or parameters are changed in between the different steps of
  the CNA calling.

#### Documentation

- Documentation added to
  [`clusterMuscadet()`](https://icagen.github.io/muscadet/reference/clusterMuscadet.md)
  and
  [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  for Leiden clustering algorithm recommendation when using the Seurat
  method.

## muscadet 0.1.1

#### Bug Fixes

- Fix issue of
  [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  crashing corner case:
  [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  no longer crashes when `averages=TRUE` for one cluster only
  (monoclonal samples).
- Fix issue of
  [`heatmapMuscadet()`](https://icagen.github.io/muscadet/reference/heatmapMuscadet.md)
  with a wrong cluster order after Seurat clustering:
  [`cluster_seurat()`](https://icagen.github.io/muscadet/reference/cluster_seurat.md)
  and
  [`imputeClusters()`](https://icagen.github.io/muscadet/reference/imputeClusters.md)
  now return sorted clusters vectors, from 1 to n, cells with imputed
  clusters are therefore placed at the end of each cluster group.
- Fix issue of wrong silhouette computation (due to clusters sorting
  changes): clusters are ordered according to the distance matrix order
  of cells used for silhouette computation, giving the correct
  silhouette scores.
- The `SeuratObject` package is now attached along with `muscadet`to
  easily use
  [`Cells()`](https://satijalab.github.io/seurat-object/reference/Cells.html)
  and
  [`Features()`](https://satijalab.github.io/seurat-object/reference/Cells.html)
  methods on muscadet objects.

## muscadet 0.1.0

### New Features

- Initial release of the muscadet package.
