# muscadet 0.1.3

### Bug Fixes
- Fixed issue on methods applied to one-omic muscadet objects to get back to methods returning a list of objects instead of returning the object itself in the one-element list case.
- Fixed the table of features coordinates for ATAC in `computeLogRatioATAC()`, now the windows not covered by peaks are dropped.
- Fixed issue about nX number of chromosomes in `cnaCalling()` to take the maximum number of chromosomes depending on genome and not on the number of chromosomes in the object (case with missing chromosomes).

# muscadet 0.1.2

### Bug Fixes
- Seurat dependency set to `5.3.0` at minimum to avoid the installation of `leidenalg` python package when selecting the Leiden method for clustering.
- The `clusters` and `ncells` elements of the `@cnacalling` slot in muscadet objects are now updated during the `cnaCalling()` function right after some clusters need to be skipped due to not enough data. This resolves discordance between `clusters` and `ncells` elements with other `@cnacalling` elements, and avoids issues in plotting functions.
- The CNA calling functions (`assignClusters()`, `mergeCounts()` and `cnaCalling()` now all in `cnaCalling.R`) now include a reset of the `@cnacalling` slot of the muscadet object, to avoid discordance between results if clusters or parameters are changed in between the different steps of the CNA calling.

### Documentation
- Documentation added to `clusterMuscadet()` and `cluster_seurat()` for Leiden clustering algorithm recommendation when using the Seurat method.

# muscadet 0.1.1

### Bug Fixes
- Fix issue of `heatmapMuscadet()` crashing corner case: `heatmapMuscadet()` no longer crashes when `averages=TRUE` for one cluster only (monoclonal samples).
- Fix issue of `heatmapMuscadet()` with a wrong cluster order after Seurat clustering: `cluster_seurat()` and `imputeClusters()` now return sorted clusters vectors, from 1 to n, cells with imputed clusters are therefore placed at the end of each cluster group.
- Fix issue of wrong silhouette computation (due to clusters sorting changes): clusters are ordered according to the distance matrix order of cells used for silhouette computation, giving the correct silhouette scores.
- The `SeuratObject` package is now attached along with `muscadet`to easily use `Cells()` and `Features()` methods on muscadet objects.

# muscadet 0.1.0

## New Features
- Initial release of the muscadet package.
