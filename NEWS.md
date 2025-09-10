# muscadet 0.1.1

### Bug Fixes

- Fix issue of `heatmapMuscadet()` crashing corner case: `heatmapMuscadet()` no longer crashes when `averages=TRUE` for one cluster only (monoclonal samples).
- Fix issue of `heatmapMuscadet()` with a wrong cluster order after Seurat clustering: `cluster_seurat()` and `imputeClusters()` now return sorted clusters vectors, from 1 to n, cells with imputed clusters are therefore placed at the end of each cluster group.
- Fix issue of wrong silhouette computation (due to clusters sorting changes): clusters are ordered according to the distance matrix order of cells used for silhouette computation, giving the correct silhouette scores.
- The `SeuratObject` package is now attached along with `muscadet`to easily use `Cells()` and `Features()` methods on muscadet objects.

# muscadet 0.1.0

## New Features
- Initial release of the muscadet package.
