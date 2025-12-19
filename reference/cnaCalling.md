# Copy Number Alteration (CNA) Calling from muscadet object

Performs copy number alteration (CNA) analysis on a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object by processing allelic and coverage counts across clusters and
evaluating cell fractions and copy numbers.

## Usage

``` r
cnaCalling(
  x,
  omics.coverage = NULL,
  depthmin.a.clusters = 30,
  depthmin.c.clusters = 50,
  depthmin.a.allcells = 30,
  depthmin.c.allcells = 50,
  depthmin.c.nor = 1000,
  depthmax.nor = NULL,
  het.thresh = 0.25,
  snp.nbhd = 250,
  hetscale = TRUE,
  cval1 = 25,
  cval2 = 150,
  min.nhet = 5,
  clonal.thresh = 0.9,
  cf.thresh = 0.5,
  dist.breakpoints = 1e+06,
  minoverlap = 1e+06,
  ploidy = "auto",
  dipLogR = NULL,
  quiet = FALSE
)
```

## Source

This function uses several functions from the
[`facets`](https://rdrr.io/pkg/facets/man/facets-package.html) package,
including:

- [`facets::clustersegs()`](https://rdrr.io/pkg/facets/man/facets-internal.html)

- [`facets::emcncf()`](https://rdrr.io/pkg/facets/man/emcncf.html)

- [`facets::findDiploidLogR()`](https://rdrr.io/pkg/facets/man/facets-internal.html)

- [`facets::fitcncf()`](https://rdrr.io/pkg/facets/man/fitcncf.html)

- [`facets::procSample()`](https://rdrr.io/pkg/facets/man/procSample.html)

- [`facets::procSnps()`](https://rdrr.io/pkg/facets/man/facets-internal.html)

- adapted function
  [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)

Seshan VE, Shen R (2021). *facets: Cellular Fraction and Copy Numbers
from Tumor Sequencing*. R package version 0.6.2,
<https://github.com/mskcc/facets>.

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object. Must contain:

  - Clustering assignments in the `cnacalling$clusters` slot (use
    [`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md)).

  - Combined allelic and coverage counts per cluster in the
    `cnacalling$combined.counts` slot (use
    [`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md)).

- omics.coverage:

  A vector of omics names to select for coverage log R ratio data.
  RECOMMENDED: select "ATAC" when ATAC and RNA omics are available, the
  ATAC coverage (DNA) signal is less noisy than RNA signal. By default,
  `NULL` selects all available data.

- depthmin.a.clusters:

  Minimum allelic depth per clusters in tumor cells (default: 30).

- depthmin.c.clusters:

  Minimum coverage depth per clusters in tumor cells (default: 50).

- depthmin.a.allcells:

  Minimum allelic depth for all tumor cells (default: 30).

- depthmin.c.allcells:

  Minimum coverage depth for all tumor cells (default: 50).

- depthmin.c.nor:

  Minimum coverage depth for normal sample (default: 1000).

- depthmax.nor:

  Optional. Maximum depth for normal sample (default: `NULL`).

- het.thresh:

  VAF (Variant Allele Frequency) threshold to call variant positions
  heterozygous for
  [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)
  (default: 0.25).

- snp.nbhd:

  Window size for selecting SNP loci to reduce serial correlation for
  [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)
  (default: 250).

- hetscale:

  Logical value indicating whether log odds ratio (logOR) should be
  scaled to give more weight in the test statistics for segmentation and
  clustering
  [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md).
  (default: `TRUE`)

- cval1:

  Critical value for segmentation for
  [`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)
  (default: 25).

- cval2:

  Critical value for segmentation for
  [`facets::procSample()`](https://rdrr.io/pkg/facets/man/procSample.html)
  (default: 150).

- min.nhet:

  Minimum number of heterozygous positions in a segment for
  [`facets::procSample()`](https://rdrr.io/pkg/facets/man/procSample.html)
  and [`facets::emcncf()`](https://rdrr.io/pkg/facets/man/emcncf.html)
  (default: 5).

- clonal.thresh:

  Threshold of minimum cell proportion to label a segment as clonal
  (default: 0.9).

- cf.thresh:

  Numeric threshold to set the minimum cell fraction (CF) allowed. CNA
  segments with a CF below this value are considered neutral (i.e.,
  2:1). The default is `0.5`. Can be set to `0` or `NULL` to disable
  this threshold and retain all segments regardless of CF.

- dist.breakpoints:

  Minimum distance between breakpoints to define distinct segments (see
  [`getSegConsensus()`](https://icagen.github.io/muscadet/reference/getSegConsensus.md))
  (default: 1e6).

- minoverlap:

  Minimum distance between a cluster specific segment and a consensus
  segment for them to overlap (see
  [`annotateSegments()`](https://icagen.github.io/muscadet/reference/annotateSegments.md))
  (default: 1e6).

- ploidy:

  Specifies ploidy assumption: `"auto"`, `"median"`, or numeric value
  (default: `"auto"`).

- dipLogR:

  Optional. Numeric value specifying an expected log ratio for diploid
  regions to use for the analysis. If `NULL`, by default, the diploid
  log ratio is estimated automatically from the data.

- quiet:

  Logical. If `TRUE`, suppresses informative messages during execution.
  Default is `FALSE`.

## Value

A modified
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object with added CNA analysis results in the `cnacalling` slot,
including: filtered counts and positions, segmentation data for clusters
and all cells, consensus segments across clusters based on breakpoints,
diploid log R ratio, purity and ploidy.

Details of the cnacalling slot:

- `combined.counts.filtered`: Filtered counts per clusters.

- `combined.counts.allcells`: Counts summed for all cells (no cluster
  distinction).

- `combined.counts.allcells.filtered`: Filtered counts summed for all
  cells (no cluster distinction).

- `positions`: Data frame of positions from the per cluster analysis.
  Positions in rows and associated data in columns: `chrom`, `maploc`
  (position), `rCountT` (read count in tumor), `rCountN` (read count in
  normal), `vafT` (variant allele frequency in tumor), `vafN` (variant
  allele frequency in normal), `cluster` (cluster id), `signal` (whether
  the counts come from coverage or allelic data), `het` (heterozygous
  status), `keep` (whether to keep position), `gcpct` (GC percentage),
  `gcbias` (GC bias correction), `cnlr` (log R ratio), `valor` (log odds
  ratio), `lorvar` (variance of log odds ratio), `seg0`, `seg_ori`
  (segment original id within each cluster), `seg` (segment id),
  `segclust` (cluster of segments id), `vafT.allcells` (vairiant allele
  frequency in all tumor cells), `colVAR` (integer for allelic position
  color depending on `vafT.allcells`).

- `segments`: Data frame of segments from the per cluster analysis.
  Segments in rows and associated data in columns: `chrom`, `seg`
  (segment id), `num.mark` (number of positions in segment), `nhet`
  (number of heterezygous positions in segment), `cnlr.median` (segment
  log R ratio median), `mafR` (segment square of expected log odds
  ratio), vafT.median (segment variant allele frequency median),
  `cluster` (cluster id), `seg_ori` (segment original id within each
  cluster), `segclust` (cluster of segments id), `cnlr.median.clust`
  (segment cluster log R ratio median), `mafR.clust` (segment cluster
  square of expected log odds ratio), `cf` (cell fraction), `tcn` (total
  copy number), `lcn` (lower copy number), `start`, `end`, `cf.em` (cell
  fraction computed with EM algorithm), `tcn.em`, (total copy number
  computed with EM algorithm), `lcn.em` (lower copy number computed with
  EM algorithm).

- `positions.allcells`: Same as `positions` but from the all cells
  analysis.

- `segments.allcells`: Same as `segments` but from the all cells
  analysis.

- `consensus.segs`: Data frame of unique consensus segments across
  clusters, with the `cna` (`logical`) and `cna_clonal` (`logical`)
  information.

- `table`: Data frame of consensus segments across clusters with
  associated information per cluster in columns: `chrom`, `start`,
  `end`, `id`, `cluster`, `cf.em` (cell fraction computed with EM
  algorithm), `tcn.em` (total copy number computed with EM algorithm),
  `lcn.em` (lower copy number computed with EM algorithm), `ncells`
  (number of cells in cluster), `prop.cluster` (proportion of cells per
  cluster), `gnl` (gain;neutral;loss : 1;0;-1), `loh` (loss of
  heterozygosity status), `state` (state of segments), `cna` (whether
  the segment is a CNA), `cna_state` (state of CNA segments), `prop.tot`
  (proportion of cells with the same state per segment), `state_clonal`
  (state of the segment if its `prop.tot` is above `clonal.thresh`),
  `cna_clonal` (whether the segment is a clonal CNA), `cna_clonal_state`
  (state of clonal CNA segments).

- `ncells`: Vector of number of cells per cluster.

- `dipLogR.clusters`: Diploid log R ratio estimated during the per
  cluster analysis.

- `dipLogR.allcells`: Diploid log R ratio estimated during the all cells
  analysis.

- `purity.clusters`: Purity estimated during the per cluster analysis.

- `purity.allcells`: Purity estimated during the all cells analysis.

- `ploidy.clusters`: Ploidy estimated during the per cluster analysis.

- `ploidy.allcells`: Ploidy estimated during the all cells analysis.

- `ploidy`: Ploidy used for the CNA analysis.

## References

- [`facets`](https://rdrr.io/pkg/facets/man/facets-package.html)
  package:

  Shen R, Seshan VE. FACETS: allele-specific copy number and clonal
  heterogeneity analysis tool for high-throughput DNA sequencing.
  Nucleic Acids Res. 2016 Sep 19;44(16):e131. doi:
  [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
  27270079; PMCID: PMC5027494.

## See also

[`assignClusters()`](https://icagen.github.io/muscadet/reference/assignClusters.md),
[`aggregateCounts()`](https://icagen.github.io/muscadet/reference/aggregateCounts.md),[`preProcSample2()`](https://icagen.github.io/muscadet/reference/preProcSample2.md)

## Examples

``` r
library("facets")
#> Loading required package: pctGCdata

# Load example muscadet object
# data("exdata_muscadet")

exdata_muscadet <- cnaCalling(
    exdata_muscadet,
    omics.coverage = "ATAC",
    depthmin.a.clusters = 3, # set low thresholds for example data
    depthmin.c.clusters = 5,
    depthmin.a.allcells = 3,
    depthmin.c.allcells = 5,
    depthmin.c.nor = 0
)
#> - Analysis per cluster -
#> Initial number of positions: 2455
#> Initial number of allelic positions: 1313
#> Initial number of coverage positions: 1142
#> Integrating omics...
#> Filtering allelic positions: tumor depth >= 3 reads
#> Selecting coverage data from omic(s): ATAC
#> Filtering coverage positions: tumor depth >= 5 reads
#> Filtering coverage positions: normal depth >= 0 reads
#> Allelic positions kept: 126
#> Coverage positions kept: 154
#> Final number of positions: 280
#> Performing segmentation per cluster...
#> Finding diploid log R ratio on clusters...
#> Diploid log R ratio = -0.273
#> Computing cell fractions and copy numbers on clusters...
#> - Analysis on all cells -
#> Aggregating allelic counts of all cells...
#> Filtering allelic positions: tumor depth >= 3 reads
#> Allelic positions kept: 186
#> Aggregating coverage counts of all cells...
#> Filtering coverage positions: tumor depth >= 5 reads
#> Filtering coverage positions: normal depth >= 0 reads
#> Coverage positions kept: 150
#> Final number of positions: 336
#> Performing segmentation on all cells...
#> Computing cell fractions and copy numbers on all cells...
#> - Consensus segments accross clusters -
#> Finding consensus segments...
#> 5 consensus segments identified, 0 CNA segments identified
#> Done.
```
