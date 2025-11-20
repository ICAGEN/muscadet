# Process read count matrix and segmentation

Function adapted from
[`facets::preProcSample()`](https://rdrr.io/pkg/facets/man/preProcSample.html)
to process read count matrix and generates a segmentation tree. The
modifications from the original function includes:

- Incorporation of `cluster` and `signal` columns into the final result.

- Change in the log ratio (`cnlr`) for allelic data, it is computed as
  the log ratio (`cnlr`) mean between the previous and next coverage
  positions.

## Usage

``` r
preProcSample2(
  rcmat,
  het.thresh = 0.25,
  snp.nbhd = 250,
  cval = 25,
  gbuild = "hg38",
  hetscale = TRUE,
  ndepth = 5,
  ndepthmax = 1000
)
```

## Source

This function is derived from the
[`facets::preProcSample()`](https://rdrr.io/pkg/facets/man/preProcSample.html)
function, with the modifications to fit its use in muscadet.

Seshan VE, Shen R (2021). *facets: Cellular Fraction and Copy Numbers
from Tumor Sequencing*. R package version 0.6.2,
<https://github.com/mskcc/facets>.

## Arguments

- rcmat:

  A data frame with 8 required columns: `Chrom`, `Pos`, `NOR.DP`,
  `NOR.RD`, `TUM.DP`, `TUM.RD`, `cluster`, and `signal` (`data.frame`).

- het.thresh:

  VAF (Variant Allele Frequency) threshold to call variant positions
  heterozygous (`numeric`). Default: 0.25.

- snp.nbhd:

  Window size for selecting loci to reduce serial correlation
  (`numeric`). Default: 250.

- cval:

  Critical value for segmentation (`numeric`). Default: 25.

- gbuild:

  Genome build used for alignment. One of `"hg19"`, `"hg38"`, or
  `"mm10"` (`character`). Default: `"hg19"`.

- hetscale:

  Logical value indicating whether log odds ratio (logOR) should be
  scaled to give more weight in the test statistics for segmentation and
  clustering (`logical`). Usually only 10% of snps are hets and hetscale
  gives the logOR contribution to T-square as 0.25/proportion of hets.
  Default: `TRUE`.

- ndepth:

  Minimum depth in normal reference to keep (`numeric`). Default: 5.

- ndepthmax:

  Maximum normal coverage threshold for filtering loci (`numeric`).
  Default: 1000.

## Value

A list containing:

- pmat:

  Read counts and other elements of all the loci.

- gbuild:

  Genome build used for the analysis.

- nX:

  Chromosome number for X (e.g., 23 for human, 20 for mouse).

- clusters:

  Unique clusters from the processed data.

- seg.tree:

  Segmentation tree for each chromosome.

- jointseg:

  Segmented variant positions data.

- hscl:

  Scaling factor for logOR data.

## Details

The function processes variant positions data to generate a segmentation
tree. It uses
[`procSnps`](https://rdrr.io/pkg/facets/man/facets-internal.html) to
compute initial values, adjusts the log ratio for allelic signals, and
computes segmentation using
[`segsnps`](https://rdrr.io/pkg/facets/man/facets-internal.html).

SNPs (or variants) in a genome are not evenly spaced, and loci are
sampled within the specified window (`snp.nbhd`) to reduce serial
correlation.

## References

- [`facets`](https://rdrr.io/pkg/facets/man/facets-package.html)
  package:

  Shen R, Seshan VE. FACETS: allele-specific copy number and clonal
  heterogeneity analysis tool for high-throughput DNA sequencing.
  Nucleic Acids Res. 2016 Sep 19;44(16):e131. doi:
  [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
  27270079; PMCID: PMC5027494.

## See also

[`preProcSample`](https://rdrr.io/pkg/facets/man/preProcSample.html),
[`procSnps`](https://rdrr.io/pkg/facets/man/facets-internal.html),
[`segsnps`](https://rdrr.io/pkg/facets/man/facets-internal.html)

## Examples

``` r
library("facets")

# Load example muscadet object
# data("muscadet_obj")

counts <- muscadet_obj$cnacalling$combined.counts
counts <- counts[complete.cases(counts),]
counts_clus <- counts[which(counts$cluster == 1),]
result <- preProcSample2(counts_clus)
```
