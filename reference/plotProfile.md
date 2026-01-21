# Plot CNA profiles from muscadet object

This function generates a multi-panel plot of copy number alteration
(CNA) profiles from a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object, including: log R ratios values, log odds ratio (or variant
allele frequency), copy numbers, CNA status and cell fractions.

## Usage

``` r
plotProfile(
  x,
  data,
  title = NULL,
  allelic.type = "vaf",
  point.cex = c(0.4, 0.5),
  chrom.colors = c("slategrey", "skyblue"),
  var.colors = c("peachpuff2", "paleturquoise3"),
  cn.colors = c("grey20", "brown2"),
  cna.colors = c(gain = "#EF6F6AFF", loss = "#6699CCFF", cnloh = "#44AA99FF"),
  cf.colors = c("white", "grey20", "bisque2"),
  dipLogR.color = c("magenta4"),
  seg.color = c("brown2")
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing CNA calling data to be visualized (generated using
  [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)).

- data:

  Either a cluster identifier to plot data of a cluster or "allcells" to
  plot data on all cells.

- title:

  An optional title for the plot. Default is `NULL`.

- allelic.type:

  A character string indicating the allelic metric to plot: "lor" for
  log odds ratio or "vaf" for variant allele frequency. Default is
  "vaf".

- point.cex:

  Numeric vector of length 1 or 2 specifying the size of points in the
  plots. If a single value is provided, it will be replicated for both
  plots. Default is `c(0.4, 0.5)`.

- chrom.colors:

  A character vector of length 2 defining alternating chromosome colors.
  Default is `c("slategrey", "skyblue")`.

- var.colors:

  A character vector of length 2 for variant positions point colors
  depending of variant allele frequency in all cells. Use "none" to use
  the alternating chromosome colors (defined by `chrom.colors`). Default
  is `c("peachpuff2", "paleturquoise3")`.

- cn.colors:

  A character vector of length 2 for total copy number and minor allele
  copy number segment colors. Default is `c("black", "brown2")`.

- cna.colors:

  A vector of 3 colors for CNA states: gain, loss, and cnloh (or named
  vector where names are "gain", "loss", and "cnloh" and the values are
  their respective colors). Default is
  `c("gain" = "#EF6F6AFF", "loss" = "#6699CCFF", "cnloh" = "#44AA99FF")`.

- cf.colors:

  A character vector of length 3 for cellular fraction gradient (of 10
  values): start color of the gradient, end color of the gradient, and
  color for normal diploid (depending on the ploidy). Default is
  `c("white", "steelblue", "bisque2")`.

- dipLogR.color:

  A character string for the diploid log R ratio line color. Default is
  "magenta4".

- seg.color:

  A character string for the color of segment medians. Default is
  "brown2".

## Value

A multi-panel plot of CNA profiles is produced.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example muscadet object
# data("exdata_muscadet")

# Plot profile for first cluster
pdf("CNAprofile_allcells.pdf", width = 15, height = 7.5) # Save as PDF
plotProfile(exdata_muscadet, data = "1", title = "Example dataset - cluster 1")
dev.off()

# Plot profile for all cells
pdf("CNAprofile_allcells.pdf", width = 15, height = 7.5) # Save as PDF
plotProfile(exdata_muscadet, data = "allcells", title = "Example dataset - all cells")
dev.off()
} # }
```
