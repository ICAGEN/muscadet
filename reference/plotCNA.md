# Plot CNA segments across clusters from a muscadet object

This function visualizes copy number alteration (CNA) segments across
clusters based on data stored in a
[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
object. It displays CNAs for each clusters and scales the y-axis based
on the proportion of cells in each cluster.

## Usage

``` r
plotCNA(
  x,
  title = NULL,
  cna.colors = c(gain = "#EF6F6AFF", loss = "#6699CCFF", cnloh = "#44AA99FF"),
  cf.gradient = TRUE
)
```

## Arguments

- x:

  A
  [`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
  object containing CNA calling data to be visualized (generated using
  [`cnaCalling()`](https://icagen.github.io/muscadet/reference/cnaCalling.md)).

- title:

  An optional title for the plot. Default is `NULL`.

- cna.colors:

  A vector of 3 colors for CNA states: gain, loss, and cnloh (or named
  vector where names are "gain", "loss", and "cnloh" and the values are
  their respective colors). Default is
  `c("gain" = "#EF6F6AFF", "loss" = "#6699CCFF", "cnloh" = "#44AA99FF")`.

- cf.gradient:

  Logical. If `TRUE` adds a alpha transparency gradient on CNA color
  block based on the cell fraction. Default is `TRUE`.

## Value

A ggplot object representing the CNA segments plot.

## Examples

``` r
if (FALSE) { # \dontrun{
library("ggplot2")

# Load example muscadet object
# data("exdata_muscadet")

# Plot CNA segments
p <- plotCNA(exdata_muscadet, title = "Copy Number Alterations in Example Data")
p
ggsave(
    filename = file.path("CNAplot.png"),
    plot = p,
    width = 3000,
    height = 800,
    units = "px"
)
} # }
```
