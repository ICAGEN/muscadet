# Add Labels to a ggplot Object

This function adds labels at the median position of each group to a
ggplot object. Labels can be added either as plain text or as label
boxes, with optional repelling to avoid overlaps (using the
[`ggprepel`](https://ggrepel.slowkow.com/reference/ggrepel.html) package
if installed).

## Usage

``` r
add_labels(
  p,
  labels,
  color = "black",
  repel = TRUE,
  label.box = FALSE,
  size = 3,
  ...
)
```

## Arguments

- p:

  A ggplot object with mapping and data (`ggplot`).

- labels:

  Column name (unquoted) indicating the group label to display.

- color:

  Color of the label text (`character`). Default is `"black"`.

- repel:

  Logical. If `TRUE`(default), overlapping labels are repelled using the
  [`ggprepel`](https://ggrepel.slowkow.com/reference/ggrepel.html)
  package.

- label.box:

  Logical. If `TRUE` it uses a boxed label (`geom_label`) instead of
  plain text (`geom_text`). Default is `FALSE`.

- size:

  Size of the label text (`numeric`). Default is `2`.

- ...:

  Additional arguments passed to the corresponding underlying geom:

  - [`ggplot2::geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
    (`repel`= `FALSE` and `label.box` = `FALSE`)

  - [`ggplot2::geom_label()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
    (`repel`= `FALSE` and `label.box` = `TRUE`)

  - [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
    (`repel`= `TRUE` and `label.box` = `FALSE`)

  - [`ggrepel::geom_label_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
    (`repel`= `TRUE` and `label.box` = `TRUE`)

## Value

A ggplot2 layer
([`geom_text`](https://ggplot2.tidyverse.org/reference/geom_text.html),
[`geom_label`](https://ggplot2.tidyverse.org/reference/geom_text.html),
[`geom_text_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html),
or
[`geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)).

## Details

The function summarizes the data by computing the median x and y
positions for each label group.

If `repel = TRUE`, it uses
[`geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
(`label.box = FALSE`) or
[`geom_label_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
(`label.box = TRUE`).

If `repel = FALSE`, it uses
[`geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
(`label.box = FALSE`) or
[`geom_label()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
(`label.box = TRUE`).

If `repel = TRUE`, the
[`ggprepel`](https://ggrepel.slowkow.com/reference/ggrepel.html) package
must be installed.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
library(ggrepel)

p <- ggplot(mtcars, aes(x = wt, y = mpg, color = as.factor(cyl))) +
  geom_point()

p2 <- add_labels(
  p,
  labels = cyl,
  repel = TRUE,
  label.box = FALSE,
  size = 5,
  min.segment.length = 0
)
p2

p3 <- add_labels(
    p,
    labels = cyl,
    repel = TRUE,
    label.box = TRUE,
    size = 3,
    box.padding = 1
)
p3

p4 <- add_labels(
    p,
    labels = cyl,
    repel = FALSE,
    label.box = TRUE,
    size = 4,
)
p4
} # }
```
