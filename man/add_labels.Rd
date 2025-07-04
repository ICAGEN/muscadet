% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{add_labels}
\alias{add_labels}
\title{Add Labels to a ggplot Object}
\usage{
add_labels(
  p,
  labels,
  color = "black",
  repel = TRUE,
  label.box = FALSE,
  size = 3,
  ...
)
}
\arguments{
\item{p}{A ggplot object with mapping and data (\code{ggplot}).}

\item{labels}{Column name (unquoted) indicating the group label to display.}

\item{color}{Color of the label text (\code{character}). Default is \code{"black"}.}

\item{repel}{Logical. If \code{TRUE}(default), overlapping labels are repelled
using the \code{\link[ggrepel:ggrepel]{ggprepel}} package.}

\item{label.box}{Logical. If \code{TRUE} it uses a boxed label (\code{geom_label})
instead of plain text (\code{geom_text}). Default is \code{FALSE}.}

\item{size}{Size of the label text (\code{numeric}). Default is \code{2}.}

\item{...}{Additional arguments passed to the corresponding underlying geom:
\itemize{
\item \code{\link[ggplot2:geom_text]{ggplot2::geom_text()}} (\code{repel}= \code{FALSE} and \code{label.box} = \code{FALSE})
\item \code{\link[ggplot2:geom_text]{ggplot2::geom_label()}} (\code{repel}= \code{FALSE} and \code{label.box} = \code{TRUE})
\item \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_text_repel()}} (\code{repel}= \code{TRUE} and \code{label.box} = \code{FALSE})
\item \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_label_repel()}} (\code{repel}= \code{TRUE} and \code{label.box} = \code{TRUE})
}}
}
\value{
A ggplot2 layer  (\code{\link[ggplot2:geom_text]{geom_text}},
\code{\link[ggplot2:geom_text]{geom_label}},
\code{\link[ggrepel:geom_text_repel]{geom_text_repel}}, or
\code{\link[ggrepel:geom_text_repel]{geom_label_repel}}).
}
\description{
This function adds labels at the median position of each group to a ggplot
object. Labels can be added either as plain text or as label boxes, with
optional repelling to avoid overlaps (using the
\code{\link[ggrepel:ggrepel]{ggprepel}} package if installed).
}
\details{
The function summarizes the data by computing the median x and y
positions for each label group.

If \code{repel = TRUE}, it uses
\code{\link[ggrepel:geom_text_repel]{geom_text_repel()}} (\code{label.box = FALSE})
or \code{\link[ggrepel:geom_text_repel]{geom_label_repel()}} (\code{label.box =
TRUE}).

If \code{repel = FALSE}, it uses \code{\link[ggplot2:geom_text]{geom_text()}} (\code{label.box
= FALSE}) or \code{\link[ggplot2:geom_text]{geom_label()}} (\code{label.box = TRUE}).

If \code{repel = TRUE}, the \code{\link[ggrepel:ggrepel]{ggprepel}} package must be
installed.
}
\examples{
\dontrun{
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
}

}
