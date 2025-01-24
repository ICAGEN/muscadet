% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mergeCounts}
\alias{mergeCounts}
\title{Merge counts for \code{muscadet} objects}
\usage{
mergeCounts(x, reference, nor.het = TRUE)
}
\arguments{
\item{x}{A \code{\link{muscadet}} object containing sample data (\code{muscadet}).
This object must include clustering assignments in the
\code{cnacalling$clusters} slot.}

\item{reference}{A \code{\link{muscadet}} object containing reference data
(\code{muscadet}).}

\item{nor.het}{A logical value to specify if normal reference allele counts
are modified to: total depth counts divided by 2, to set positions as
heterozygous in allelic data (e.g. when heterozygous SNP positions are retrieve based on
matched WGS data) before combining coverage and allelic data. Default is
\code{TRUE}.}
}
\value{
A modified \code{\link{muscadet}} object corresponding to the \code{x} muscadet object,
with updated \code{cnacalling} slot containing:
\itemize{
\item \code{allelic.counts}: Processed allelic counts on SNP positions, for all omics.
\item \code{coverage.counts}: Processed coverage counts merged with the reference.
\item \code{combined.counts}: Combined data for allelic and coverage counts.
}
Abbreviations:
\itemize{
\item RD = Reference allele read depth
\item AD = Alternative allele read depth
\item DP = Total read depth
\item TUM = tumor sample
\item NOR = normal reference
\item omic = omic specific (\code{omic} column)
\item all = for all omics
}
}
\description{
This function combines allelic (counts at single nucleotide polymorphisms
(SNPs)) and coverage counts (counts on features) from all omics per cluster
for both sample and reference. The resulting merged data is stored in the
\code{cnacalling} slot of the sample \code{muscadet} object.
}
\examples{
# Load example muscadet objects
data(muscadet_obj)
data(muscadet_obj_ref)

# Merge counts from all omics from both sample and reference
muscadet_obj <- mergeCounts(muscadet_obj, muscadet_obj_ref)

}
