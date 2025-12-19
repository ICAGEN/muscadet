# Example data: Feature coordinates

Data frames of features (peaks, genes...) coordinates on genome. Only
chromosomes 3, 4 and 8 are present in the `exdata` example dataset.

## Usage

``` r
exdata_genes

exdata_peaks
```

## Format

A data frame with the following columns:

- `CHROM`:

  Chromosome names, e.g. "3", "X".

- `start`:

  Start positions.

- `end`:

  End positions.

- `id`:

  Unique identifiers, e.g. gene name "CDH1" or peak identifier
  CHROM_start_end "1_1600338_1600838".
