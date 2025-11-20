# Example data: Feature coordinates

Data frames of features (peaks, genes...) coordinates on genome.

## Usage

``` r
genes

peaks
```

## Format

A data frame with the following columns:

- `CHROM`:

  Chromosome names in character format, e.g. "15", "X" (`character`).

- `start`:

  Start positions (`integer`).

- `end`:

  End positions (`character`).

- `id`:

  Unique identifiers, e.g. gene name "CDH1" or peak identifier
  CHROM_start_end "1_1600338_1600838" (`character`).
