# Example data: muscadet objects

[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects, containing two single-cell omic datasets: scATAC-seq and
scRNA-seq.

- `exdata_muscadet` with tumor cells data: sample cells

- `exdata_muscadet_ref` with normal cells data: reference cells

## Usage

``` r
exdata_muscadet

exdata_muscadet_ref
```

## Format

[`muscadet`](https://icagen.github.io/muscadet/reference/muscadet-class.md)
objects with the following slots:

- `omics`:

  List of
  [`muscomic`](https://icagen.github.io/muscadet/reference/muscomic-class.md)
  objects, one per single-cell omic (`list`).

- `bulk.data`:

  List of data from paired bulk sequencing (`list`).

- `clustering`:

  List of data associated with the clustering of cells based on coverage
  log R ratio values (`list`).

- `cnacalling`:

  List of data associated with the calling of copy number alterations
  (CNAs) (`list`).

- `genome`:

  Reference genome name among: "hg38", "hg19" and "mm10" (`character`).
