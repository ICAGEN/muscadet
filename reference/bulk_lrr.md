# Example data: Log R ratio from bulk sequencing data

Data frame containing log R ratio values per genomic segments from bulk
sequencing data.

## Usage

``` r
bulk_lrr
```

## Format

A data frame with the following columns:

- `CHROM`:

  Chromosome in integer format, e.g. 15, 23 (for X chromosome)
  (`character`).

- `start`:

  Start position of the segment (`integer`).

- `end`:

  End position of the segment (`character`).

- `lrr`:

  Log R ratio of the segment ("cnlr.median" column from
  [`facets::fitcncf()`](https://rdrr.io/pkg/facets/man/fitcncf.html)
  `$cncf` data frame) (`numeric`).

## Note

Data obtained from whole genome sequencing (WGS) after using
[`facets::fitcncf()`](https://rdrr.io/pkg/facets/man/fitcncf.html) from
[`facets`](https://rdrr.io/pkg/facets/man/facets-package.html) "Cellular
Fraction and Copy Numbers from Tumor Sequencing" version `0.6.2`:
`$cncf` data frame columns `chrom`, `start`, `end`, and `cnlr.median`.

## References

- [`facets`](https://rdrr.io/pkg/facets/man/facets-package.html)
  package:

  Shen R, Seshan VE. FACETS: allele-specific copy number and clonal
  heterogeneity analysis tool for high-throughput DNA sequencing.
  Nucleic Acids Res. 2016 Sep 19;44(16):e131. doi:
  [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520). PMID:
  27270079; PMCID: PMC5027494.
