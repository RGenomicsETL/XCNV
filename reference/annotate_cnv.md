# Annotate CNV records with X-CNV features

Compute the coding, genome-wide, regulatory, haploinsufficiency, and
population-frequency features used by X-CNV. Overlap is computed
in-process using BED half-open interval arithmetic; the reciprocal 0.5
thresholds used by the upstream frequency step are retained. CNV lengths
use the upstream inclusive `end - start + 1` convention.

## Usage

``` r
annotate_cnv(x, resources, overlap_backend = c("duckdb", "reference"))
```

## Arguments

- x:

  A path, data frame, or matrix containing CNV records.

- resources:

  An `xcnv_resources` object or an explicit resource directory.

- overlap_backend:

  `"duckdb"` for the production inequality-join path, or `"reference"`
  for the small brute-force correctness oracle.

## Value

A data frame in legacy feature order, with the input columns first, the
model and auxiliary annotation columns next, and `Type.1` and `Length`
at the end.
