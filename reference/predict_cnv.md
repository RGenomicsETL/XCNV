# Predict X-CNV MVP scores

Annotate CNVs and compute the X-CNV MVP score. With an upstream model
this evaluates a staged `xcnv-tree-v1` artifact in ordinary R.
Installing and running the package does not require XGBoost.

## Usage

``` r
predict_cnv(
  x,
  resources,
  model = NULL,
  overlap_backend = c("duckdb", "reference")
)
```

## Arguments

- x:

  A path, data frame, or matrix containing CNV records.

- resources:

  An `xcnv_resources` object or an explicit resource directory.

- model:

  A model object/path. If `NULL`, a model in `resources` is used when
  present, otherwise the bundled published model is used.

- overlap_backend:

  Passed to
  [`annotate_cnv()`](https://rgenomicsetl.github.io/XCNV/reference/annotate_cnv.md).

## Value

A data frame containing the legacy annotation columns and a final
`MVP_score` column.
