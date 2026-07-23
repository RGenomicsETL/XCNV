# Run X-CNV with a file-compatible interface

Read CNVs, compute annotations and MVP scores, and optionally write the
output table. Unlike the historical executable, this function never
writes beside the input by default and never installs software or
downloads data.

## Usage

``` r
run_xcnv(
  input,
  resources,
  model = NULL,
  output = NULL,
  overlap_backend = c("duckdb", "reference")
)

xcnv(
  input,
  resources,
  model = NULL,
  output = NULL,
  overlap_backend = c("duckdb", "reference")
)
```

## Arguments

- input:

  A CNV path, data frame, or matrix.

- resources:

  An `xcnv_resources` object or explicit resource directory.

- model:

  A model object/path, or `NULL` to use a resource-bundle override when
  present and otherwise the bundled published model.

- output:

  `NULL` to return results without writing, or one explicit CSV path.

- overlap_backend:

  Passed to
  [`annotate_cnv()`](https://rgenomicsetl.github.io/XCNV/reference/annotate_cnv.md).

## Value

A prediction data frame; the same object is returned invisibly when
`output` is non-`NULL`.
