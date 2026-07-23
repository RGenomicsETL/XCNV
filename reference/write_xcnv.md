# Write X-CNV-compatible output

Write a prediction data frame as a comma-separated X-CNV output table. A
file is written only when an explicit output path is supplied.

## Usage

``` r
write_xcnv(x, output)
```

## Arguments

- x:

  A data frame returned by
  [`predict_cnv()`](https://rgenomicsetl.github.io/XCNV/reference/predict_cnv.md).

- output:

  Explicit output path.

## Value

The normalized output path, invisibly.
