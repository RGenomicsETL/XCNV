# Load an X-CNV model

Load a portable, tabular X-CNV tree model. The portable artifact is
exported from the pinned model during resource staging and is evaluated
without an XGBoost runtime. Serialized `.Rdata`, `.rds`, and native
XGBoost files are deliberately rejected; use `tools/stage_xcnv_model.R`
to convert them in an environment where XGBoost is installed.

## Usage

``` r
load_xcnv_model(model = NULL)
```

## Arguments

- model:

  A portable model object or path to an `xcnv-tree-v1` TSV file. The
  bundled published model is used when omitted or `NULL`.

## Value

A model object suitable for
[`predict_cnv()`](https://rgenomicsetl.github.io/XCNV/reference/predict_cnv.md).
