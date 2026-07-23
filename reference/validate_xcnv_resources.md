# Validate an X-CNV resource bundle

Check that a local directory contains the named X-CNV annotation and/or
model resources. Validation never downloads files and never changes the
working directory.

## Usage

``` r
validate_xcnv_resources(path = NULL, require = "annotations")
```

## Arguments

- path:

  Resource directory. If omitted, `XCNV_RESOURCE_DIR` must be set.

- require:

  Resource groups to require: `"annotations"` (the default), `"model"`,
  or `"all"`. A model file is required only when `"model"` or `"all"` is
  requested; otherwise prediction uses the bundled model.

## Value

An object of class `xcnv_resources`, containing the normalized root and
resolved resource paths.
