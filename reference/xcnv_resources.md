# Discover an X-CNV resource bundle

Resolve and validate a resource directory. With no `path`, only the
explicit `XCNV_RESOURCE_DIR` environment variable is consulted; the
current working directory is never searched implicitly.

## Usage

``` r
xcnv_resources(path = NULL, require = "annotations")
```

## Arguments

- path:

  Resource directory, or `NULL` to use `XCNV_RESOURCE_DIR`.

- require:

  Resource groups to require.

## Value

An object of class `xcnv_resources`.
