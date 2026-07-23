# Read and normalize a CNV table

Read a tab-separated BED-like CNV table or normalize an in-memory table.
The first four columns are used, matching the published command-line
tool; additional columns are ignored.

## Usage

``` r
read_cnv(x)
```

## Arguments

- x:

  A path, data frame, or matrix containing chromosome, start, end, and
  CNV type in its first four columns.

## Value

A four-column data frame with columns `chr`, `start`, `end`, and `type`,
classed as `xcnv_cnv`.
