# Validate CNV records

Validate the first four columns of a CNV table using the input rules of
the pinned X-CNV executable. Chromosome names may use the executable's
leading `chr` convention, and coordinates are retained as integer
values.

## Usage

``` r
validate_cnv(x)
```

## Arguments

- x:

  A path, data frame, or matrix containing chromosome, start, end, and
  CNV type in its first four columns.

## Value

`TRUE` invisibly when validation succeeds; otherwise an ordinary R error
is raised.
