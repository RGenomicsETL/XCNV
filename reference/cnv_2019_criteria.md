# ACMG/ClinGen 2019 constitutional CNV criterion contract

Return the machine-readable criterion identifiers and numeric ranges
used by the independently authored scoring lane. This is a scoring
contract, not an automatic clinical evidence adjudicator. Dynamic
section 5 criteria have an `NA` default and require an explicit score
selected by the curator.

## Usage

``` r
cnv_2019_criteria(cnv_type = c("loss", "gain"))
```

## Arguments

- cnv_type:

  One or both of `"loss"` and `"gain"`.

## Value

A data frame with allowed per-evidence ranges and aggregate caps.
