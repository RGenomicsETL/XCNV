# Score transparent ACMG/ClinGen 2019 constitutional CNV evidence

The input is an auditable evidence relation. Every row names its CNV,
deletion/gain context, criterion, numeric score, source, source release,
and evidence record. Repeated case evidence is capped by the published
criterion group maximum/minimum. Evidence selection remains explicit and
is never inferred from free text or phenotype similarity by this
function.

## Usage

``` r
score_cnv_2019(evidence, detail = c("summary", "groups"))
```

## Arguments

- evidence:

  Data frame with `cnv_id`, `cnv_type`, `criterion`, `score`,
  `source_id`, `source_release`, and `evidence_id`; optional logical
  `selected` defaults to `TRUE`.

- detail:

  `"summary"` for one row per CNV or `"groups"` for the capped
  criterion-group audit relation.

## Value

A data frame.
