# Derive objective ACMG/ClinGen 2019 CNV content evidence

Derive only the criteria determined mechanically by gene content and
established ClinGen dosage records: section 1 content, section 2A
complete dosage-feature containment, and section 3 protein-coding gene
count. Case, phenotype, segregation, population, partial-gene, and
inheritance evidence are intentionally not inferred. Inputs are BED
half-open relations and must carry source identifiers and releases.

## Usage

``` r
derive_cnv_2019_evidence(
  cnvs,
  genes,
  dosage,
  functional_biotypes = "protein_coding",
  established_score = 3
)
```

## Arguments

- cnvs:

  Data frame with `cnv_id`, `chr`, `start`, `end`, and `cnv_type`
  (`"loss"` or `"gain"`).

- genes:

  Gene relation with `chr`, `start`, `end`, `gene_id`, `biotype`,
  `source_id`, and `source_release`.

- dosage:

  ClinGen dosage relation with `chr`, `start`, `end`, `record_id`,
  `hi_score`, `ts_score`, `source_id`, and `source_release`.

- functional_biotypes:

  Biotypes counted as functionally important for sections 1 and 3;
  defaults to `"protein_coding"`.

- established_score:

  ClinGen score treated as established; defaults to 3.

## Value

A provenance-complete evidence data frame accepted by
[`score_cnv_2019()`](https://rgenomicsetl.github.io/XCNV/reference/score_cnv_2019.md).
