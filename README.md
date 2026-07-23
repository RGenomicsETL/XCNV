XCNV
================

- [XCNV](#xcnv)
  - [Installation](#installation)
  - [A real ClinVar deletion](#a-real-clinvar-deletion)
  - [Inspectable ACMG/ClinGen
    evidence](#inspectable-acmgclingen-evidence)
  - [Published model without XGBoost](#published-model-without-xgboost)
  - [Full annotation resources](#full-annotation-resources)
  - [Citation](#citation)

<!-- README.md is rendered from this file. -->

# XCNV

XCNV annotates copy-number variants and evaluates the published X-CNV
classifier without bedtools or a runtime XGBoost dependency. It also
provides a separate ACMG/ClinGen 2019 lane in which every selected
criterion remains a row that can be inspected, changed, and rescored.

The implementation is table-oriented. DuckDB evaluates the genomic
overlap and reciprocal-overlap joins; ordinary R owns input validation,
feature reduction, the portable tree evaluator, and evidence scoring.

## Installation

``` r
install.packages(
  "XCNV",
  repos = c(
    "https://rgenomicsetl.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
```

## A real ClinVar deletion

The package includes a compact, case-limited annotation set for ClinVar
allele 72461: a pathogenic deletion containing `SOX2`. This is a real
scientific example, not the synthetic fixture used by the test suite.

``` r
library(XCNV)

example <- xcnv_real_example()
example$clinvar[, c(
  "clinvar_allele_id", "clinical_significance",
  "grch38_chrom", "grch38_start1", "grch38_end1",
  "grch37_chrom", "grch37_start0", "grch37_end0"
)]
#>   clinvar_allele_id clinical_significance grch38_chrom grch38_start1
#> 1             72461            Pathogenic            3     181710887
#>   grch38_end1 grch37_chrom grch37_start0 grch37_end0
#> 1   181748657            3     181428674   181466445
```

ClinVar publishes this record in both assemblies. X-CNV uses its GRCh37
half-open interval. When a source provides only one assembly, lift it
during resource staging—for example with Rduckhts’
`rduckhts_liftover()`—and retain the source and destination coordinates
beside the case.

The following call evaluates all five X-CNV annotation families and the
bundled published model:

``` r
prediction <- predict_cnv(example$cnv, example$resources)
prediction[, c(
  "Chr", "Start", "End", "Type",
  "CDTS_1st", "CDTS_5th",
  "gain.freq", "loss.freq",
  "pLI", "Episcore", "GHIS", "MVP_score"
)]
#>   Chr     Start       End Type   CDTS_1st  CDTS_5th gain.freq loss.freq
#> 1   3 181428674 181466445 loss 0.01032511 0.0765117         0         0
#>         pLI  Episcore      GHIS MVP_score
#> 1 0.3966965 0.9209667 0.6169144 0.8537559
```

`MVP_score` is the X-CNV model output; XCNV does not relabel it as an
ACMG criterion or clinical class.

## Inspectable ACMG/ClinGen evidence

The same real deletion completely contains the GENCODE v19 `SOX2` gene
and its ClinGen haploinsufficiency record. The automatic provider
derives only the objective content criteria it can establish from those
relations:

``` r
cnv <- transform(
  example$cnv,
  cnv_id = example$clinvar$cnv_id,
  cnv_type = type
)
evidence <- derive_cnv_2019_evidence(
  cnv, example$genes, example$dosage
)
evidence[, c(
  "criterion", "score", "source_id",
  "source_release", "evidence_id"
)]
#>   criterion score                       source_id source_release
#> 1        1A     0                         GENCODE            v19
#> 2        3A     0                         GENCODE            v19
#> 3        2A     1 ClinGen Gene Dosage Sensitivity     2026-07-03
#>               evidence_id
#> 1 genes:ENSG00000181449.2
#> 2            gene-count:1
#> 3       dosage:HGNC:11195
```

The rows can be scored directly:

``` r
score_cnv_2019(evidence)
#>          cnv_id cnv_type score classification criterion_groups
#> 1 clinvar-72461     loss     1     pathogenic                3
```

Here criterion 2A supplies one point because the deletion completely
contains an established haploinsufficient gene. The automatic provider
does not invent phenotype, inheritance, segregation, population,
case-report, or partial-gene evidence. Those require explicit additional
evidence rows.

## Published model without XGBoost

The published classifier is one tree with 41 nodes. XCNV stores the tree
as a small, validated tabular object and evaluates it directly:

``` r
model <- xcnv_bundled_model()
data.frame(
  trees = length(unique(model$nodes$tree_id)),
  nodes = nrow(model$nodes),
  features = length(model$feature_names)
)
#>   trees nodes features
#> 1     1    41       30

model_data <- new.env(parent = emptyenv())
data("xcnv_portable_model", package = "XCNV", envir = model_data)
identical(model_data$xcnv_portable_model, model)
#> [1] TRUE
```

`tools/stage_xcnv_model.R` rebuilds this portable form from the original
serialized booster and compares both evaluators on 1,024 deterministic
feature rows. XGBoost is needed for that maintainer-only conversion, not
for package installation or prediction.

## Full annotation resources

The bundled real example is deliberately valid only for ClinVar allele
72461. A general run uses a full GRCh37 resource directory:

``` r
resources <- xcnv_resources("/path/to/xcnv-resources", require = "annotations")
cnvs <- read_cnv("input.bed")
predictions <- predict_cnv(cnvs, resources)
write_xcnv(predictions, "predictions.output.csv")
```

`xcnv_resource_manifest()` lists the required files. `run_xcnv()`
provides the same workflow as one call and writes only when an output
path is supplied.

The compatibility contract, coordinate conventions, and current gaps are
in [`docs/upstream-compatibility.md`](docs/upstream-compatibility.md).

## Citation

Please cite Zhang et al. (2021), *Genome Medicine* 13:132,
[doi:10.1186/s13073-021-00945-4](https://doi.org/10.1186/s13073-021-00945-4).
