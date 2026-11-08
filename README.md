XCNV
================

<!-- README.md is rendered from this file. -->

XCNV annotates copy-number variants (CNVs) and scores them in two
separate lanes:

- **X-CNV**: the published X-CNV model (Zhang et al. 2021) gives each
  CNV an `MVP_score`, a model probability of pathogenicity. It is not a
  clinical class.
- **ACMG/ClinGen**: the 2019 constitutional CNV standard. Each criterion
  is one evidence row that you can inspect, change and rescore.

DuckDB does the overlap joins. R does validation, feature reduction, the
tree evaluator and scoring. No bedtools and no XGBoost are needed at run
time. Ownership, output contracts and validation rules are in
[`ARCHITECTURE.md`](ARCHITECTURE.md).

## Installation

``` r
install.packages(
  "XCNV",
  repos = c("https://rgenomicsetl.r-universe.dev", "https://cloud.r-project.org")
)
```

## Example: a ClinVar deletion

The package includes the annotations for one real case, ClinVar allele
72461: a pathogenic deletion that contains `SOX2`.

``` r
library(XCNV)
example <- xcnv_real_example()
example$clinvar[, c("clinvar_allele_id", "clinical_significance",
                    "grch37_chrom", "grch37_start0", "grch37_end0")]
#>   clinvar_allele_id clinical_significance grch37_chrom grch37_start0
#> 1             72461            Pathogenic            3     181428674
#>   grch37_end0
#> 1   181466445
```

### X-CNV score

``` r
prediction <- predict_cnv(example$cnv, example$resources)
prediction[, c("Chr", "Start", "End", "Type", "gain.freq", "loss.freq",
               "pLI", "Episcore", "GHIS", "MVP_score")]
#>   Chr     Start       End Type gain.freq loss.freq       pLI  Episcore
#> 1   3 181428674 181466445 loss         0         0 0.3966965 0.9209667
#>        GHIS MVP_score
#> 1 0.6169144 0.8537559
```

### ACMG/ClinGen evidence and class

The deletion contains `SOX2`, an established haploinsufficient gene, so
criterion 2A gives one point.

``` r
cnv <- transform(example$cnv, cnv_id = example$clinvar$cnv_id,
                 cnv_type = type, assembly = "GRCh37")
genes <- transform(example$genes, assembly = "GRCh37")
dosage <- transform(example$dosage, assembly = "GRCh37")
evidence <- derive_cnv_2019_evidence(cnv, genes, dosage)
evidence[, c("criterion", "score", "source_id", "evidence_id")]
#>   criterion score                       source_id             evidence_id
#> 1        1A     0                         GENCODE genes:ENSG00000181449.2
#> 2        3A     0                         GENCODE            gene-count:1
#> 3        2A     1 ClinGen Gene Dosage Sensitivity       dosage:HGNC:11195
score_cnv_2019(evidence)[, c("cnv_id", "cnv_type", "score", "classification", "status")]
#>          cnv_id cnv_type score classification    status
#> 1 clinvar-72461     loss     1     pathogenic supported
```

The automatic provider derives sections 1, 2A and 3 from gene and dosage
tables. Phenotype, inheritance, segregation, population, case and
partial-gene evidence are extra rows that you supply.

## GRCh38 calls

The X-CNV model and its resources are GRCh37. Declare a liftover to
score GRCh38 calls. XCNV lifts each call with DuckHTS and scores the
GRCh37 interval:

``` r
liftover <- xcnv_liftover(
  chain = "hg38ToHg19.over.chain",
  source_fasta = "GRCh38.fa", destination_fasta = "GRCh37.fa",
  source_assembly = "GRCh38", destination_assembly = "GRCh37",
  release = "ucsc-hg38ToHg19"
)
predict_cnv(grch38_calls, resources, liftover = liftover)
```

Each row keeps its GRCh38 coordinates and adds the lifted interval. A
call that does not lift to one interval on the same contig and strand is
`unsupported`, with a `liftover_reason` and no `MVP_score`. The
ACMG/ClinGen lane works on GRCh38 directly.

## Resources

The resources are on Hugging Face at
[`RgenomicsETL/xcnv-resources`](https://huggingface.co/datasets/RgenomicsETL/xcnv-resources).
`SHA256SUMS` covers every file.

| Folder | Contents | Used by |
|----|----|----|
| `grch37/` | X-CNV annotations and model | X-CNV lane |
| `grch38/` | GENCODE 50 genes and ClinGen dosage as Parquet, with `sources.tsv` | ACMG/ClinGen lane |

Score calls with the X-CNV lane:

``` r
resources <- xcnv_resources("xcnv-resources/grch37", require = "annotations")
predictions <- predict_cnv(read_cnv("input.bed"), resources)
write_xcnv(predictions, "predictions.output.csv")
```

`xcnv_resource_manifest()` lists the required files. `run_xcnv()` does
the same in one call.

Derive ACMG/ClinGen evidence for GRCh38 calls:

``` r
con <- DBI::dbConnect(duckdb::duckdb())
genes <- DBI::dbGetQuery(con, "FROM 'xcnv-resources/grch38/genes.parquet'")
dosage <- DBI::dbGetQuery(con, "FROM 'xcnv-resources/grch38/dosage.parquet'")
calls <- data.frame(
  cnv_id = "sox2-del", chr = "chr3", start = 181700000, end = 181800000,
  cnv_type = "loss", assembly = "GRCh38"
)
evidence <- derive_cnv_2019_evidence(calls, genes, dosage)
score_cnv_2019(evidence)
```

On these relations the call above scores 1 point (2A, `SOX2`):
pathogenic.

`tools/build_grch38_resources.R` rebuilds `grch38/` from SHA-256-pinned
sources.

## The model

The published classifier is one tree with 41 nodes. XCNV stores it as a
small table and evaluates it in R. `tools/stage_xcnv_model.R` rebuilds
this table from the original XGBoost booster and checks both evaluators
on 1,024 rows; only that maintainer step needs XGBoost.

## Citation

Zhang et al. (2021), *Genome Medicine* 13:132,
[doi:10.1186/s13073-021-00945-4](https://doi.org/10.1186/s13073-021-00945-4).
