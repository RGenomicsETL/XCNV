# XCNV 0.0.0.9000

* Corrected ACMG/ClinGen contract resolution to be per CNV. Evidence rows now
  retain their source provider identity (including gene/dosage providers), while
  summaries use only the policy-kernel provider. Group and summary provenance
  is scalar only when unique, with counts otherwise; clinical classes emit only
  for supported status, and unrelated releases are not concatenated.

* Added an explicit relational contract for provider, assembly, resource
  release, model/policy identity, and supported/unsupported/conflict status.
  The X-CNV compatibility lane remains GRCh37-calibrated; the ACMG/ClinGen
  2019 lane validates coherent declared assemblies and keeps model scores
  separate from clinical classes. DECIPHER is documented only as an optional,
  access-controlled evidence source and is not downloaded or bundled.

* Bundled the published 41-node X-CNV classifier as a validated
  runtime-neutral R data object, and made it the default when an annotation
  resource directory does not supply an override. Runtime prediction remains
  independent of XGBoost.
* Bundled a compact real ClinVar 72461/SOX2 example with its case-local X-CNV
  annotations and objective GENCODE/ClinGen content evidence.
* Corrected the published-resource readers: LJB26 `START == END` records are
  one-base points rather than empty BED intervals, LJB scores are selected by
  their explicit site index, and the headerless merged-sample file retains its
  first observation.
* Added an installable R package API for reading CNV tables, validating resource
  bundles, computing X-CNV annotations, predicting MVP scores, and writing the
  legacy table format.
* Added a small deterministic fixture and tinytest coverage, including an
  independent R overlap oracle for the production DuckDB inequality join.
* Removed bedtools and XGBoost from the installed runtime. A deterministic
  staging script converts a separately obtained XGBoost model to a validated,
  tabular `xcnv-tree-v1` artifact; the package evaluates that artifact in R.
* Added an independently authored ACMG/ClinGen 2019 constitutional-CNV scoring
  contract with provenance-required evidence rows, published score ranges,
  criterion-group caps, and five-tier classification. A DuckDB-backed content
  provider derives only sections 1, 2A, and 3 from explicit gene and ClinGen
  dosage relations; case, phenotype, inheritance, partial-gene, and population
  evidence remain curator/provider inputs.
* Recorded the pinned upstream compatibility scope and current resource
  limitations.
