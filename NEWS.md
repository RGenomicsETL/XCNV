# XCNV 0.0.0.9000

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
* Recorded the pinned upstream compatibility scope and the unresolved upstream,
  resource, model, and tool licensing blockers. ClassifyCNV is an external
  differential reference only because its license prohibits modification and
  redistribution; no ClassifyCNV code or databases are copied.
