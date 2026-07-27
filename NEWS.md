# XCNV 0.0.0.9000

* Selected DuckDB's explicit non-shared home mode when supported, keeping
  internal overlap connections silent and README rendering deterministic across
  DuckDB releases.

* Tightened ACMG/ClinGen relation labels: contigs now require non-missing,
  nonempty character/factor labels without surrounding whitespace. Gene biotypes
  reject surrounding whitespace and case-insensitive placeholder tokens, and
  `functional_biotypes` must be a nonempty unique explicit vocabulary. Legitimate
  explicit nonfunctional Ensembl-style biotypes remain valid and produce the
  contract's 1B/3A absence criteria when appropriate.

* Hardened the native ACMG/ClinGen scorer contract: an omitted `policy_id`
  defaults to the native policy, but a supplied value must be nonempty and
  exactly `acmg-clingen-2019-constitutional-cnv` on every row. Foreign or
  mixed policies now fail before a score or classification is emitted; evidence
  provider identity remains provenance rather than a voting or rewritten field.

* Hardened native ACMG/ClinGen input evidence and BED geometry: required
  provenance/identifiers and supplied contract labels now require character or
  factor labels, rejecting logical, numeric, list, and matrix values before
  conversion; score labels are parsed without factor-code coercion. Model
  attribution now uses only selected rows contributing to each group or summary
  score, emits a scalar only for a complete common selected model, and never
  broadcasts an unselected model. CNV, gene, and dosage coordinates must be
  exactly representable integers in `[0, 2^53)` before overlap; strings
  beyond that range fail before conversion. A zero-row CNV relation is neutral,
  while nonempty CNVs require nonempty, provenance-labelled gene and dosage
  resources and therefore cannot fabricate negative evidence. Missing, empty,
  logical, or list-valued gene biotypes and functional-biotype vocabularies now
  fail before filtering, so they cannot produce 1B/3A absence evidence. Contig names
  remain exact caller/provider identifiers, and a valid contig with no matching
  resource remains a no-overlap result rather than an error.

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
