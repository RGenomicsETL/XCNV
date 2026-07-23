# Upstream compatibility and design note

## Authority and scope

This package is being developed against the Zhang et al. (2021) X-CNV
executable in the upstream repository `https://github.com/kbvstmd/XCNV` at
pinned commit `2237e765fa9c31555dd8f4591b7f2693cff1ff8a` (19 October 2021).
The historical files under `tools/upstream/data/` remain a reviewable
compatibility reference and are excluded from built packages. They are kept
outside the package `data/` directory so a source install does not mistake raw
annotations or legacy serializations for R datasets.

The supported subset in this first vertical slice is the published four-column
BED-like input (`chromosome`, `start`, `end`, `gain`/`loss`), GRCh37/hg19
resources, all five annotation families used by the executable (LJB26, CDTS,
merged CNV frequency, ENCODE ccRE, and GENCODE haploinsufficiency), and MVP
prediction with the published classifier. The package evaluates an equivalent
runtime-neutral tree rather than loading XGBoost.

## Preserved semantics

* Chromosomes are the strings `1` through `22`, `X`, and `Y`; the executable's
  first-row-triggered literal `chr` removal is retained.
* Input records use the first four tab-separated fields. Records are sorted by
  chromosome order, start, and end before annotation.
* The BED overlap arithmetic is half-open. Frequency annotations retain both
  reciprocal overlap thresholds `-f 0.5 -F 0.5`.
* CNV length and regulatory/coverage denominators retain `end - start + 1`,
  while overlap lengths retain BED `end - start` arithmetic.
* LJB26 missing values use `0, 0, -12.30, -11.958, -20.000, 0.0003`; absent
  LJB26 positive calls contribute to the original fraction denominator. The
  converted sites file's fourth field is the score-table row index, matching
  the upstream `awk $5,$9` selection.
* CDTS missing values are 100, each qualifying record contributes 10 bases, and
  the two coverage scores are capped at 1.
* ENCODE labels, GENCODE max aggregation, feature ordering, model feature names,
  gain coding (`Type.1 = 1`, loss coding `0`), and output ordering are retained.
  The R result uses the unique name `Type.1`; the CSV compatibility wrapper
  emits the upstream duplicate `Type` header.
* The historical single-row duplicate/50-base workaround and sorted-row removal
  are retained, including its short-interval edge behavior.

The R API is table-oriented: `read_cnv()`, `annotate_cnv()`, `predict_cnv()`,
and `run_xcnv()`/`write_xcnv()`. It does not create one object per CNV.

## Deliberate implementation choices

The package uses DuckDB inequality joins. The historical 2016 bedtools archive
has been removed. A small ordinary-R overlap implementation is the independent
test oracle. The upstream reference scripts still show their original shell
pipeline, but those scripts are excluded from the package and are not called by
this API.

The installed package does not depend on XGBoost. The maintainer staging tool
`tools/stage_xcnv_model.R` converts the serialized booster to the stable
`xcnv-tree-v1` TSV contract and differentially validates 1,024 probe rows
against the source booster. The package bundles and evaluates the resulting
one-tree, 41-node object.

No full-genome annotation or training resource is bundled. A user places
full-genome resources in a directory and calls `xcnv_resources(path)`. A
compact ClinVar 72461 resource subset is bundled solely for the evaluated
example and regression test.

## Parity limits

The real ClinVar 72461 example now executes LJB26, CDTS, population-frequency,
ENCODE, GENCODE, and model evaluation together. A whole-genome differential
against the historical executable remains outstanding. The output header count
in the old README also does not match the number of population columns produced
by `sample.info.txt`; the implementation follows the executable's `PAF`
construction rather than hard-coding the README's count.

## Independent ACMG/ClinGen 2019 lane

The package's constitutional-CNV scoring lane is authored from Riggs et al.
(2020), DOI `10.1038/s41436-019-0686-8`. Its evidence relation requires source
identifiers, source releases, and evidence record identifiers. Numeric
criterion ranges and aggregate caps are machine readable through
`cnv_2019_criteria()`. The scorer keeps raw and capped group totals separately
and never lets an X-CNV model score select or admit evidence.

`derive_cnv_2019_evidence()` implements the current automatic subset: section 1
functional content, section 2A complete containment of an established ClinGen
HI/TS feature, and section 3 protein-coding gene count. This provider is backed
by the same DuckDB interval join as X-CNV annotation. Partial-gene mechanisms,
PVS1/NMD, phenotype specificity, published cases, segregation, inheritance,
and population evidence are deliberately outside this automatic subset until a
separately tested provider supplies them.

Current full annotation packs are external and release-labelled. As of July 2026,
the source-authority manifest identifies ClinGen's nightly dosage sensitivity
downloads, GENCODE 50 for GRCh38/Ensembl 116, gnomAD SV 4.1, assembly-specific
NCBI RefSeq, and dated DGV snapshots. The manifest intentionally points to
official landing pages; staging requires an explicit artifact URL, release, and
assembly.
