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
and never lets an X-CNV model score select or admit evidence. Evidence and
group rows retain contributing provider/source identities and linked evidence
counts; provider identity is provenance, not a voting or rewritten field.
`derive_cnv_2019_evidence()` uses the gene or dosage `source_id` as that row's
evidence `provider_id`. Every required evidence identifier, provenance, and
supplied contract value is a non-missing, nonempty character/factor label;
logical, numeric, list, and matrix values fail before conversion. Score labels
are parsed through their character labels, never factor codes. A scalar evidence
`model_id` is taken only from selected rows that contribute to its group or
CNV/type score, and is emitted only when every such row has the same nonmissing
model; unselected rows cannot supply it. Conflicting nonmissing selected model
identities for one CNV/type are rejected. An omitted evidence
`policy_id` defaults to `acmg-clingen-2019-constitutional-cnv`; when supplied,
it must be that exact, nonempty identity on every row. Foreign or mixed policies
are rejected before scoring rather than rewritten in a native-policy summary. Classification
contracts are resolved per `(cnv_id, cnv_type)`; the summary provider is the
XCNV ACMG/ClinGen policy kernel. A clinical class is emitted only for
`supported` status. Assembly-unspecified, unsupported, or conflicting
assembly/status contracts retain scores and provenance but return `NA`
classification, and multiple releases are not concatenated into one misleading
scalar.

`derive_cnv_2019_evidence()` implements the current automatic subset: section 1
functional content, section 2A complete containment of an established ClinGen
HI/TS feature, and section 3 protein-coding gene count. Its CNV, gene, and
dosage inputs must have nonempty contigs and exactly representable nonnegative
integer half-open BED coordinates in `[0, 2^53)` with `start < end`; character
coordinates are range-checked before conversion, so values at or beyond `2^53`
fail before the DuckDB interval join. A zero-row CNV relation is neutral after
nonempty resources agree, while nonempty CNVs require nonempty gene and dosage
relations with explicit provenance; an empty resource is never treated as
negative evidence. Contig identifiers are non-missing, nonempty
character/factor labels without surrounding whitespace. Gene `biotype` values
are explicit, nonempty character/factor labels, and `functional_biotypes` is an
explicit, nonempty, unique character/factor vocabulary. Both reject surrounding
whitespace or placeholder tokens before functional filtering. Gene biotype
values need not be unique across rows. Legitimate explicit nonfunctional
Ensembl-style biotypes remain valid, while incomplete annotation cannot
establish 1B/3A absence evidence. Contig identifiers are exact caller/provider
contracts: XCNV
does not infer aliases, and a valid contig with no matching resource row is a
no-overlap result. This provider is backed by the same DuckDB
interval join as X-CNV annotation. Partial-gene mechanisms,
PVS1/NMD, phenotype specificity, published cases, segregation, inheritance,
and population evidence are deliberately outside this automatic subset until a
separately tested provider supplies them.

Current full annotation packs are external and release-labelled. As of July 2026,
the source-authority manifest identifies ClinGen's nightly dosage sensitivity
downloads, GENCODE 50 for GRCh38/Ensembl 116, gnomAD SV 4.1, assembly-specific
NCBI RefSeq, and dated DGV snapshots. The manifest intentionally points to
official landing pages; staging requires an explicit artifact URL, release, and
assembly.
