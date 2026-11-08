# XCNV architecture authority

This is the current boundary and output contract for XCNV. The package is
relation-first: public computation returns ordinary tables, not one object per
CNV and not a VariantStory-specific wrapper.

## Ownership

- **XCNV** owns table-oriented X-CNV annotation/model output and the separate
  ACMG/ClinGen 2019 constitutional-CNV evidence and classification lane.
- **VariantStory** may consume these tables. It does not wrap, reimplement, or
  reinterpret either XCNV model behavior or the ACMG/ClinGen policy.
- **VariantStoryBench** generates and evaluates cases and records. It is not a
  source of XCNV policy or a replacement classifier.

## Relational contracts

X-CNV prediction preserves the legacy feature columns and `MVP_score`, then
adds `provider_id`, `assembly`, `resource_release`, `model_id`, `policy_id`,
and `status`. The published compatibility lane is identified as XCNV's
GRCh37 resource/model contract (`upstream-2237e765fa9c`) and has no clinical
policy id. `MVP_score` is a model result; it is never a clinical class.

ACMG/ClinGen evidence rows retain `provider_id`, `source_id`,
`source_release`, `evidence_id`, and optional evidence-model `model_id`; required
identifier, provenance, and contract labels are non-missing, nonempty
character/factor labels, never logical, numeric, list, or matrix values. The
automatic content provider uses the gene or dosage `source_id` as that evidence
row's provider identity. Group and summary rows carry a scalar evidence identity
only when it is unique, plus counts; otherwise that scalar is `NA`. A scalar
`model_id` is based only on selected rows that contribute to that group or
CNV/type score, and is emitted only when every such row has the same nonmissing
model; unselected rows cannot supply model attribution. Conflicting nonmissing
selected model identities for one CNV/type are rejected rather than silently
erased. Evidence provider identity is provenance, not a voting input, and is
never relabelled. Classification contracts are
resolved independently for each `(cnv_id, cnv_type)`; the summary `provider_id`
is always the XCNV ACMG/ClinGen policy kernel, never a provider conflict. An
omitted evidence `policy_id` defaults to
`acmg-clingen-2019-constitutional-cnv`; when supplied, every row must carry
that exact nonempty identity. Foreign or mixed policy identities are rejected
before scoring, rather than being rewritten into a native-policy summary. A
scalar `resource_release` is emitted only when one release applies, without
semicolon-concatenating unrelated provenance. Classification output keeps the
policy `score` and `classification` separate from X-CNV's `MVP_score`.
`status` is explicit: `supported`, `unspecified`, `unsupported`, or `conflict`;
only `supported` emits a clinical classification.

## Assembly and resources

- The published X-CNV compatibility model is calibrated against GRCh37/hg19
  resources. `xcnv_resources(assembly = "GRCh38")` is marked unsupported and
  prediction refuses it; XCNV does not claim GRCh38 compatibility.
- CNV calls on another assembly are scored only through an explicit
  `xcnv_liftover()` declaration: the caller stages the chain and the source and
  destination FASTA files and declares both assemblies and the chain release.
  XCNV never downloads, bundles, or infers a liftover, and the destination
  assembly must equal the resource assembly. DuckHTS `bcftools_liftover()`
  lifts each CNV on XCNV's DuckDB connection as one symbolic `<DEL>` (loss) or
  `<DUP>` (gain) record at POS `start + 1` with `END = end`; the lifted record
  is read back as the BED interval `[POS - 1, END)`. Rows keep their input
  coordinates and `assembly`, and add `lifted_assembly`, `lifted_chr`,
  `lifted_start`, `lifted_end`, `liftover_chain`, `liftover_release`,
  `liftover_length_ratio`, and `liftover_reason`. A call whose start or end
  does not lift, whose end lifts before its start, or that lifts to another
  contig or the reverse strand is `unsupported`, names that reason, and has no
  features or `MVP_score`. The length ratio is reported, not gated.
- The native ACMG/ClinGen lane accepts coherent CNV, gene, and dosage relations;
  nonempty relations declare one assembly. Before any overlap, all rows must have
  nonempty contigs and integer BED coordinates with `start < end` in the
  exactly representable range `[0, 2^53)`; character coordinates are
  range-checked before conversion. GRCh38 is supported when all nonempty
  relations declare GRCh38; mixed assemblies are rejected and other assemblies
  remain explicitly unsupported. A zero-row CNV relation is neutral and returns
  typed empty evidence after validating nonempty resources; a nonempty CNV
  relation requires nonempty gene and dosage resources, so an empty resource
  cannot fabricate negative evidence or provenance. Contigs are non-missing,
  nonempty character/factor labels without surrounding whitespace. Gene
  `biotype` values are explicit, nonempty labels, and the
  `functional_biotypes` vocabulary is additionally unique; both reject
  surrounding whitespace and placeholder tokens before functional filtering.
  Explicit nonfunctional Ensembl-style biotypes remain valid, while incomplete
  annotation cannot establish 1B/3A absence evidence.
- Native ACMG/ClinGen contig matching is exact. XCNV does not infer aliases
  such as `chr1`/`1`, and a valid CNV contig with no matching gene or dosage
  row is an admissible no-overlap result, not a resource-validation error.
- Source/resource release identity is carried in the tables. Large resources
  are staged by the user; XCNV does not download or bundle them.

## Evidence boundary

DECIPHER is an optional permitted evidence or case source only when the user
has the required access and redistribution controls. XCNV does not download,
bundle, or infer protected DECIPHER case evidence. Case, phenotype,
segregation, inheritance, and similar evidence remain explicit provider rows;
the automatic ACMG provider only derives the documented objective content
subset.
