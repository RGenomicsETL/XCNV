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
`source_release`, and `evidence_id`; the automatic content provider uses the
gene or dosage `source_id` as that evidence row's provider identity. Group and
summary rows carry a scalar
evidence identity only when it is unique, plus counts; otherwise that scalar is
`NA`. Classification contracts are resolved independently for each
`(cnv_id, cnv_type)`; the summary `provider_id` is always the XCNV
ACMG/ClinGen policy kernel, never a provider conflict. A scalar
`resource_release` is emitted only when one release applies, without
semicolon-concatenating unrelated provenance. Classification output keeps the policy `score` and `classification`
separate from X-CNV's `MVP_score`. `status` is explicit: `supported`,
`unspecified`, `unsupported`, or `conflict`; only `supported` emits a clinical
classification.

## Assembly and resources

- The published X-CNV compatibility model is calibrated against GRCh37/hg19
  resources. `xcnv_resources(assembly = "GRCh38")` is marked unsupported and
  prediction refuses it; XCNV does not silently lift coordinates or claim
  GRCh38 compatibility.
- The native ACMG/ClinGen lane accepts coherent CNV, gene, and dosage relations
  declaring one assembly. GRCh38 is supported when all three relations declare
  GRCh38; mixed assemblies are rejected and other assemblies remain explicitly
  unsupported.
- Source/resource release identity is carried in the tables. Large resources
  are staged by the user; XCNV does not download or bundle them.

## Evidence boundary

DECIPHER is an optional permitted evidence or case source only when the user
has the required access and redistribution controls. XCNV does not download,
bundle, or infer protected DECIPHER case evidence. Case, phenotype,
segregation, inheritance, and similar evidence remain explicit provider rows;
the automatic ACMG provider only derives the documented objective content
subset.
