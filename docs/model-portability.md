# Portable X-CNV model validation

Status: current validation record for the bundled published model.

The historical `tools/upstream/data/xcnv.model.Rdata` at XCNV commit
`2237e765fa9c31555dd8f4591b7f2693cff1ff8a` contains an XGBoost 1.2.0
`reg:logistic` classifier with 30 named features, one tree, and 41 nodes.
Because old XGBoost serialized objects are not stable interchange files,
maintainer conversion uses the matching R package, XGBoost 1.2.0.1.

`tools/stage_xcnv_model.R` converts the booster to the reviewable
`xcnv-tree-v1` table and compares predictions from both evaluators over 1,024
deterministic rows spanning all features, missing values, and the all-zero and
all-one cases. The observed maximum absolute difference is
`2.6508216177667521e-08`; the acceptance threshold is `1e-6`.

`tools/xcnv.portable.model.tsv` is the readable model source.
`tools/build_bundled_model_data.R` validates it and generates
`data/xcnv_portable_model.rda`, the model object exported by the installed
package. Package tests assert the 41-node structure, missing-value routing,
real ClinVar prediction, and identity of the public data object and default
model.
