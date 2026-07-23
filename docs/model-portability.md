# Portable X-CNV model validation

Status: current validation record for the unredistributed pinned upstream model.

The historical `data/xcnv.model.Rdata` at XCNV commit
`2237e765fa9c31555dd8f4591b7f2693cff1ff8a` has SHA-256
`76f64547a74ea478047dabc48751ec7e7d5e8fbeeac4e9e2a84941702d00ba8f`.
Its serialized configuration reports XGBoost 1.2.0, objective `reg:logistic`,
30 named features, one tree, and 41 nodes. Because XGBoost serialized objects
are not stable model files, conversion must use the matching R package
`xgboost` 1.2.0.1. Current XGBoost correctly rejects or warns about direct
deserialization of this object.

The `xcnv-tree-v1` evaluator was independently checked on 21 July 2026 against
an XGBoost compatibility load of that pinned serialized object. The probe used
4,096 rows over all 30 features, included all-zero and all-one rows plus missing
values, and observed a maximum absolute prediction difference of
`2.6508216177667521e-08`. The portable tree contained the same one tree and 41
nodes. The checked evaluator threshold is `1e-6`.

`tools/stage_xcnv_model.R` makes this differential mandatory for every staged
artifact, using 1,024 deterministic probe rows, and writes the observed maximum
difference into the sibling receipt. The package fixture separately tests
branch direction, missing-value routing, and logistic transformation without
XGBoost.

No portable derivative of the upstream scientific model is committed or
distributed. The upstream repository provides no license grant, and changing
serialization does not change the model's copyright or redistribution status.
