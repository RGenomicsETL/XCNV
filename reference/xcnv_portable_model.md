# Published X-CNV portable model

The fitted X-CNV classifier from Zhang et al. (2021), converted once
from its serialized XGBoost representation to the package's
runtime-neutral `xcnv-tree-v1` relation.

## Format

An `xcnv_portable_model` list containing the logistic base score, 30
feature names, and a 41-row tree relation.

## Source

Zhang et al. (2021),
[doi:10.1186/s13073-021-00945-4](https://doi.org/10.1186/s13073-021-00945-4)
.
