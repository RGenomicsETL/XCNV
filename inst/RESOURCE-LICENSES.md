# External resource and tool inventory

The following items are required for a scientific X-CNV run but are not package
payload. Their provenance and redistribution terms must be established before
release:

* `hg19_ljb26_all_converted_sites.vcf` and
  `hg19_ljb26_all_converted_scores.txt` (LJB26-derived annotations);
* `CDTS_percentile.txt` (CDTS annotation);
* `merged.cnv.sites.bed`, `merged.cnv.sample.info.txt`, `sample.info.txt`, and
  the training data referenced by the upstream repository;
* `hg19-ccREs.bed` (ENCODE candidate cis-regulatory elements);
* `gencode_v19_features.bed` (GENCODE v19-derived features);
* `xcnv.model.Rdata` (serialized XGBoost model) and any derived portable model;
* current ClinGen dosage-sensitivity, GENCODE/RefSeq gene, and population-SV
  releases staged for the independent ACMG/ClinGen 2019 scoring lane.

The historical `tools/bedtools2.tar.gz` archive has been removed. The package
does not compile, install, or execute bedtools. The synthetic files under
`inst/extdata/fixture/` are test-only and do not carry the scientific resource
claims above.
