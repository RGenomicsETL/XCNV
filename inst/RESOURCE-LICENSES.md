# External resource and tool inventory

Full-genome X-CNV annotation uses the following external scientific resources:

* `hg19_ljb26_all_converted_sites.vcf` and
  `hg19_ljb26_all_converted_scores.txt` (LJB26-derived annotations);
* `CDTS_percentile.txt` (CDTS annotation);
* `merged.cnv.sites.bed`, `merged.cnv.sample.info.txt`, `sample.info.txt`, and
  the training data referenced by the upstream repository;
* `hg19-ccREs.bed` (ENCODE candidate cis-regulatory elements);
* `gencode_v19_features.bed` (GENCODE v19-derived features);
* current ClinGen dosage-sensitivity, GENCODE/RefSeq gene, and population-SV
  releases staged for the independent ACMG/ClinGen 2019 scoring lane.

The historical `tools/bedtools2.tar.gz` archive has been removed. The package
does not compile, install, or execute bedtools. The package includes the
published 41-node classifier in a runtime-neutral tabular form. The files under
`inst/extdata/real_case/` are case-limited subsets for ClinVar allele 72461 and
must not be used to annotate unrelated CNVs.
