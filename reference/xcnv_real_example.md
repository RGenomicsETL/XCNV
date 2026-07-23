# Load the bundled real X-CNV example

Return ClinVar allele 72461, a pathogenic deletion containing `SOX2`,
with the case-limited GRCh37 annotation relations needed to run X-CNV
and the objective ACMG/ClinGen 2019 content-evidence lane. The
annotation files are subsets for this case and must not be used for
unrelated CNVs.

## Usage

``` r
xcnv_real_example()
```

## Value

A list containing `cnv`, `clinvar`, `resources`, `genes`, and `dosage`.
