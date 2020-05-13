# Release 2.7.0 [2020-05-13]

### Changes
  * The value None is no longer supported for VCFRecord.filter in `anacore.vcf`.
  The field takes a list in three possible states:
    * If no filter was applied, the field contains an empty list ("." in VCF file)
    * If filters were applied but the record passes filters, the field should
    contain ["PASS"]
    * If filters were applied and the record does not pass filters, the field
    should contain ["filter_name", ...]

### Improvements
  * Add a classes to manage fusions detected by Arriba and STAR-Fusion in `anacore.fusion`.
  * Add classes to manage VCF containing fusions in `anacore.fusion`.

### Bug fixes
  * Fix bug in iterOverlappedByRegion from `anacore.region` when all chromosomes
  of queries are not in subjects.

# Release 2.6.0 [2020-01-03]

### Improvements
  * Add a Mutalyzer Batch manager in `anacore.hgvs`.
  * Add management for RTAComplete.txt V2 in `anacore.illumina.RTAComplete`.
  * Add a conda recipe.

### Bug fixes
  * Add automatic management for multiple date formats in `anacore.illumina.RunInfo`.

# Release 2.5.0 [2019-10-25]
First public release

# Release 1.0.0 [2017‑10‑07]
First release
