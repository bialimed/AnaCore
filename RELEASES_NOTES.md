# Release 2.12.1 [2022-10-01]

### Improvements
  * Add MiSeq invalid RFID markup management in `anacore.illumina.RunParameters`.

### Bug fixes
  * Fix bug in `anacore.msi.sample.setStatusByInstabilityRatio()`: the value
  voting_loci was invalid.
  * Fix bug in `anacore.msi.sample.setScore()`: the value locus_weight_is_score
  was ignored.

# Release 2.12.0 [2022-09-06]

### Changes
  * `anacore.vcf.VCFIO`: Manage None value in a vcf INFO field as missing key.
  Example: for an INFO field equal to `AF=0.5;DP=.` the record.info is
  `{"AF": 0.5}`.
  * `anacore.msi.msings`:
    * Remove `record.results["mSINGS"].data["peaks"]` from returned records in parser
    `anacore.msi.msings.MSINGSAnalysisIO`: data can be retrieved from
    `record.results["mSINGS"].data["lengths"]`.
    * Rename `MSINGSAnalysis` to `MSINGSAnalysisIO`.
    * Rename method name `MSINGS` to `mSINGS`.
  * Change parameter behaviour for `min_voting_loci` in
  `anacore.msi.sample.MSISample.setStatusByInstabilityRatio()`: from number of
  loci to rate of loci.
  * Move `anacore.msi.MSIReport` to `anacore.msi.reportIO.ReportIO`.
  * Refactor `anacore.msi.LocusRes*` to create `anacore.msi.locus.LocusDataDistrib`.
  This class store length distribution and is linked to `anacore.msi.locus.LocusRes`
  in `data["lengths"]`. The class `anacore.msi.LocusResDistrib` and children are
  removed.
  * Move MSI libraries in the new subpackage `msi`:
    * `anacore.msi` to `anacore.msi.base`
    * `anacore.msiannot` to `anacore.msi.annot`
    * `anacore.msings` to `anacore.msi.msings`

### Improvements
  * Add `anacore.illumina.DemultStat` to read demultiplex statistics from
  bcl2fastq.
  * Add `anacore.illumina.Run` getters to known run information and status from
  the run folder.
  * Add `anacore.msi.hubble` to manage results from Hubble software.
  * Add `anacore.msi.msisensorpro` to manage results from MSIsensor-pro software.

# Release 2.11.0 [2022-03-10]

### Improvements
  * Update pysam from `0.15.3` to `0.18.0`.
  * Update numpy from `1.6.0` (pypi) or `1.16.5` (conda) to `1.20.1`.

# Release 2.10.0 [2021-09-28]

### Improvements
  * Add `anacore.illumina.Bcl2fastqLog` to read bcl2fastq log file.
  * Add alias "Description" for "Sample_Description" in `samplesheet.samples`
  (`anacore.illumina.SampleSheetIO`).
  * Add empty value for `samplesheet.header["Description"]` when "Description"
  is not present in SampleSheet (`anacore.illumina.SampleSheetIO`).
  * Add utilities to manage Homo sapiens genome accessions in `anacore.db.homo_sapiens.accession`.

### Bug fixes
  * Fix bug in `anacore.annotVcf.AnnotVCFIO` when parsing ANN declaration from
  SnpEff.
  * Fix no casting itemRGB as list in `anacore.bed.BEDIO`.

# Release 2.9.0 [2020-04-28]

### Changes
  * `anacore.genomicRegion.Protein.setTranscript` and `anacore.genomicRegion.Transcript.setProteins`
  replaced by setter declaration.
  * Move `anacore.sequenceIO.Sequence` to `anacore.sequence.Sequence`.

### Improvements
  * Implement `getPosOnRegion()` in `anacore.genomicRegion.Transcript`.
  * Add functions to get information about codon in `anacore.genomicRegion.Protein`:
  `getCodonRefPos()`, `getCodonSeqFromProtPos()` and `getCodonInfo()`.
  * Add `AA3LettersAlphabet`, `CodonAlphabet`, `DNAAlphabet` and `RNAAlphabet` in
  `anacore.sequence` to validate sequences and provide translation and reverse
  complement utilities.
  * Add `anacore.vcf.VCFRecord.fastDownstreamed` to get quickly the most downstream
  version of the variant.
  * Add management of metadata in SV files (`anacore.sv`). Metadata must be
  present before title and/or data. They starts with a particular string: "##"
  by default.
  * Add `anacore.hgvs.HGVSProtChange` to manage change part of proteic HGVS (ex:
  "Val600Glu").
  * Add `getSub` and open mode "i" in `anacore.vcf.VCFIO` to return records
  overlapping the specified region in file with tabix index.
  * Add specific management for samples description rows in VCF header
  (`anacore.vcf.VCFIO`).
  * Add reader for picard tools outputs in `anacore.picardIO`.
  * Manage UMI in Illumina's sequence ID with `anacore.illumina.getInfFromSeqID`.

# Release 2.8.0 [2020-10-09]

### Improvements
  * Add management for None value in `anacore.filters.Filter`.
  * Increase speed to read the VCF in `anacore.vcf.VCFIO`.

### Bug fixes
  * Fix bug in GTFIO from `anacore.gtf` when an attribute value contains
  semicolon.
  * Fix bug with empty list in an INFO field from VCF (`anacore.vcf`).
  Previously, the reader returned a list containing an empty string. For example,
  for the INFO field containing "AF=0.5;DB=;DP=100" where DB is a list, the reader
  returned: {..., "DB": [""], ...}. Now, the reader return: {..., "DB": [],
  ...}.

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
