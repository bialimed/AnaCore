## Description

********************************************************************************

## Intallation
### 1. Download code

********************************************************************************

The application folder has the following structure:

    <APP_DIR>/
    ├── bin/              # Scripts and softwares called used ADIVaR
    ├── jflow/            # Scripts and libraries for workflow management
    │   ├── bin/          # Scripts workflow management (run, rerun, monitor)
    │   ├── ...
    │   └── workflows/    # Workflow code
    ├── lib/              # Internal librairies
    ├── README.md
    └── ressources/       # Pre-loaded panels files


### 2. Install dependencies
The following programs must be accessible in PATH or added in softwares section
in jflow/application.properties.

* BWA
 Version: >=0.7.10
 Named as: bwa
 Download: http://bio-bwa.sourceforge.net/
* cutadapt
 Version: >=1.13
 Named as: cutadapt
 Download: http://cutadapt.readthedocs.io/en/stable/installation.html
* FastQC
 Version: >=0.11.3
 Named as: fastqc
 Download: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* InterOp
 Version: >=1.0.18
 Named as: dumptext
 Download: https://github.com/Illumina/interop
* jFlow
 Version: >=1.3
 Download: https://mulcyber.toulouse.inra.fr/scm/?group_id=186
* Python interpreter
 Version: >=3.5
 Additional libraries: numpy and pysam
* SAMtools
 Version: >=1.3.1
 Named as: samtools
 Download: http://www.htslib.org/
* VarDictJava
 Version: >=1.5.1
 Named as: bwa
 Download: https://github.com/AstraZeneca-NGS/VarDictJava
* VEP
 Version: >=89
 Named as: vep
 Download: http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

### 3. Test install
Use the folowing command to launch ADIVaR on a reduced dataset.

    <APP_DIR>/jflow/bin/jflow_cli.py adivar \
     @workflows/ADIVaR/data/test.cfg \
     --output-dir <OUT_DIR>

The file `<OUT_DIR>/variants/17T033348_filtered.tsv` must contain only one
variant on KRAS: C becomes A at position 25398285 on chromosome 12.

## Panel preparation
This step is only required once by panel. It prepares configuration files used
in workflow for a particular panel.

    # Library A
    manifestToBED.py \
      --without-primers \
      --input-manifest libA_IlluminaManifest.txt \
      --input-genome hg19.fasta \
      --output-BED libA_woutPrimers.bed
    manifestToBED.py \
      --input-manifest libA_IlluminaManifest.txt \
      --input-genome hg19.fasta \
      --output-BED libA_withPrimers.bed
    nonOverlappingDesign.py \
      --margin 10 \
      --input-panel libA_withPrimers.bed \
      --output-design libA_nonOvelapGroups.tsv
    # Library B
    manifestToBED.py \
      --without-primers \
      --input-manifest libB_IlluminaManifest.txt \
      --input-genome hg19.fasta \
      --output-BED libB_woutPrimers.bed
    manifestToBED.py \
      --input-manifest libB_IlluminaManifest.txt \
      --input-genome hg19.fasta \
      --output-BED libB_withPrimers.bed
    nonOverlappingDesign.py \
      --margin 10 \
      --input-panel libB_withPrimers.bed \
      --output-design libB_nonOvelapGroups.tsv

For INCa-V1 panel and INCa-V2 the files has already be processed and are
available in `<APP_DIR>/ressources/amplicon_design/INCa_V1/<ASSEMBLY>` or
`<APP_DIR>/ressources/amplicon_design/INCa_V2/<ASSEMBLY>`.

## Workflow management
### Launch
#### Input data

A `SampleSheet.csv` must be found in the folder where fastq for R1 and R2 are located.
This is the typically structure presents after sequencing in folder
`<RUN_FOLDER>/Data/Intensities/BaseCalls`.
The SampleSheet.csv has the following format:

    [Header]
    IEMFileVersion,4
    Investigator Name,IV1
    Experiment Name,INCa-V1_20170719CB
    Date,19/07/2017
    Workflow,Amplicon - DS
    Application,Amplicon - DS
    Assay,TruSight Tumor
    Description,
    Chemistry,Amplicon

    [Manifests]
    A,INCa_V1_A.txt
    B,INCa_V1_B.txt

    [Reads]
    151
    151

    [Settings]

    [Data]
    Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Manifest,GenomeFolder,Sample_Project,Description
    HORIZON_A,HORIZON,20170718_INCaV1CB,B01,A712,AGGAGTGG,A502,TGCTAAGT,A,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFasta,,
    HORIZON_B,HORIZON,20170718_INCaV1CB,B03,A702,ACAGTGGT,A502,TGCTAAGT,B,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFasta,,
    17T033348_A,17T033348,20170718_INCaV1CB,C01,A712,AGGAGTGG,A503,TGTTCTCT,A,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFasta,,
    17T033348_B,17T033348,20170718_INCaV1CB,C03,A702,ACAGTGGT,A503,TGTTCTCT,B,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFasta,,

#### Command line
The following command is the example used in installation test:

    <APP_DIR>/jflow/bin/jflow_bin.py adivar \
      --samplesheet <APP_DIR>/jflow/workflows/ADIVaR/data/samples/SampleSheet.csv \
      --filters <APP_DIR>/jflow/workflows/ADIVaR/data/ADIVAR_filters.json \
      # Reference
      --assembly-version GRCh37 \
      --genome-seq <APP_DIR>/jflow/workflows/ADIVaR/data/reference/Homo_sapiens_GRCh37_subset.fa \
      # Panel
      --RNA-selection <APP_DIR>/ressources/amplicon_design/INCa_V1/reference_RNA.tsv \
      --libA-folder <APP_DIR>/ressources/amplicon_design/INCa_V1/GRCh37/libA \
      --libB-folder <APP_DIR>/ressources/amplicon_design/INCa_V1/GRCh37/libB \
      # Protocol
      --R1-end-adapter <APP_DIR>/ressources/amplicon_design/adapters/Illumina_3prim_adapter.fasta \
      --R2-end-adapter <APP_DIR>/ressources/amplicon_design/adapters/Illumina_5prim_adapter_rvc.fasta \
      # Control samples
      --pos-ctrl-names HORIZON \
      --pos-ctrl-expected  <APP_DIR>/jflow/workflows/ADIVaR/data/panel_on_GRCh37/pos_ctrl_expected.vcf

Use `<APP_DIR>/jflow/bin/jflow_bin.py adivar --help` for more information about
parameters.

### Monitor
#### For monitoring all workflows

Use the following command:

    <APP_DIR>/jflow/bin/jflow_admin.py status

This command has the following output:

    ID      NAME    STATUS  ELAPSED_TIME    START_TIME      END_TIME
    000002  adivar  completed       0:14:30 Tue Sep 26 13:55:15 2017        Tue Sep 26 14:09:46 2017
    000001  adivar  completed       0:11:02 Mon Sep 25 13:39:31 2017        Mon Sep 25 13:50:34 2017

In this example two workflows has been processed and completed without errors.

#### For monitoring a specific workflow

Use the following command:

    <APP_DIR>/jflow/bin/jflow_manager.py status \
      --workflow-id <YOUR_WF_ID> \
      --errors

This command has the following output:

    ********************* example de sortie ************************************

### Rerun
You can rerun failed/incomplete steps with the following command:

    <APP_DIR>/jflow/bin/jflow_manager.py rerun \
      --workflow-id <YOUR_WF_ID>

### Troubleshooting
#### Triallelic variants

********************************************************************************

#### Deletions with same homopolymer upstream and downstream the variant

********************************************************************************

## License
GNU GPL v3

## Copyright
2017 Laboratoire d'Anatomo-Cytopathologie de l'Institut Universitaire du Cancer
Toulouse - Oncopole

## Contact
escudie.frederic@iuct-oncopole.fr
