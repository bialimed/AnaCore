## Description

## Intallation
### 1. Download code
### 2. Install dependencies
The following programs must be accessible in PATH or added in softwares section
in jflow/application.properties.

* BWA
 Version: >=*****************************
 Named as: bwa
 Download: http://bio-bwa.sourceforge.net/
* cutadapt
 Version: >=1.8******************************************
 Named as: cutadapt
 Download: http://cutadapt.readthedocs.io/en/stable/installation.html
* FastQC
 Version: >=******************************
 Named as: fastqc
 Download: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* InterOp
 Version: >=1.0.25************************************
 Named as: dumptext
 Download: https://github.com/Illumina/interop
* jFlow
 Version: >=1.3
 Download: https://mulcyber.toulouse.inra.fr/scm/?group_id=186
* Python interpreter
 Version: >=3.*****************************
 Additional libraries: numpy and pysam
* SAMtools
 Version: >=1.*****************************
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
available ********************************************** in folder or at the following address

## Workflow management
### Launch
    <JFLOW_DIR>/bin/jflow_bin.py adivar **********************************

### Monitor
* For monitoring all workflows
    <JFLOW_DIR>/bin/jflow_manager.py status
    ********************* example de sortie
* For monitoring a specific workflow
    <JFLOW_DIR>/bin/jflow_manager.py status --workflow-id <YOUR_WF_ID> --errors
    ********************* example de sortie

### Rerun
You can rerun failed/incomplete steps with the following command.
    <JFLOW_DIR>/bin/jflow_manager.py rerun --workflow-id <YOUR_WF_ID>
