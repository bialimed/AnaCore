#!/bin/bash

mkdir /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37


# Lib A
echo "libA"

mkdir /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libA

manifestToBED.py \
 --input-manifest /save/fescudie/data/amplicon_design/Lymphome_B/TSCA_Manifest_117223.txt \
 --input-genome /work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa \
 --output-BED /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libA/withPrimers.bed

manifestToBED.py \
 --without-primers \
 --input-manifest /save/fescudie/data/amplicon_design/Lymphome_B/TSCA_Manifest_117223.txt \
 --input-genome /work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa \
 --output-BED /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libA/woutPrimers.bed

nonOverlappingDesign.py \
 --input-panel /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libA/withPrimers.bed \
 --output-design /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libA/nonOverlappingGroups.tsv

 
# Lib B
echo "libB"

mkdir /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libB

manifestToBED.py \
 --input-manifest /save/fescudie/data/amplicon_design/Lymphome_B/TSCA_Manifest_117224.txt \
 --input-genome /work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa \
 --output-BED /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libB/withPrimers.bed

manifestToBED.py \
 --without-primers \
 --input-manifest /save/fescudie/data/amplicon_design/Lymphome_B/TSCA_Manifest_117224.txt \
 --input-genome /work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa \
 --output-BED /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libB/woutPrimers.bed

nonOverlappingDesign.py \
 --input-panel /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libB/withPrimers.bed \
 --output-design /save/fescudie/data/amplicon_design/Lymphome_B/GRCh37/libB/nonOverlappingGroups.tsv
