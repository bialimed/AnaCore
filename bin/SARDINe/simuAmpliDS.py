#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import os
import sys
import time
import random
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from sequenceIO import *


########################################################################
#
# FUNCTIONS
#
########################################################################
def revcom( seq ):
    """
    @summary : Reverse complement the sequence.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}

    return( "".join([complement_rules[base] for base in seq[::-1]]) )

def getAmpliconsFromManifest( manifest_path ):
    amplicons = list()
    with open(manifest_path) as FH_manifest:
        section_probe = False
        probes_header = list()
        for line in FH_manifest:
            if line.strip() != "":
                if re.search('^\[\w+\]$', line.strip()) is not None:
                    if line.strip() != "[Probes]":
                        section_probe = False
                    else:
                        section_probe = True
                        probes_header = [field.strip().lower().replace(" ", "_") for field in FH_manifest.readline().split("\t")]
                elif section_probe:
                    fields = { probes_header[idx]:field.strip() for idx, field in enumerate(line.split("\t"))}
                    amplicons.append( fields )
    return amplicons

def getAmplicons( reference_path, manifest_path ):
    # Get amplicons
    amplicons = getAmpliconsFromManifest( manifest_path )
    amnplicons_by_chr = dict()
    for ampli in amplicons:
        chr = ampli["chromosome"]
        if chr not in amnplicons_by_chr:
            amnplicons_by_chr[chr] = list()
        amnplicons_by_chr[chr].append({
            "strand": ampli["probe_strand"],
            "up_primer": ampli["ulso_sequence"].upper(),
            "down_primer": ampli["dlso_sequence"].upper(),
            "start": None,
            "end": None
        })

    # Find amplicons coord
    FH_ref = SequenceFileReader.factory( reference_path )
    try:
        for record in FH_ref:
            if record.id in amnplicons_by_chr:
                chr_str = record.string.upper()
                for ampli in amnplicons_by_chr[record.id]:
                    up_primer = ampli["up_primer"]
                    down_primer = ampli["down_primer"]
                    if ampli["strand"] == "-":
                        up_primer = revcom(ampli["down_primer"])
                        down_primer = revcom(ampli["up_primer"])
                    try:
                        ampli["start"] = chr_str.index( up_primer ) + 1
                        ampli["end"] = chr_str.index( down_primer ) + len(down_primer)########## +1                            
                    except:
                        raise Exception("The primers '" + up_primer + "' and '" + down_primer + "' cannot be found in " + record.id )
                    # Check multiple target in chr
                    if chr_str.count(up_primer) > 1:
                        raise Exception("The primer '" + up_primer + "' is found multiple twice in " + record.id )
                    if chr_str.count(down_primer) > 1:
                        raise Exception("The primer '" + down_primer + "' is found multiple twice in " + record.id )
    finally:
        FH_ref.close()
    return( amnplicons_by_chr )

def getVariantsProfile( profile_path ):
    profiles = list()
    with open(profile_path) as FH_profile:
        header = FH_profile.readline()
        for line in FH_profile:
            if line.strip() != "":
                fields = [field.strip() for field in line.split("\t")]
                profiles.append(
                    { "type": fields[0],
                      "occurence": float(fields[1]),
                      "freq": {"min":float(fields[2]), "max":float(fields[3])},
                      "length": {"min":int(fields[4]), "max":int(fields[5])} }
                )
    return profiles

def getCoveredLength( area_by_chr ):
    covered_len = 0
    for chr in area_by_chr:
        for area in area_by_chr[chr]:
            covered_len += area["end"] - area["start"] + 1
    return( covered_len )

def setAmpliconSeq( ref_sequences, area_by_chr ):
    FH_ref = SequenceFileReader.factory( ref_sequences )
    for record in FH_ref:
        if record.id in area_by_chr:
            for area in area_by_chr[record.id]:
                area["ref"] = record.string[(area["start"]-1):area["end"]]
    FH_ref.close()

def getCoveredPos( area_by_chr ):
    covered_pos = list()
    for chr in area_by_chr:
        for area in area_by_chr[chr]:
            for pos in range( area["start"], (area["end"] + 1) ):
                covered_pos.append( chr + ":" + str(pos) )
    return( covered_pos )

def setAlt( variant, region_sequence, start_idx ):
    variant["ref"] = ""
    variant["alt"] = ""
    if variant["type"] == "substitution": # Substitution
        for idx in range(variant["length"]):
            # Ref
            variant["ref"] += region_sequence[start_idx + idx]
            # Alt
            choice_opt = ["A", "T", "G", "C"]
            choice_opt.remove(variant["ref"][-1].upper())
            variant["alt"] += random.choice(choice_opt)
    elif variant["type"] == "deletion": # Deletion
        for idx in range(variant["length"]):
            variant["ref"] += region_sequence[start_idx + idx]
    else: # Insertion
        variant["ref"] = region_sequence[start_idx]
        variant["alt"] = variant["ref"]
        for idx in range(variant["length"]):
            variant["alt"] += random.choice(["A", "T", "G", "C"])

def getFastaPath(out_dir, spl_name, lib_name, read_number):
    return os.path.join( out_dir, spl_name + "-" + lib_name + "_R" + str(read_number) + ".fasta" )

def writeVariants( out_path, variants_by_chr ):
    with open( out_path, "w" ) as FH_out_vcf:
        FH_out_vcf.write( '##fileformat=VCFv4.0\n' )
        FH_out_vcf.write( '##source=simuAmpliDS_IUCT\n' )
        FH_out_vcf.write( '##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">\n' )
        FH_out_vcf.write( '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n' )
        for chr in variants_by_chr:
            for start in sorted(variants_by_chr[chr]):
                variant = variants_by_chr[chr][start]
                FH_out_vcf.write( chr + "\t" + str(start) + "\t.\t" + variant["ref"] + "\t" + variant["alt"] + "\t.\t.\t" + "AF=" + str(variant["freq"]) + "\t.\n" )

def getCoveredFromAmplicons( amplicons_by_chr, remove_primers=True ):
    cover_region_by_chr = dict()
    for chr in amplicons_by_chr:
        cover_region_by_chr[chr] = dict()
        for ampli in amplicons_by_chr[chr]:
            start = ampli["start"]            
            end = ampli["end"]
            if remove_primers:
                up_primer = ampli["up_primer"]
                down_primer = ampli["down_primer"]
                if ampli["strand"] == "-":
                    up_primer = revcom(ampli["down_primer"])
                    down_primer = revcom(ampli["up_primer"])
                start += len(up_primer)
                end -= len(down_primer)
            if start not in cover_region_by_chr[chr]:
                cover_region_by_chr[chr][start] = end
            else:
                if cover_region_by_chr[chr][start] < end:
                    cover_region_by_chr[chr][start] = end
                else:
                    print("inclusion")
    area_by_chr = dict()
    for chr in cover_region_by_chr:
        area_by_chr[chr] = list()
        prev_end = -1
        for start in sorted(cover_region_by_chr[chr]):
            end = cover_region_by_chr[chr][start]
            if start <= prev_end:
                area_by_chr[chr][-1]["end"] = end
            else:
                area_by_chr[chr].append({"start": start, "end": end, "ref":None, "alt":None})
            prev_end = end
    return( area_by_chr )

def writeReadsPair( pair_id, matrix, rvc_matrix, reads_length, FH_R1, FH_R2 ):
    R1 = Sequence( pair_id, matrix[0:reads_length] )
    R2 = Sequence( pair_id, rvc_matrix[0:reads_length] )
    if len(R1.string) < reads_length:
        rvc_p5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTACTTAGCAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        rvc_p7 = "GCGAATTTCGACGATCGTTGCTTAATCGCGACACCACACATCCCGTATGCCGTCTTCTG"
        R1.string += rvc_p7[:(reads_length - len(R1.string))]
        R2.string += rvc_p5[:(reads_length - len(R2.string))]
    FH_R1.write( R1 )
    FH_R2.write( R2 )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='**************************************************************.' )
    parser.add_argument( '-s', '--random-seed', type=int, help="The *****************************************. [Default: auto]" )
    parser.add_argument( '-d', '--min-distance', type=int, default=2, help="The *****************************************. [Default: %(default)s]" )
    parser.add_argument( '-n', '--sample-name', required=True, help="The *****************************************." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--input-manifests', required=True, help='*********************************** (format: *****************).' )
    group_input.add_argument( '-g', '--input-genome', required=True, help='*************************************** (format: fasta).' )
    group_input.add_argument( '-p', '--input-profile', required=True, help='Path to the variants profile file (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--output-variants', required=True, help='*********************************** (format: VCF).' )
    group_output.add_argument( '--output-folder', default=os.getcwd(), help='*************************************** (format: fasta). [Default: %(default)s]' )
    args = parser.parse_args()

    # Get reference sequences
    print( "[GET amplicons]\tStart" )
    amplicons_by_chr = getAmplicons( args.input_genome, args.input_manifests )
    #~ amplicons_by_chr = {'chr15': [{'strand': '+', 'down_primer': 'GCTGGCAATGGCGGTGTGGTGTTCAA', 'up_primer': 'TGATGAGCAGCAGCGAAAGCGCCTTGA', 'start': 66727410, 'end': 66727535}, {'strand': '-', 'down_primer': 'CAAAGTGGGGAGCACAAGTCAATAC', 'up_primer': 'ACCTTCTGCTTCTGGGTAAGAAAGGC', 'start': 66727334, 'end': 66727463}, {'strand': '-', 'down_primer': 'CCCAGCTCACTGATCTTCTCAAAGT', 'up_primer': 'GGCTTCTAAGTACCCTGAGAAATAATCC', 'start': 66727484, 'end': 66727631}], 'chr14': [{'strand': '+', 'down_primer': 'CGCCAGGTCTTGATGTACTCCCCTACA', 'up_primer': 'ACAGAGAAGTTGTTGAGGGGAGCCTCA', 'start': 105246430, 'end': 105246558}, {'strand': '-', 'down_primer': 'TCTTTGCCCACCAGGCCCCACGATG', 'up_primer': 'AAGGAGCGGCCGCAGGATGTGGACCAAC', 'start': 105246346, 'end': 105246485}, {'strand': '-', 'down_primer': 'CCACGCTACTTCCTCCTCAAGAATG', 'up_primer': 'AGTTCCTGCCTGGCTGCCTGGCGAGG', 'start': 105246506, 'end': 105246627}], 'chr7': [{'strand': '+', 'down_primer': 'TCACCCTAAGTACATTACCTTATGCCTG', 'up_primer': 'GTAATCTGCCCATCAGGAATCTCCCAA', 'start': 140481435, 'end': 140481580}, {'strand': '-', 'down_primer': 'ATGTAATGTGGTGACATTGTGACAAGT', 'up_primer': 'ACTTGGTAGACGGGACTCGAGTGATG', 'start': 140481343, 'end': 140481488}, {'strand': '+', 'down_primer': 'GATCTCATTTTCCTATCAGAGCAAGCA', 'up_primer': 'CATCGAGATTTCACTGTAGCTAGACCA', 'start': 140453124, 'end': 140453250}, {'strand': '-', 'down_primer': 'ATTGAGGCTATTTTTCCACTGATTAAATT', 'up_primer': 'GAAGACCTCACAGTAAAAATAGGTGATT', 'start': 140453040, 'end': 140453179}, {'strand': '+', 'down_primer': 'TCTCATGGTCTGGTGGGGAGCCCAGA', 'up_primer': 'AGCTCTCTTGAGGATCTTGAAGGAAAC', 'start': 55241655, 'end': 55241804}, {'strand': '-', 'down_primer': 'GGACAAGAACACAGAGACAAGGGT', 'up_primer': 'GCCCAGCACTTTGATCTTTTTGAATTC', 'start': 55241583, 'end': 55241709}, {'strand': '+', 'down_primer': 'CTGAACCTCAGGCCCACCTTTTCTCAT', 'up_primer': 'GTCGCTATCAAGGAATTAAGAGAAGCA', 'start': 55242454, 'end': 55242576}, {'strand': '-', 'down_primer': 'CAGAGAGAGAAGGAAGACGTTAACTGG', 'up_primer': 'GATTTCCTTGTTGGCTTTCGGAGATG', 'start': 55242382, 'end': 55242507}, {'strand': '+', 'down_primer': 'TGGTGTGTGCAGATCGCAAAGGTAAT', 'up_primer': 'GCTGGGCATCTGCCTCACCTCCA', 'start': 55249033, 'end': 55249176}, {'strand': '-', 'down_primer': 'GAGGGAGAGGCACGTCAGTGTG', 'up_primer': 'ATGAGCTGCGTGATGAGCTGCACG', 'start': 55248957, 'end': 55249080}, {'strand': '-', 'down_primer': 'TTGAGCAGGTACTGGGAGCCAATATTG', 'up_primer': 'CATATCCCCATGGCAAACTCTTGCTAT', 'start': 55249123, 'end': 55249271}, {'strand': '+', 'down_primer': 'AAAGTAAGGAGGTGGCTTTAGGTCAG', 'up_primer': 'TGGCAGCCAGGAACGTACTGGTGAAAA', 'start': 55259455, 'end': 55259590}, {'strand': '-', 'down_primer': 'AGAGAAGACCCTGCTGTGAGGGACAG', 'up_primer': 'AATCTGTGATCTTGACATGCTGCGGT', 'start': 55259379, 'end': 55259508}, {'strand': '-', 'down_primer': 'CCTCCTTCTGCATGGTATTCTTTCTCT', 'up_primer': 'CCCCTGCATGTGTTAAACAATACAGCT', 'start': 55259537, 'end': 55259659}, {'strand': '+', 'down_primer': 'GCTACTTTTCCAGAAGGTATATTTCAGT', 'up_primer': 'GAGTACACACTCCTCATTTGGATAGGC', 'start': 116411936, 'end': 116412055}, {'strand': '-', 'down_primer': 'ACAGAGAGAAAGAAAGAGCTTGTTAAA', 'up_primer': 'TTGGGCTTACACTTCGGGCACTTACA', 'start': 116411870, 'end': 116411989}, {'strand': '-', 'down_primer': 'CGGTAGTCTACAGATTCATTTGAAACCAT', 'up_primer': 'CCAAAAATAAACAACAATGTCACAACC', 'start': 116411998, 'end': 116412120}, {'strand': '+', 'down_primer': 'CTAAATCCAGAGCTGGTCCAGGCAGT', 'up_primer': 'AACGGTTCATGCCGACAAGTGCAGTAT', 'start': 116414958, 'end': 116415103}, {'strand': '+', 'down_primer': 'CAACTGACAGAGCAGTGATAACAAGTG', 'up_primer': 'TTCAATGAAGTCATAGGAAGAGGTAAG', 'start': 116415144, 'end': 116415288}, {'strand': '-', 'down_primer': 'AGAACAACAGTAAAACCTCATTTAATGG', 'up_primer': 'TTAGGATGGGGGACATGTCTGTCAGAG', 'start': 116414882, 'end': 116415012}, {'strand': '-', 'down_primer': 'GCACTGAGGTCAATGTGGACAGTATTT', 'up_primer': 'CGTATTTAACAAAAAGCTGAGTGGAAAT', 'start': 116415050, 'end': 116415199}, {'strand': '+', 'down_primer': 'AACGTTGATTTACACTTTCCCCTTGTG', 'up_primer': 'CAATGATGGCAAGAAAATTCACTGTGC', 'start': 116417480, 'end': 116417620}, {'strand': '-', 'down_primer': 'AAGAAAGAACTTGGTTAGCACTGCG', 'up_primer': 'ATGCCACTTACTGTTCAAGGATTTCAC', 'start': 116417410, 'end': 116417534}, {'strand': '+', 'down_primer': 'TCTTCGAAATTTCATTCGAAATGAGACTC', 'up_primer': 'TGACCGAGGGAATCATCATGAAAGA', 'start': 116418863, 'end': 116419009}, {'strand': '-', 'down_primer': 'TCCAGTGGTGGAGACATTAACTTCATT', 'up_primer': 'AGCGAGAGGACATTGGGATGACTAAA', 'start': 116418795, 'end': 116418914}, {'strand': '-', 'down_primer': 'CTCCATGTTTCATGTATGGTAGGACCA', 'up_primer': 'TTTTGAAGGGATGGCTGGCTTACAGCT', 'start': 116418953, 'end': 116419073}, {'strand': '+', 'down_primer': 'ACAATCCAAATTAAGTGACAAGGAGGA', 'up_primer': 'GTAGCCAAAGGCATGAAATATCTTGCA', 'start': 116422081, 'end': 116422201}, {'strand': '-', 'down_primer': 'CAAAGAGAGTTAGAAAAGCATAAACTAAGC', 'up_primer': 'CCAAGTCTCTGTGGACAAACTTTTTGC', 'start': 116422005, 'end': 116422135}, {'strand': '+', 'down_primer': 'GTGGTAATGTATTGGTTATCTCTGAGT', 'up_primer': 'TTGCCAGAGACATGTATGATAAAGAATAC', 'start': 116423399, 'end': 116423547}, {'strand': '-', 'down_primer': 'GGACAAAAATTATTACCCGTGGCTGA', 'up_primer': 'TTGCACCTGTTTTGTTGTGTACACTAT', 'start': 116423323, 'end': 116423455}, {'strand': '-', 'down_primer': 'TCTGACTTGGTGGTAAACTTTTGAGTT', 'up_primer': 'GAAAACTGGAATTGGTGGTGTTGAATT', 'start': 116423493, 'end': 116423639}, {'strand': '+', 'down_primer': 'TTTACAGAAATGCCTGCCTTCAAAGGG', 'up_primer': 'CGTAAACACCTTTGATATAACTGTTTACT', 'start': 116435768, 'end': 116435908}, {'strand': '-', 'down_primer': 'ACAAGAAACAGGACAGATGAGGTGA', 'up_primer': 'TTGTAGGAGTCTTCTCCCTTGCAAC', 'start': 116435676, 'end': 116435822}, {'strand': '+', 'down_primer': 'CCCGCTGTGCTTGCACCTGGCAT', 'up_primer': 'AAATGGTATAGGTCTTTCAGTTTTCTCTTC', 'start': 116339050, 'end': 116339170}, {'strand': '+', 'down_primer': 'TTCACCGCGGAAACACCCATCCAGAAT', 'up_primer': 'CCTGTTTACCTTGGTGCAGAGGAGCAAT', 'start': 116339180, 'end': 116339300}, {'strand': '+', 'down_primer': 'CCTGTGCTGGAACACCCAGATTGTTT', 'up_primer': 'TCTACATGAGCATCACATTTTCCTTGG', 'start': 116339306, 'end': 116339425}, {'strand': '+', 'down_primer': 'TTAGCTGTGGCAGCGTCAACAGA', 'up_primer': 'TCAGGACTGCAGCAGCAAAGCCAATTTA', 'start': 116339432, 'end': 116339552}, {'strand': '+', 'down_primer': 'GACTGTGTGGTGAGCGCCCTGGGA', 'up_primer': 'CAGCGACATGTCTTTCCCCACAATCATA', 'start': 116339562, 'end': 116339681}, {'strand': '+', 'down_primer': 'GTGAGAAGGCTAAAGGAAACGAAAGATG', 'up_primer': 'TTTCATCTGTAAAGGACCGGTTCATCA', 'start': 116339692, 'end': 116339811}, {'strand': '+', 'down_primer': 'TTTATTTACTTCTTGACGGTCCAAAGGG', 'up_primer': 'TGTTTTTGACGGACCAGTCCTACATTG', 'start': 116339818, 'end': 116339937}, {'strand': '+', 'down_primer': 'AAAAGAGATCCACAAAGAAGGAAGTGT', 'up_primer': 'CAGACTTTTCACACAAGAATAATCAGG', 'start': 116339952, 'end': 116340078}, {'strand': '+', 'down_primer': 'CCTATCAAATATGTCAACGACTTCTTC', 'up_primer': 'CCCAGCTTGCTAGACAAATAGGAGC', 'start': 116340118, 'end': 116340260}, {'strand': '+', 'down_primer': 'GGCTTCTTTTGTGCTTTGTAAATGGTG', 'up_primer': 'CATTTTTACGGACCCAATCATGAGCAC', 'start': 116340300, 'end': 116340444}, {'strand': '-', 'down_primer': 'GGACTTTTACATAGAAGGAAAACAAAAACT', 'up_primer': 'TTCAAGGCGAGAGCAGTTCAGTTGTCAG', 'start': 116338992, 'end': 116339111}, {'strand': '-', 'down_primer': 'CCTTCATTATGAGAGGTTTATCTGCCAAAA', 'up_primer': 'ACTTTGCTAGTGCCTCTTTACACTCC', 'start': 116339116, 'end': 116339235}, {'strand': '-', 'down_primer': 'GGGAAGCTGATACTTCATATTCACATT', 'up_primer': 'CTCATTTAAAACATAAATGTAGTTAGTGGC', 'start': 116339244, 'end': 116339363}, {'strand': '-', 'down_primer': 'CAGTCTTGTACTCAGCAACCTTCTGAAG', 'up_primer': 'GTTGATGTTATCTTTCCAAACACCTCCT', 'start': 116339370, 'end': 116339489}, {'strand': '-', 'down_primer': 'GTTGATCATCATAGTAGGTGTCGACAA', 'up_primer': 'TATGCAGTGAACCTCCGACTGTATGTCA', 'start': 116339500, 'end': 116339621}, {'strand': '-', 'down_primer': 'ACACTGGCTGGGCTCTTCTATCTGTGG', 'up_primer': 'GAAGAATTTATGGTATTGCCTACAAAGAAG', 'start': 116339628, 'end': 116339749}, {'strand': '-', 'down_primer': 'GATATCGAATGCAATGGATGATCTGGG', 'up_primer': 'GGGTAAGAATCTCTGAACTCAGGTAAAACA', 'start': 116339756, 'end': 116339875}, {'strand': '-', 'down_primer': 'TTGTTGCTTTCAAAGGCATGGACAT', 'up_primer': 'GGAATGCAATCCAGAGTTTATGGAACAGA', 'start': 116339884, 'end': 116340008}, {'strand': '-', 'down_primer': 'TCTCTTTTCTGTGAGAATACACTCCAG', 'up_primer': 'CCCGAAAAGAATGTCATCATTCAGGCT', 'start': 116340024, 'end': 116340170}, {'strand': '-', 'down_primer': 'AATGCACACATGGCAGATCGATC', 'up_primer': 'GGAACTGATGTGACTTACCCTATTAAAGC', 'start': 116340210, 'end': 116340356}], 'chr3': [{'strand': '+', 'down_primer': 'GTCACAGGTAAGTGCTAAAATGGAGAT', 'up_primer': 'AGTAACAGACTAGCTAGAGACAATGAA', 'start': 178935998, 'end': 178936142}, {'strand': '-', 'down_primer': 'TCACAGATGATTTACAGAAAAAGCAAAT', 'up_primer': 'TGAGCTGTTCTTTGTCATTTTCCCTTA', 'start': 178935912, 'end': 178936052}, {'strand': '-', 'down_primer': 'CCATAGAAAATCTTTCTCCTGCTCAGT', 'up_primer': 'GGTATGGTAAAAACATGCTGAGATCAG', 'start': 178936088, 'end': 178936214}, {'strand': '+', 'down_primer': 'GCATGCCAATCTCTTCATAAATCTTTTCTC', 'up_primer': 'TTAACATCATTTGCTCCAAACTGACCA', 'start': 178951834, 'end': 178951953}, {'strand': '+', 'down_primer': 'ATTTCATGAAACAAATGAATGATGCAC', 'up_primer': 'ATGCTTGGCTCTGGAATGCCAGAACTA', 'start': 178951958, 'end': 178952084}, {'strand': '+', 'down_primer': 'CACTGCACTGTTAATAACTCTCAGCAG', 'up_primer': 'GCTGGACAACAAAAATGGATTGGATCT', 'start': 178952094, 'end': 178952218}, {'strand': '+', 'down_primer': 'TGATAGCACTTAAACTAGTTCATTTCAAA', 'up_primer': 'GATTGCATAGGAATTGCACAATCCATG', 'start': 178952228, 'end': 178952358}, {'strand': '+', 'down_primer': 'TATAGAAATGATGGAGAAGGAAAAAGTGA', 'up_primer': 'ATAATGCGCAATTTCATGTTATGCCTT', 'start': 178952372, 'end': 178952504}, {'strand': '-', 'down_primer': 'GAGAGGCTTTCGTAGACAAAGCTATAA', 'up_primer': 'CTCCTGAAACCTATAAGTAATAAGAACAGT', 'start': 178951772, 'end': 178951891}, {'strand': '-', 'down_primer': 'GTCGAATAGCTAGATAAGCCTTGTAACAC', 'up_primer': 'CGAATGTATGCAATGTCATCAAAAGATT', 'start': 178951894, 'end': 178952013}, {'strand': '-', 'down_primer': 'CTCCAAAGCCTCTTGCTCAGTTTTATC', 'up_primer': 'GTTCAATGCATGCTGTTTAATTGTGTGG', 'start': 178952030, 'end': 178952149}, {'strand': '-', 'down_primer': 'GGAATCCAGAGTGAGCTTTCATTTTCT', 'up_primer': 'CTGTTCTTGCTGTAAATTCTAATGCTGT', 'start': 178952164, 'end': 178952283}, {'strand': '-', 'down_primer': 'ACCCTGTTTGCGTTTACATTATTTAAA', 'up_primer': 'CAATCTTCAAAGTTTACCTTTTTGGACT', 'start': 178952302, 'end': 178952427}, {'strand': '-', 'down_primer': 'ATTTTGGGGATTTTTGTTTTGTTTTGT', 'up_primer': 'CAATAACAGCCTTTGTTGTGTCCACAT', 'start': 178952448, 'end': 178952568}], 'chr11': [{'strand': '+', 'down_primer': 'AGGGTCTCCTGCCCCACCTGCCAA', 'up_primer': 'GATCAGCTGGATGGTCAGCGCACTCTT', 'start': 534251, 'end': 534382}, {'strand': '-', 'down_primer': 'TGGCGCCGCCGTCCAGGTGCCAGC', 'up_primer': 'TGGTGGTGGTGGGCGCCGGCGGTGTGGG', 'start': 534181, 'end': 534306}, {'strand': '-', 'down_primer': 'TAGGAGGACCCCGGGCCGCAGGCCCCT', 'up_primer': 'TTTGCCCTTCAGATGGCCCTGCCAGC', 'start': 534331, 'end': 534462}, {'strand': '+', 'down_primer': 'GGCGAGGTCCTGAGCCTGCCGAGATT', 'up_primer': 'GGCGGGGCGGGTCCCTGGCTAGCTGT', 'start': 533407, 'end': 533526}, {'strand': '+', 'down_primer': 'GGGGCTGCAGGCGCAGCGGCATCCAG', 'up_primer': 'CAGGTCACACTTGTTCCCCACCA', 'start': 533543, 'end': 533668}, {'strand': '+', 'down_primer': 'TCCTCAAAAGACTTGGTGTTGTTGATG', 'up_primer': 'ACAGGAGGCCCCTGCCTGGACGCA', 'start': 533685, 'end': 533807}, {'strand': '+', 'down_primer': 'CTGCAGGAGGACAGGGCTCAGGGA', 'up_primer': 'GTCCCGCATGGCGCTGTACTCCT', 'start': 533849, 'end': 533968}, {'strand': '-', 'down_primer': 'ACAGGGCAGCCGCTCTGGCTCTAGCT', 'up_primer': 'GCAGGTGAGGCAGCTCTCCACC', 'start': 533337, 'end': 533456}, {'strand': '-', 'down_primer': 'AAGCTACGGCATCCCCTACATCGAGA', 'up_primer': 'AAGGACTCGGATGACGTGCCCATGGTG', 'start': 533473, 'end': 533593}, {'strand': '-', 'down_primer': 'TAGCCAGCTCTCGCTTTCCACCTCTC', 'up_primer': 'ACGCCGCACAGGTGGGGCCAGGCCG', 'start': 533615, 'end': 533735}, {'strand': '-', 'down_primer': 'ATCCACCAGTACAGGTGAACCCCGTGA', 'up_primer': 'TTGGACATCCTGGATACCGCCGGCCAG', 'start': 533753, 'end': 533899}, {'strand': '-', 'down_primer': 'ATTCCTACCGGAAGCAGGTGGTCATTG', 'up_primer': 'ATGAGGGGCATGAGAGGTACCAGGGA', 'start': 533917, 'end': 534038}], 'chr1': [{'strand': '+', 'down_primer': 'CGAGCCACATCTACAGTACTTTAAAGC', 'up_primer': 'CGCTTTTCCCAACACCACCTGCTCCAA', 'start': 115258730, 'end': 115258856}, {'strand': '-', 'down_primer': 'CCAGTGGTAGCCCGCTGACCTG', 'up_primer': 'GAAATGACTGAGTACAAACTGGTGGTG', 'start': 115258642, 'end': 115258784}, {'strand': '+', 'down_primer': 'TTTCTGTAAGAATCCTGGGGGTGTG', 'up_primer': 'GGAAGCCTTCGCCTGTCCTCATGTATT', 'start': 115256476, 'end': 115256610}, {'strand': '-', 'down_primer': 'GGAGCATTATTTTCTCTGAAAGGATGA', 'up_primer': 'CAAGAAGAGTACAGTGCCATGAGAGAC', 'start': 115256388, 'end': 115256530}, {'strand': '-', 'down_primer': 'CAAGTGGTTATAGATGGTGAAACCTGT', 'up_primer': 'AGATGCTTATTTAACCTTGGCAATAGC', 'start': 115256558, 'end': 115256700}, {'strand': '+', 'down_primer': 'GAGTCTTTTACTCGCTTAATCTGCTCCC', 'up_primer': 'AACTCTTGGCCAGTTCGTGGGCTTGTTTT', 'start': 115252231, 'end': 115252350}, {'strand': '-', 'down_primer': 'CAGCTTTCAGCATTTGTGCAAGAGTTT', 'up_primer': 'GTGATTTGCCAACAAGGACAGTTGATA', 'start': 115252155, 'end': 115252287}, {'strand': '-', 'down_primer': 'ATGATGTACCTATGGTGCTAGTGGGA', 'up_primer': 'ATGAGCCACTGTACCCAGCCTAATCT', 'start': 115252295, 'end': 115252414}], 'chr4': [{'strand': '+', 'down_primer': 'ATGCCAGTAGGACGCCTGGCGCCAACA', 'up_primer': 'CGAGGACAACGTGATGAAGATCGCAGA', 'start': 1807819, 'end': 1807961}, {'strand': '-', 'down_primer': 'AAGGGGTGGGAGGCAGGGCTGAAGC', 'up_primer': 'TTGTGCACGTCCCGGGCCAGCCCGAA', 'start': 1807747, 'end': 1807872}, {'strand': '+', 'down_primer': 'CGTGCTCAAGGTGGGCCACCGTGT', 'up_primer': 'CTGGGCAGCGACGTGGAGTTCCACT', 'start': 1803621, 'end': 1803766}, {'strand': '-', 'down_primer': 'GGCAGATGACGCTCAGGGGCCA', 'up_primer': 'CTGTGCGTCACTGTACACCTTG', 'start': 1803533, 'end': 1803668}, {'strand': '-', 'down_primer': 'TAACGTAGGGTGTGCCGTCCGGGCCCA', 'up_primer': 'AAGCTCCAACCCCTAGACCCAAATCCT', 'start': 1803715, 'end': 1803861}, {'strand': '+', 'down_primer': 'GCGGCTGTGACGCTCTGCCGCCT', 'up_primer': 'TTTGCAGCCGAGGAGGAGCTGGT', 'start': 1806050, 'end': 1806174}, {'strand': '+', 'down_primer': 'GGTCCTGGGCTGTGTGAGCCCTCTCT', 'up_primer': 'GCACAAGATCTCCCGCTTCCCGCTCAA', 'start': 1806214, 'end': 1806347}, {'strand': '-', 'down_primer': 'ATGGATGCCCCTGGCCCAGAGCCGCA', 'up_primer': 'ATACACACTGCCCGCCTCGTCAGCCTC', 'start': 1805976, 'end': 1806100}, {'strand': '-', 'down_primer': 'ACCACCAGGATGAACAGGAAGAAGC', 'up_primer': 'CCTGGTATCTACTTTCTGTTACCTGTCG', 'start': 1806126, 'end': 1806269}, {'strand': '+', 'down_primer': 'CCTTTTTGGGTACACATAACAGTGACT', 'up_primer': 'GTTTACATAGACCCAACACAACTTCCT', 'start': 55593639, 'end': 55593765}, {'strand': '-', 'down_primer': 'GAGAAAGGGAAAAATAGATCACCTTTTAAT', 'up_primer': 'TCTGGGAAACTCCCATTTGTGATCAT', 'start': 55593545, 'end': 55593692}, {'strand': '+', 'down_primer': 'TGTCTGTCAGGTTATCAAAACATGACA', 'up_primer': 'AAGTCCTGAGTTACCTTGGTAATCACA', 'start': 55594222, 'end': 55594342}, {'strand': '-', 'down_primer': 'TTGGAAACATGCATTTTAGCAAAAAGC', 'up_primer': 'GCAGGCTCCAAGTAGATTCACAATATTC', 'start': 55594144, 'end': 55594277}, {'strand': '+', 'down_primer': 'CTTGACAGTCCTGCAAAGGATTTTTAG', 'up_primer': 'CTCCTTACTCATGGTCGGATCACAAAG', 'start': 55599269, 'end': 55599403}, {'strand': '-', 'down_primer': 'GGTTGGAGGAGAAAAGAAAACCATTTA', 'up_primer': 'TGTCTCTGGCTAGACCAAAATCACAAA', 'start': 55599203, 'end': 55599323}, {'strand': '+', 'down_primer': 'CTTCTTGAAGTTTCATTGGTGTCCTGC', 'up_primer': 'TGTGTATACACGTTTGAAAGTGACGTC', 'start': 55602709, 'end': 55602847}, {'strand': '-', 'down_primer': 'TAGAAGCACAACAGAGTCAATAATGTT', 'up_primer': 'GCTCCCAAAGAAAAATCCCATAGGACC', 'start': 55602631, 'end': 55602763}, {'strand': '+', 'down_primer': 'GAAGGACTGCAATTCACTTGAATTTCA', 'up_primer': 'ATGCTCCAATGTGTGGCAGCAGGATT', 'start': 55589791, 'end': 55589931}, {'strand': '-', 'down_primer': 'ACAGAAATGGCCATATGTCAGAGTG', 'up_primer': 'CAAAAATACCAATCTATTGTGGGCTCTG', 'start': 55589715, 'end': 55589846}, {'strand': '+', 'down_primer': 'CAAGACTTCTGCCTATTTTAACTTTGC', 'up_primer': 'ATCTGGGCCACCGTTTGGAAAGCTAGT', 'start': 55592068, 'end': 55592196}, {'strand': '-', 'down_primer': 'GGAAGAAAACAAAAGCCCTGGCTTACT', 'up_primer': 'TGCACTAGAATCTATAGAACTCTGAAC', 'start': 55591990, 'end': 55592122}, {'strand': '-', 'down_primer': 'CCACATCGTTGTAAGCCTTACATTCAA', 'up_primer': 'ATCATGACTGATATGGTAGACAGAGCC', 'start': 55592142, 'end': 55592280}, {'strand': '+', 'down_primer': 'CTTTTCCCTTGCACACAACTTTACAAT', 'up_primer': 'GGACATGAATATATTTATGTGGACCCG', 'start': 55141059, 'end': 55141199}, {'strand': '-', 'down_primer': 'TGACCAGGACAGGTAACTGGTGAATT', 'up_primer': 'CATCTTGAGTCATAAGGCAGCTGCA', 'start': 55140975, 'end': 55141111}, {'strand': '+', 'down_primer': 'GATTTTCACTGGACACATGTGGTTGTG', 'up_primer': 'TGAACTGAAGATAATGACTCACCTGGG', 'start': 55144100, 'end': 55144224}, {'strand': '-', 'down_primer': 'AGAATAAATCACATATCAGTCCAGCTG', 'up_primer': 'CCAGCAAGTTTACAATGTTCAAATGTGG', 'start': 55144028, 'end': 55144155}, {'strand': '+', 'down_primer': 'TCAGGCTCATCCTCCTTCACTTTAATC', 'up_primer': 'AAATTGTGAAGATCTGTGACTTTGGCC', 'start': 55152057, 'end': 55152181}, {'strand': '-', 'down_primer': 'GGAAAAGGAAGAAATGACTCAGGATCA', 'up_primer': 'TTCGAATCATGCATGATGTCTCTGGCC', 'start': 55151975, 'end': 55152111}], 'chr10': [{'strand': '+', 'down_primer': 'TCAGTTCATTTCCTCTAACTCATGGGC', 'up_primer': 'TGTTTCCCAATCATCTTCATCATCTCC', 'start': 123258050, 'end': 123258172}, {'strand': '-', 'down_primer': 'GGAGGAAAAACTGCATTCGCCCAAAT', 'up_primer': 'AGACCTTTCTGATCTGGTGTCAGAGA', 'start': 123257976, 'end': 123258103}, {'strand': '+', 'down_primer': 'CAGAACAAAAAGGAAATATGTTCATTTCT', 'up_primer': 'AAAGTCTGCTATTTTCATCACATTGTT', 'start': 123247556, 'end': 123247704}, {'strand': '-', 'down_primer': 'ATGGCAGTAACACAGTGGGCAGGGG', 'up_primer': 'GCAGCCAGAAATGTTTTGGTAACAGA', 'start': 123247472, 'end': 123247609}, {'strand': '+', 'down_primer': 'CCGGCAGTCCGGCTTGGAGGAT', 'up_primer': 'CGTTCTTTTCCACGTGCTTGATCCA', 'start': 123279540, 'end': 123279663}, {'strand': '-', 'down_primer': 'TGAATCTAAAGGTACCCACAACTGGGG', 'up_primer': 'TTTACAGTGATGCCCAGCCCCACATCCA', 'start': 123279454, 'end': 123279593}, {'strand': '-', 'down_primer': 'AAATGCCTCCACAGTGGTCGGAGGAGA', 'up_primer': 'TTTATTGGTCTCTCATTCTCCCATCCC', 'start': 123279614, 'end': 123279740}], 'chr2': [{'strand': '+', 'down_primer': 'ACTGGGAGGAACAGGATACAAAGTTAC', 'up_primer': 'AACCGGGGCAGGGATTGCAGGCTCA', 'start': 29443639, 'end': 29443787}, {'strand': '-', 'down_primer': 'AGAACCAGTCTTTGCTGCAGTTGTTG', 'up_primer': 'AGAACATTGTTCGCTGCATTGGG', 'start': 29443539, 'end': 29443687}, {'strand': '+', 'down_primer': 'TGAAATGCTGGCTTCCAGTGCTCACAA', 'up_primer': 'AGAAGGTCCAGCATGGCCAGGGA', 'start': 29436913, 'end': 29437033}, {'strand': '-', 'down_primer': 'AAAGTGACTGTGCTCTTCCTGTCATC', 'up_primer': 'TTCTGTCTCCCCACAGAGCCAGCCCTC', 'start': 29436817, 'end': 29436963}, {'strand': '+', 'down_primer': 'AACGGCCATCACTAGGATTTTATCTCC', 'up_primer': 'GCCTGGACAGGTCAAGAGGCAGTTTCT', 'start': 29432705, 'end': 29432825}, {'strand': '-', 'down_primer': 'CTCACCCCTCCGGGCCTGTCTC', 'up_primer': 'TTCCTTTCTTCCCAGAGACATTGCTGC', 'start': 29432615, 'end': 29432759}, {'strand': '+', 'down_primer': 'ATTAGCAGTAGGGGTGTTATCAATGCA', 'up_primer': 'ATGGCTTCAATTGCATTGTAAGGGTCC', 'start': 212570090, 'end': 212570238}, {'strand': '-', 'down_primer': 'TTTTAAACTGCAACTGTTGCGCTAT', 'up_primer': 'TTGGCTACATCTCTTCTTGATTTTCTA', 'start': 212570012, 'end': 212570144}, {'strand': '+', 'down_primer': 'CAGGGACTGGAACTGTAGAGAGGTGAT', 'up_primer': 'CTGTTGTCCCGGATTACTATTCTCTGG', 'start': 212566709, 'end': 212566857}, {'strand': '-', 'down_primer': 'TTTCCCATTCTTCTTTGGACCAAACAA', 'up_primer': 'CTGGACAACACTCTTCAGCACAATCA', 'start': 212566621, 'end': 212566762}, {'strand': '-', 'down_primer': 'AGGAAATCAGCGCAGGAAACATCTATA', 'up_primer': 'CCTTATTTTACTCCGGAATGCGTTTCA', 'start': 212566803, 'end': 212566949}], 'chr17': [{'strand': '+', 'down_primer': 'CCAAGGTATGCACCTGGGCTCTTT', 'up_primer': 'TGCCTGACATCCACGGTGCAGCT', 'start': 37881036, 'end': 37881183}, {'strand': '-', 'down_primer': 'ACAAGGGTACGCTGAGAGGGTATGGGA', 'up_primer': 'AGCCATAGGGCATAAGCTGTGTCAC', 'start': 37880946, 'end': 37881084}, {'strand': '-', 'down_primer': 'AATCTGCATACACCAGTTCAGCAGGTCCT', 'up_primer': 'ACATGCCCAGCAAGAGTCCCCATCCTA', 'start': 37881130, 'end': 37881250}], 'chr12': [{'strand': '+', 'down_primer': 'GAACATGTCACACATAAGGTTAATACAC', 'up_primer': 'TCTGAATTAGCTGTATCGTCAAGGCAC', 'start': 25398243, 'end': 25398392}, {'strand': '-', 'down_primer': 'CTTGTTTTAATATGCATATTACTGGTGCA', 'up_primer': 'TAGTTGGAGCTGGTGGCGTAGGCAAG', 'start': 25398173, 'end': 25398296}, {'strand': '+', 'down_primer': 'CTGTGTCGAGAATATCCAAGAGACAGG', 'up_primer': 'CACCTATAATGGTGAATATCTTCAAATGA', 'start': 25380165, 'end': 25380309}, {'strand': '-', 'down_primer': 'TTGAGTCTTTGCTAATGCCATGC', 'up_primer': 'GCTTTCTTTGTGTATTTGCCATAAATAATA', 'start': 25380087, 'end': 25380228}, {'strand': '-', 'down_primer': 'AGGTCAAGAGGAGTACAGTGCAATGAG', 'up_primer': 'GGTGCACTGTAATAATCCAGACTGTGT', 'start': 25380255, 'end': 25380387}, {'strand': '+', 'down_primer': 'CCTTAACTCTTTTAATTTGTTCTCTGGG', 'up_primer': 'CCATAACTTCTTGCTAAGTCCTGAGCC', 'start': 25378585, 'end': 25378712}, {'strand': '-', 'down_primer': 'ACACTGAAATAAATACAGATCTGTTTTCTG', 'up_primer': 'TGCCTTCTAGAACAGTAGACACAAAAC', 'start': 25378511, 'end': 25378639}, {'strand': '-', 'down_primer': 'CTCTGAAGATGTACCTATGGTCCTAGT', 'up_primer': 'ATGACAAAAGTTGTGGACAGGTTTTGA', 'start': 25378657, 'end': 25378785}]}
    
    print( "[GET amplicons]\tEnd" ) ############################### Sortir ce code avec la recup des sequences dans un autre script

    # Get covered positions
    print( "[GET coverage]\tStart" )
    area_by_chr = getCoveredFromAmplicons( amplicons_by_chr, False )
    # Get covered length
    #~ covered_len = getCoveredLength(area_by_chr)
    #~ print( "[GET coverage]\t\t|- Theoriticall covered length:\t" + str(covered_len) )
    print( "[GET coverage]\tEnd" )
    #################################

    # Get reference sequences
    print( "[GET sequences]\tStart" ) ############################################### Ne faire la recherche qu une fois pour les amplicon puis fusionner les region a ce moment du code
    setAmpliconSeq( args.input_genome, area_by_chr )
    #~ area_by_chr = {'chr14': [{'alt': None, 'ref': 'CATCGTGGGGCCTGGTGGGCAAAGAGGGCTCCAGCCAACCCCCCAAATCTGAATCCCGAGAGGCCAAGGGGATACTTACGCGCCACAGAGAAGTTGTTGAGGGGAGCCTCACGTTGGTCCACATCCTGCGGCCGCTCCTTGTAGCCAATGAAGGTGCCATCATTCTTGAGGAGGAAGTAGCGTGGCCGCCAGGTCTTGATGTACTCCCCTACAGACGTGCGGGTGGTGAGAGCCACGCACACTCTACCCGTCAGACCCTCGCCAGGCAGCCAGGCAGGAACT', 'end': 105246627, 'start': 105246346}], 'chr10': [{'alt': None, 'ref': 'CCCCTGCCCACTGTGTTACTGCCATCGACTTACATTGGTGGTCTTTTTGTAATAGTCTATATTGTTGATATCTCTGGCGAGTCCAAAGTCTGCTATTTTCATCACATTGTTTTCTGTTACCAAAACATTTCTGGCTGCTAAATCTCGATGAATACACTGAAATCAAGAAAGAAGCAAGAGAAATAACTAATTTCAAAACACCGCCAGAACAAAAAGGAAATATGTTCATTTCT', 'end': 123247704, 'start': 123247472}, {'alt': None, 'ref': 'ATTTGGGCGAATGCAGTTTTTCCTCCTACTCACCATCCTGTGTGCAGGCTCCAAGAAGATTTATGATATTCTTGTGTTTCCCAATCATCTTCATCATCTCCATCTCTGACACCAGATCAGAAAGGTCTTTCTCTGTGGCATCATCTATGAACAGTAGGCATATTCACAAATCAGTTCATTTCCTCTAACTCATGGGC', 'end': 123258172, 'start': 123257976}, {'alt': None, 'ref': 'CCCCAGTTGTGGGTACCTTTAGATTCAGAAAGTCCTCACCTTGAGAACCTTGAGGTAGGGCAGCCCGTCGGGCCCGTATTTACTGCCGTTCTTTTCCACGTGCTTGATCCACTGGATGTGGGGCTGGGCATCACTGTAAACCTTGCAGACAAACTCTACGTCTCCTCCGACCACTGTGGAGGCATTTGCCGGCAGTCCGGCTTGGAGGATGGGCCGGTGAGGCGATCGCTCTGGTGGAGAGAGGGAAGAAAGGAGGAGTGGGGATGGGAGAATGAGAGACCAATAAA', 'end': 123279740, 'start': 123279454}], 'chr12': [{'alt': None, 'ref': 'CAGAAAACAGATCTGTATTTATTTCAGTGTTACTTACCTGTCTTGTCTTTGCTGATGTTTCAATAAAAGGAATTCCATAACTTCTTGCTAAGTCCTGAGCCTGTTTTGTGTCTACTGTTCTAGAAGGCAAATCACATTTATTTCCTACTAGGACCATAGGTACATCTTCAGAGTCCTTAACTCTTTTAATTTGTTCTCTGGGAAAGAAAAAAAAGTTATAGCACAGTCATTAGTAACACAAATATCTTTCAAAACCTGTCCACAACTTTTGTCAT', 'end': 25378785, 'start': 25378511}, {'alt': None, 'ref': 'gcatggcattagcaaagACTCAAAAAATAAAAACTATAATTACTCCTTAATGTCAGCTTATTATATTCAATTTAAACCCACCTATAATGGTGAATATCTTCAAATGATTTAGTATTATTTATGGCAAATACACAAAGAAAGCCCTCCCCAGTCCTCATGTACTGGTCCCTCATTGCACTGTACTCCTCTTGACCTGCTGTGTCGAGAATATCCAAGAGACAGGTTTCTCCATCAATTACTACTTGCTTCCTGTAGGAATCCTGAGAAGGGAGAAACACAGTCTGGATTATTACAGTGCACC', 'end': 25380387, 'start': 25380087}, {'alt': None, 'ref': 'TGCACCAGTAATATGCATATTAAAACAAGATTTACCTCTATTGTTGGATCATATTCGTCCACAAAATGATTCTGAATTAGCTGTATCGTCAAGGCACTCTTGCCTACGCCACCAGCTCCAACTACCACAAGTTTATATTCAGTCATTTTCAGCAGGCCTTATAATAAAAATAATGAAAATGTGACTATATTAGAACATGTCACACATAAGGTTAATACAC', 'end': 25398392, 'start': 25398173}], 'chr1': [{'alt': None, 'ref': 'AAACTCTTGCACAAATGCTGAAAGCTGTACCATACCTGTCTGGTCTTGGCTGAGGTTTCAATGAATGGAATCCCGTAACTCTTGGCCAGTTCGTGGGCTTGTTTTGTATCAACTGTCCTTGTTGGCAAATCACACTTGTTTCCCACTAGCACCATAGGTACATCATCCGAGTCTTTTACTCGCTTAATCTGCTCCCTAAAAACGGGAATATATTATCAGAACATAAGAAAAACAAGATTAggctgggtacagtggctcat', 'end': 115252414, 'start': 115252155}, {'alt': None, 'ref': 'TCATCCTTTCAGAGAAAATAATGCTCCTAGTACCTGTAGAGGTTAATATCCGCAAATGACTTGCTATTATTGATGGCAAATACACAGAGGAAGCCTTCGCCTGTCCTCATGTATTGGTCTCTCATGGCACTGTACTCTTCTTGTCCAGCTGTATCCAGTATGTCCAACAAACAGGTTTCACCATCTATAACCACTTGTTTTCTGTAAGAATCCTGGGGGTGTGGAGGGTAAGGGGGCAGGGAGGGAGGGAAGTTCAATTTTTATTAAAAACCACAGGGAATGCAATGCTATTGCCAAGGTTAAATAAGCATCT', 'end': 115256700, 'start': 115256388}, {'alt': None, 'ref': 'CAGGTCAGCGGGCTACCACTGGGCCTCACCTCTATGGTGGGATCATATTCATCTACAAAGTGGTTCTGGATTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTTCACACCAGCAAGAACCTGTTGGAAACCAGTAATCAGGGTTAATTGGCGAGCCACATCTACAGTACTTTAAAGC', 'end': 115258856, 'start': 115258642}], 'chr15': [{'alt': None, 'ref': 'GTATTGACTTGTGCTCCCCACTTTGGAACAGGACCAACTTGGAGGCCTTGCAGAAGAAGCTGGAGGAGCTAGAGCTTGATGAGCAGCAGCGAAAGCGCCTTGAGGCCTTTCTTACCCAGAAGCAGAAGGTGGGAGAACTGAAGGATGACGACTTTGAGAAGATCAGTGAGCTGGGGGCTGGCAATGGCGGTGTGGTGTTCAAGGTCTCCCACAAGCCTTCTGGCCTGGTCATGGCCAGAAAGGTGAGTTTGCCTTGATTAACAGGTAATTGGATTATTTCTCAGGGTACTTAGAAGCC', 'end': 66727631, 'start': 66727334}], 'chr2': [{'alt': None, 'ref': 'GAGACAGGCCCGGAGGGGTGAGGCAGTCTTTACTCACCTGTAGATGTCTCGGGCCATCCCGAAGTCTCCAATCTTGGCCACTCTTCCAGGGCCTGGACAGGTCAAGAGGCAGTTTCTGGCAGCAATGTCTCTGGGAAGAAAGGAAATGCATTTCCTAATTTTATCCCTAGGAAGATGAGTGTACAACGGCCATCACTAGGATTTTATCTCC', 'end': 29432825, 'start': 29432615}, {'alt': None, 'ref': 'GATGACAGGAAGAGCACAGTCACTTTGACTCACCGGTGGATGAAGTGGTTTTCCTCCAAATACTGACAGCCACAGGCAATGTCCCGAGCCACGTGCAGAAGGTCCAGCATGGCCAGGGAGGAGGGCTGGCTCTGTGGGGAGACAGAAGCGGGCCACTGACGAGGAGCTTGTCAGTGAGAGGAGGGAAATCTGAAATGCTGGCTTCCAGTGCTCACAA', 'end': 29437033, 'start': 29436817}, {'alt': None, 'ref': 'CAACAACTGCAGCAAAGACTGGTTCTCACTCACCGGGCGAGGGCGGGTCTCTCGGAGGAAGGACTTGAGGTCTCCCCCCGCCATGAGCTCCAGCAGGATGAACCGGGGCAGGGATTGCAGGCTCACCCCAATGCAGCGAACAATGTTCTGGTGGTTGAATTTGCTGCAGAGCAGAGAGGGATGTAACCAAAATTAACTGAGCTGAGTCTGGGCAAATCTTAAACTGGGAGGAACAGGATACAAAGTTAC', 'end': 29443787, 'start': 29443539}, {'alt': None, 'ref': 'TTGTTTGGTCCAAAGAAGAATGGGAAAAAATTTAAGTTTCTATGTTTTAAATGTCTGAGTAATGTACTTACTACAATTTTCAGCTTTTCTGTTGTCCCGGATTACTATTCTCTGGTTGATTGTGCTGAAGAGTGTTGTCCAGTTAATGGTATGATAATAACACAGGTTGCTGTTGTCAGTAATATAGATGTTTCCTGCGCTGATTTCCTTCAGGGACTGGAACTGTAGAGAGGTGATGCCCTGTTGCTTGAGGATAAGCAAGGACAGGCCACTAAGGAGGGGGAAGTGAGAAAACGGAACCATGAAACGCATTCCGGAGTAAAATAAGG', 'end': 212566949, 'start': 212566621}, {'alt': None, 'ref': 'ATAGCGCAACAGTTGCAGTTTAAAAAATTACCTGTTATCTCTCTGACTGTCCGAAAGACGTTCAGTTTCTCTGGGTCTATGGCTTCAATTGCATTGTAAGGGTCCCTAGAAAATCAAGAAGAGATGTAGCCAAATTTAAATTTTACTAAAGGATTGAAAATATGAGAAATGTGACaatattttattcaaatttaattttaattaGCAGTAGGGGTGTTATCAATGCA', 'end': 212570238, 'start': 212570012}], 'chr3': [{'alt': None, 'ref': 'ATTTGCTTTTTCTGTAAATCATCTGTGAATCCAGAGGGGAAAAATATGACAAAGAAAGCTATATAAGATATTATTTTATTTTACAGAGTAACAGACTAGCTAGAGACAATGAATTAAGGGAAAATGACAAAGAACAGCTCAAAGCAATTTCTACACGAGATCCTCTCTCTGAAATCACTGAGCAGGAGAAAGATTTTCTATGGAGTCACAGGTAAGTGCTAAAATGGAGATTCTCTGTTTCTTTTTCTTTATTACAGAAAAAATAACTGAATTTGGCTGATCTCAGCATGTTTTTACCATACC', 'end': 178936214, 'start': 178935912}, {'alt': None, 'ref': 'TTATAGCTTTGTCTACGAAAGCCTCTCTAATTTTGTGACATTTGAGCAAAGACCTGAAGGTATTAACATCATTTGCTCCAAACTGACCAAACTGTTCTTATTACTTATAGGTTTCAGGAGATGTGTTACAAGGCTTATCTAGCTATTCGACAGCATGCCAATCTCTTCATAAATCTTTTCTCAATGATGCTTGGCTCTGGAATGCCAGAACTACAATCTTTTGATGACATTGCATACATTCGAAAGACCCTAGCCTTAGATAAAACTGAGCAAGAGGCTTTGGAGTATTTCATGAAACAAATGAATGATGCACATCATGGTGGCTGGACAACAAAAATGGATTGGATCTTCCACACAATTAAACAGCATGCATTGAACTGAAAAGATAACTGAGAAAATGAAAGCTCACTCTGGATTCCACACTGCACTGTTAATAACTCTCAGCAGGCAAAGACCGATTGCATAGGAATTGCACAATCCATGAACAGCATTAGAATTTACAGCAAGAACAGaaataaaatactatataatttaaataatGTAAACGCAAACAGGGTTTGATAGCACTTAAACTAGTTCATTTCAAAATTAAGCTTTAGAATAATGCGCAATTTCATGTTATGCCTTAAGTCCAAAAAGGTAAACTTTGAAGATTGTTTGTATCTTTTTTTAAAAAACAAAACAAAACAAAAATCCCCAAAATATATAGAAATGATGGAGAAGGAAAAAGTGATGGTTTTTTTTGTCTTGCAAATGTTCTATGTTTTGAAATGTGGACACAACAAAGGCTGTTATTG', 'end': 178952568, 'start': 178951772}], 'chr11': [{'alt': None, 'ref': 'AGCTAGAGCCAGAGCGGCTGCCCTGTGTCAAGGGAGAGGGTCAGTGAGTGCTGCTCCCTGGCTGGGGCGGGGCGGGGCGGGTCCCTGGCTAGCTGTGGGGTGGAGAGCTGCCTCACCTGCCGGGTCTTGGCCGAGGTCTCGATGTAGGGGATGCCGTAGCTTCGGGCGAGGTCCTGAGCCTGCCGAGATTCCACAGTGCGTGCAGCCAGGTCACACTTGTTCCCCACCAGCACCATGGGCACGTCATCCGAGTCCTTCACCCGTTTGATCTGCTCCCTGAGAGGTGGAAAGCGAGAGCTGGCTACGGGGGCTGCAGGCGCAGCGGCATCCAGGACATGCGCAGAGAGGACAGGAGGCCCCTGCCTGGACGCAGCCGGCCTGGCCCCACCTGTGCGGCGTGGGCTCCCGGGCCAGCCTCACGGGGTTCACCTGTACTGGTGGATGTCCTCAAAAGACTTGGTGTTGTTGATGGCAAACACACACAGGAAGCCCTCCCCGGTGCGCATGTACTGGTCCCGCATGGCGCTGTACTCCTCCTGGCCGGCGGTATCCAGGATGTCCAACAGGCACGTCTCCCCATCAATGACCACCTGCTTCCGGTAGGAATCCTGCAGGAGGACAGGGCTCAGGGACCCCCTCAGGACCTTCCGTGGGGGGAGTTCACACAGCCAGCCTCTCCCTGGTACCTCTCATGCCCCTCAT', 'end': 534038, 'start': 533337}, {'alt': None, 'ref': 'GCTGGCACCTGGACGGCGGCGCCAGGCTCACCTCTATAGTGGGGTCGTATTCGTCCACAAAATGGTTCTGGATCAGCTGGATGGTCAGCGCACTCTTGCCCACACCGCCGGCGCCCACCACCACCAGCTTATATTCCGTCATCGCTCCTCAGGGGCCTGCGGCCCGGGGTCCTCCTACAGGGTCTCCTGCCCCACCTGCCAAGGAGGGCCCTGCTCAGCCAGGCCCAGGCCCAGCCCCAGGCCCCACAGGGCAGCTGCTGGCAGGGCCATCTGAAGGGCAAA', 'end': 534462, 'start': 534181}], 'chr17': [{'alt': None, 'ref': 'TCCCATACCCTCTCAGCGTACCCTTGTCCCCAGGAAGCATACGTGATGGCTGGTGTGGGCTCCCCATATGTCTCCCGCCTTCTGGGCATCTGCCTGACATCCACGGTGCAGCTGGTGACACAGCTTATGCCCTATGGCTGCCTCTTAGACCATGTCCGGGAAAACCGCGGACGCCTGGGCTCCCAGGACCTGCTGAACTGGTGTATGCAGATTGCCAAGGTATGCACCTGGGCTCTTTGCAGGTCTCTCCGGAGCAAACCCCTATGTCCACAAGGGGCTAGGATGGGGACTCTTGCTGGGCATGT', 'end': 37881250, 'start': 37880946}], 'chr7': [{'alt': None, 'ref': 'ACCCTTGTCTCTGTGTTCTTGTCCCCCCCAGCTTGTGGAGCCTCTTACACCCAGTGGAGAAGCTCCCAACCAAGCTCTCTTGAGGATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGTAAGGTCCCTGGCACAGGCCTCTGGGCTGGGCCGCAGGGCCTCTCATGGTCTGGTGGGGAGCCCAGA', 'end': 55241804, 'start': 55241583}, {'alt': None, 'ref': 'CCAGTTAACGTCTTCCTTCTCTCTCTGTCATAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTCGATGTGAGTTTCTGCTTTGCTGTGTGGGGGTCCATGGCTCTGAACCTCAGGCCCACCTTTTCTCAT', 'end': 55242576, 'start': 55242382}, {'alt': None, 'ref': 'CACACTGACGTGCCTCTCCCTCCCTCCAGGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGTAATCAGGGAAGGGAGATACGGGGAGGGGAGATAAGGAGCCAGGATCCTCACATGCGGTCTGCGCTCCTGGGATAGCAAGAGTTTGCCATGGGGATATG', 'end': 55249271, 'start': 55248957}, {'alt': None, 'ref': 'CTGTCCCTCACAGCAGGGTCTTCTCTGTTTCAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTAAGGAGGTGGCTTTAGGTCAGCCAGCATTTTCCTGACACCAGGGACCAGGCTGCCTTCCCACTAGCTGTATTGTTTAACACATGCAGGGG', 'end': 55259659, 'start': 55259379}, {'alt': None, 'ref': 'AGTTTTTGTTTTCCTTCTATGTAAAAGTCCAGTTGGGAAGCTTTATTTCTGATAGATTAAATGGTATAGGTCTTTCAGTTTTCTCTTCATTTCTGACAACTGAACTGCTCTCGCCTTGAACCTGTTTTGGCAGATAAACCTCTCATAATGAAGGCCCCCGCTGTGCTTGCACCTGGCATCCTCGTGCTCCTGTTTACCTTGGTGCAGAGGAGCAATGGGGAGTGTAAAGAGGCACTAGCAAAGTCCGAGATGAATGTGAATATGAAGTATCAGCTTCCCAACTTCACCGCGGAAACACCCATCCAGAATGTCATTCTACATGAGCATCACATTTTCCTTGGTGCCACTAACTACATTTATGTTTTAAATGAGGAAGACCTTCAGAAGGTTGCTGAGTACAAGACTGGGCCTGTGCTGGAACACCCAGATTGTTTCCCATGTCAGGACTGCAGCAGCAAAGCCAATTTATCAGGAGGTGTTTGGAAAGATAACATCAACATGGCTCTAGTTGTCGACACCTACTATGATGATCAACTCATTAGCTGTGGCAGCGTCAACAGAGGGACCTGCCAGCGACATGTCTTTCCCCACAATCATACTGCTGACATACAGTCGGAGGTTCACTGCATATTCTCCCCACAGATAGAAGAGCCCAGCCAGTGTCCTGACTGTGTGGTGAGCGCCCTGGGAGCCAAAGTCCTTTCATCTGTAAAGGACCGGTTCATCAACTTCTTTGTAGGCAATACCATAAATTCTTCTTATTTCCCAGATCATCCATTGCATTCGATATCAGTGAGAAGGCTAAAGGAAACGAAAGATGGTTTTATGTTTTTGACGGACCAGTCCTACATTGATGTTTTACCTGAGTTCAGAGATTCTTACCCCATTAAGTATGTCCATGCCTTTGAAAGCAACAATTTTATTTACTTCTTGACGGTCCAAAGGGAAACTCTAGATGCTCAGACTTTTCACACAAGAATAATCAGGTTCTGTTCCATAAACTCTGGATTGCATTCCTACATGGAAATGCCTCTGGAGTGTATTCTCACAGAAAAGAGAAAAAAGAGATCCACAAAGAAGGAAGTGTTTAATATACTTCAGGCTGCGTATGTCAGCAAGCCTGGGGCCCAGCTTGCTAGACAAATAGGAGCCAGCCTGAATGATGACATTCTTTTCGGGGTGTTCGCACAAAGCAAGCCAGATTCTGCCGAACCAATGGATCGATCTGCCATGTGTGCATTCCCTATCAAATATGTCAACGACTTCTTCAACAAGATCGTCAACAAAAACAATGTGAGATGTCTCCAGCATTTTTACGGACCCAATCATGAGCACTGCTTTAATAGGGTAAGTCACATCAGTTCCCCACTTATAAACTGTGAGGTATAAATTAGAAATAAGTATCAGTCTCAAAAAGAATATCCAGGGCTTCTTTTGTGCTTTGTAAATGGTG', 'end': 116340444, 'start': 116338992}, {'alt': None, 'ref': 'TTTAACAAGCTCTTTCTTTCTCTCTGTTTTAAGATCTGGGCAGTGAATTAGTTCGCTACGATGCAAGAGTACACACTCCTCATTTGGATAGGCTTGTAAGTGCCCGAAGTGTAAGCCCAACTACAGAAATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGGTATATTTCAGTTTATTGTTCTGAGAAATACCTATACATATACCTCAGTGGGTTGTGACATTGTTGTTTATTTTTGG', 'end': 116412120, 'start': 116411870}, {'alt': None, 'ref': 'CCATTAAATGAGGTTTTACTGTTGTTCTTTAATAATTTTCCTTCATCTTACAGATCAGTTTCCTAATTCATCTCAGAACGGTTCATGCCGACAAGTGCAGTATCCTCTGACAGACATGTCCCCCATCCTAACTAGTGGGGACTCTGATATATCCAGTCCATTACTGCAAAATACTGTCCACATTGACCTCAGTGCTCTAAATCCAGAGCTGGTCCAGGCAGTGCAGCATGTAGTGATTGGGCCCAGTAGCCTGATTGTGCATTTCAATGAAGTCATAGGAAGAGGTAAGTATTTCCACTCAGCTTTTTGTTAAATACGATTTTCCAGTAAGCATTTTATCTTTGGCCTTTGCAGATTAGGAACTTAGACAATGGTGAAAGCAACTGACAGAGCAGTGATAACAAGTG', 'end': 116415288, 'start': 116414882}, {'alt': None, 'ref': 'CGCAGTGCTAACCAAGTTCTTTCTTTTGCACAGGGCATTTTGGTTGTGTATATCATGGGACTTTGTTGGACAATGATGGCAAGAAAATTCACTGTGCTGTGAAATCCTTGAACAGTAAGTGGCATTTTATTTAACCATGGAGTATACTTTTGTGGTTTGCAACCTAATAAATAGCTTATAATAAAACGTTGATTTACACTTTCCCCTTGTG', 'end': 116417620, 'start': 116417410}, {'alt': None, 'ref': 'AATGAAGTTAATGTCTCCACCACTGGATTTCTCAGGAATCACTGACATAGGAGAAGTTTCCCAATTTCTGACCGAGGGAATCATCATGAAAGATTTTAGTCATCCCAATGTCCTCTCGCTCCTGGGAATCTGCCTGCGAAGTGAAGGGTCTCCGCTGGTGGTCCTACCATACATGAAACATGGAGATCTTCGAAATTTCATTCGAAATGAGACTCATGTAAGTTGACTGCCAAGCTTACTAACTGGCAAACTAGCTGTAAGCCAGCCATCCCTTCAAAA', 'end': 116419073, 'start': 116418795}, {'alt': None, 'ref': 'GCTTAGTTTATGCTTTTCTAACTCTCTTTGACTGCAGAATCCAACTGTAAAAGATCTTATTGGCTTTGGTCTTCAAGTAGCCAAAGGCATGAAATATCTTGCAAGCAAAAAGTTTGTCCACAGAGACTTGGCTGCAAGAAACTGTATGTAAGTATCAGAATCTCTGTGCCACAATCCAAATTAAGTGACAAGGAGGA', 'end': 116422201, 'start': 116422005}, {'alt': None, 'ref': 'TCAGCCACGGGTAATAATTTTTGTCCTTTCTGTAGGCTGGATGAAAAATTCACAGTCAAGGTTGCTGATTTTGGTCTTGCCAGAGACATGTATGATAAAGAATACTATAGTGTACACAACAAAACAGGTGCAAAGCTGCCAGTGAAGTGGATGGCTTTGGAAAGTCTGCAAACTCAAAAGTTTACCACCAAGTCAGATGTGGTAATGTATTGGTTATCTCTGAGTTTCTCCTCTTTTACTTTCATATCCAACTTTTTTTGAAGTTTTATCACTACTTAATTTTTTAAAAAAATTCAACACCACCAATTCCAGTTTTC', 'end': 116423639, 'start': 116423323}, {'alt': None, 'ref': 'TCACCTCATCTGTCCTGTTTCTTGTTTTACTAGTGGTCCTTTGGCGTGCTCCTCTGGGAGCTGATGACAAGAGGAGCCCCACCTTATCCTGACGTAAACACCTTTGATATAACTGTTTACTTGTTGCAAGGGAGAAGACTCCTACAACCCGAATACTGCCCAGACCCCTTGTAAGTAGTCTTTCTGTACCTCTTACGTTCTTTACTTTTACAGAAATGCCTGCCTTCAAAGGG', 'end': 116435908, 'start': 116435676}, {'alt': None, 'ref': 'AATTTAATCAGTGGAAAAATAGCCTCAATTCTTACCATCCACAAAATGGATCCAGACAACTGTTCAAACTGATGGGACCCACTCCATCGAGATTTCACTGTAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATATCTGAGGTGTAGTAAGTAAAGGAAAACAGTAGATCTCATTTTCCTATCAGAGCAAGCA', 'end': 140453250, 'start': 140453040}, {'alt': None, 'ref': 'ACTTGTCACAATGTCACCACATTACATACTTACCATGCCACTTTCCCTTGTAGACTGTTCCAAATGATCCAGATCCAATTCTTTGTCCCACTGTAATCTGCCCATCAGGAATCTCCCAATCATCACTCGAGTCCCGTCTACCAAGTGTTTTCTTGATAAAAACAGTAAAAAAGTCAAGTCAAGCCAAACAGAAAAAGAAAACCTTATGTTTCACCCTAAGTACATTACCTTATGCCTG', 'end': 140481580, 'start': 140481343}], 'chr4': [{'alt': None, 'ref': 'tggCCCCTGAGCGTCATCTGCCCCCACAGAGCGCTCCCCGCACCGGCCCATCCTGCAGGCGGGGCTGCCGGCCAACCAGACGGCGGTGCTGGGCAGCGACGTGGAGTTCCACTGCAAGGTGTACAGTGACGCACAGCCCCACATCCAGTGGCTCAAGCACGTGGAGGTGAATGGCAGCAAGGTGGGCCCGGACGGCACACCCTACGTTACCGTGCTCAAGGTGGGCCACCGTGTGCACGTGGGTGCCGCCGCTGGGGCTCCTGGGCTGGCCCCAAGGGTGCCCCTTGGCTGCGGGTTGCGTGAGGATTTGGGTCTAGGGGTTGGAGCTT', 'end': 1803861, 'start': 1803533}, {'alt': None, 'ref': 'TGCGGCTCTGGGCCAGGGGCATCCATGGGAGCCCCGTGGGGGGGGGGGCCAGGCCAGGCCTCAACGCCCATGTCTTTGCAGCCGAGGAGGAGCTGGTGGAGGCTGACGAGGCGGGCAGTGTGTATGCAGGCATCCTCAGCTACGGGGTGGGCTTCTTCCTGTTCATCCTGGTGGTGGCGGCTGTGACGCTCTGCCGCCTGCGCAGCCCCCCCAAGAAAGGCCTGGGCTCCCCCACCGTGCACAAGATCTCCCGCTTCCCGCTCAAGCGACAGGTAACAGAAAGTAGATACCAGGTTCTGAGCTGCCTGCCCGCCAGGCCTCCTGGAGCCCCACCTCGGCCCACGCTGGTCCTGGGCTGTGTGAGCCCTCTCT', 'end': 1806347, 'start': 1805976}, {'alt': None, 'ref': 'GCTTCAGCCCTGCCTCCCACCCCTTCCCCAGTGCATCCACAGGGACCTGGCTGCCCGCAATGTGCTGGTGACCGAGGACAACGTGATGAAGATCGCAGACTTCGGGCTGGCCCGGGACGTGCACAACCTCGACTACTACAAGAAGACGACCAACGTGAGCCCGGCCCTGGGGTGCGGGGGTGGGGGTCATGCCAGTAGGACGCCTGGCGCCAACA', 'end': 1807961, 'start': 1807747}, {'alt': None, 'ref': 'AATTCACCAGTTACCTGTCCTGGTCATTTATAGAAACCGAGGTATGAAATTCGCTGGAGGGTCATTGAATCAATCAGCCCAGATGGACATGAATATATTTATGTGGACCCGATGCAGCTGCCTTATGACTCAAGATGGGAGTTTCCAAGAGATGGACTAGTGCTTGGTAAGTTCCATGGGGTAACCTCCCAAGACTCCCTTTTCCCTTGCACACAACTTTACAAT', 'end': 55141199, 'start': 55140975}, {'alt': None, 'ref': 'CAGCTGGACTGATATGTGATTTATTCTTTCAACAGCCACGGCCAGATCCAGTGAAAAACAAGCTCTCATGTCTGAACTGAAGATAATGACTCACCTGGGGCCACATTTGAACATTGTAAACTTGCTGGGAGCCTGCACCAAGTCAGGTGGGCTCACTGACCTGGAGTGAGGATTTTCACTGGACACATGTGGTTGTG', 'end': 55144224, 'start': 55144028}, {'alt': None, 'ref': 'TGATCCTGAGTCATTTCTTCCTTTTCCATGCAGTGTGTCCACCGTGATCTGGCTGCTCGCAACGTCCTCCTGGCACAAGGAAAAATTGTGAAGATCTGTGACTTTGGCCTGGCCAGAGACATCATGCATGATTCGAACTATGTGTCGAAAGGCAGTGTACGTCCTCACTTCCCTCACTGGTCAGGCTCATCCTCCTTCACTTTAATC', 'end': 55152181, 'start': 55151975}, {'alt': None, 'ref': 'CACTCTGACATATGGCCATTTCTGTTTTCCTGTAGCAAAACCAGAAATCCTGACTTACGACAGGCTCGTGAATGGCATGCTCCAATGTGTGGCAGCAGGATTCCCAGAGCCCACAATAGATTGGTATTTTTGTCCAGGAACTGAGCAGAGGTGAGATGATTATTTTTGGCACTGCTTATAATGCAGAGGGGAAGGACTGCAATTCACTTGAATTTCA', 'end': 55589931, 'start': 55589715}, {'alt': None, 'ref': 'AGTAAGCCAGGGCTTTTGTTTTCTTCCCTTTAGATGCTCTGCTTCTGTACTGCCAGTGGATGTGCAGACACTAAACTCATCTGGGCCACCGTTTGGAAAGCTAGTGGTTCAGAGTTCTATAGATTCTAGTGCATTCAAGCACAATGGCACGGTTGAATGTAAGGCTTACAACGATGTGGGCAAGACTTCTGCCTATTTTAACTTTGCATTTAAAGGTAACAACAAAGGTATATTTCTTTTTAATCCAATTTAAGGGGATGTTTAGGCTCTGTCTACCATATCAGTCATGAT', 'end': 55592280, 'start': 55591990}, {'alt': None, 'ref': 'ATTAAAAGGTGATCTATTTTTCCCTTTCTCCCCACAGAAACCCATGTATGAAGTACAGTGGAAGGTTGTTGAGGAGATAAATGGAAACAATTATGTTTACATAGACCCAACACAACTTCCTTATGATCACAAATGGGAGTTTCCCAGAAACAGGCTGAGTTTTGGTCAGTATGAAACAGGGGCTTTCCATGTCACCTTTTTGGGTACACATAACAGTGACT', 'end': 55593765, 'start': 55593545}, {'alt': None, 'ref': 'GCTTTTTGCTAAAATGCATGTTTCCAATTTTAGCGAGTGCCCATTTGACAGAACGGGAAGCCCTCATGTCTGAACTCAAAGTCCTGAGTTACCTTGGTAATCACATGAATATTGTGAATCTACTTGGAGCCTGCACCATTGGAGGTAAAGCCGTGTCCAAGCTGCCTTTTATTGTCTGTCAGGTTATCAAAACATGACA', 'end': 55594342, 'start': 55594144}, {'alt': None, 'ref': 'TAAATGGTTTTCTTTTCTCCTCCAACCTAATAGTGTATTCACAGAGACTTGGCAGCCAGAAATATCCTCCTTACTCATGGTCGGATCACAAAGATTTGTGATTTTGGTCTAGCCAGAGACATCAAGAATGATTCTAATTATGTGGTTAAAGGAAACGTGAGTACCCATTCTCTGCTTGACAGTCCTGCAAAGGATTTTTAG', 'end': 55599403, 'start': 55599203}, {'alt': None, 'ref': 'AACATTATTGACTCTGTTGTGCTTCTATTACAGGCTCGACTACCTGTGAAGTGGATGGCACCTGAAAGCATTTTCAACTGTGTATACACGTTTGAAAGTGACGTCTGGTCCTATGGGATTTTTCTTTGGGAGCTGTTCTCTTTAGGTAAAATGATCCTTGCCAAAGACAACTTCATTAGACTCAGAGCATCTTCTTGAAGTTTCATTGGTGTCCTGC', 'end': 55602847, 'start': 55602631}]}
    print( "[GET sequences]\tEnd" )

    # Get variant profile
    models = getVariantsProfile( args.input_profile )

    # Apply variant profile
    #################################
    random_seed = args.random_seed if args.random_seed is not None else int(time.time())
    random.seed(random_seed)
    print( "Random seed used:\t" + str(random_seed) )

    # Find positions
    print( "[GET variants pos]\tStart" )
    area_wt_primers_by_chr = getCoveredFromAmplicons( amplicons_by_chr, True )
    mutable_pos = getCoveredPos( area_wt_primers_by_chr )
    covered_len = len(mutable_pos)
    variant_by_pos = dict()
    for model in models:
        nb_variant = int(model["occurence"] * covered_len)
        for idx in range(nb_variant):
            mutable_pos_len = len(mutable_pos)
            variant_len = random.randint(model["length"]["min"], model["length"]["max"])
            ########################################### fct
            pos_start_idx = None
            pos_end_idx = None
            variant_start_chr = None
            variant_end_chr = None
            valid_position = False
            while not valid_position:
                pos_start_idx = random.randint(0, mutable_pos_len -1)
                pos_end_idx = pos_start_idx + variant_len -1
                variant_start_chr, variant_start_pos = mutable_pos[pos_start_idx].split(":")
                variant_start_pos = int(variant_start_pos)
                if pos_end_idx < mutable_pos_len:
                    variant_end_chr, variant_end_pos = mutable_pos[pos_end_idx].split(":")
                    variant_end_pos = int(variant_end_pos)
                    if variant_end_chr == variant_start_chr:
                        if variant_end_pos == (variant_start_pos + variant_len -1):
                            valid_position = True
            variant_chr = variant_start_chr
            ###########################################
            masked_indexes = [pos_idx for pos_idx in range(pos_start_idx, pos_end_idx + 1)] # List of indexes to remove from the mutable positions
            for margin_offset in range(1, (args.min_distance + 1)):
                if pos_end_idx + margin_offset < mutable_pos_len:
                    end_chr = mutable_pos[pos_end_idx + margin_offset].split(":")[0]
                    if variant_chr == end_chr:
                        masked_indexes.append( pos_end_idx + margin_offset )
                if pos_start_idx - margin_offset >= 0:
                    start_chr = mutable_pos[pos_start_idx - margin_offset].split(":")[0]
                    if variant_chr == start_chr:
                        masked_indexes.insert(0, pos_start_idx - margin_offset )
            # Remove mutated pos and margin from mutable positions
            masked_indexes = sorted(masked_indexes, reverse=True)
            if( len(masked_indexes) < (variant_len + 2*args.min_distance) ):
                print( "[GET variants pos]\t\t|- Limited mask on: " + ", ".join([mutable_pos[idx] for idx in masked_indexes]) )
            for idx in masked_indexes:
                del mutable_pos[idx]
            variant_freq = random.randint(int(1000*model["freq"]["min"]), int(1000*model["freq"]["max"]))/float(1000)
            if variant_chr not in variant_by_pos:
                variant_by_pos[variant_chr] = dict()
            variant_by_pos[variant_chr][variant_start_pos] = {"type":model["type"], "length":variant_len, "ref":None, "alt":None, "freq":variant_freq}
    print( "[GET variants pos]\tEnd" )

    # Create variants set origin
    for chr in area_by_chr:
        if chr in variant_by_pos:
            for area in area_by_chr[chr]:
                area["alt"] = [list(area["ref"]) for idx in range(100)]

    # Create variant set
    for chr in area_by_chr:
        if chr in variant_by_pos:
            for area in area_by_chr[chr]:
                for pos in range(area["start"], (area["end"] + 1)):
                    if pos in variant_by_pos[chr]:
                        # Set alternative variant
                        setAlt( variant_by_pos[chr][pos], area["ref"], pos - area["start"] )
                        ref_allele = variant_by_pos[chr][pos]["ref"]
                        alt_allele = variant_by_pos[chr][pos]["alt"]
                        #~ print( (pos - area["start"]), ref_allele, alt_allele )
                        # Random selection of alternative sequence
                        nb_alterated_amplicons = int( (variant_by_pos[chr][pos]["freq"]*100) + 0.5 )
                        alterated_amplicons = random.sample( area["alt"], nb_alterated_amplicons )
                        #~ print("-------------------------------------")
                        # Apply variant in selected alternatives
                        for amplicon_variant in alterated_amplicons:
                            #~ print( amplicon_variant )
                            # Apply variant
                            if "." in alt_allele or len(alt_allele) < len(ref_allele): # Deletion
                                del_length = len(ref_allele)
                                for idx in range(del_length):
                                    amplicon_variant[(pos + idx) - area["start"]] = ""
                            elif "." in ref_allele or len(alt_allele) > len(ref_allele): # Insertion
                                amplicon_variant[pos - area["start"]] = alt_allele 
                            else: # Substitution
                                subst_length = len(alt_allele)
                                for idx in range(subst_length):
                                    amplicon_variant[(pos + idx) - area["start"]] = alt_allele[idx]
                            #~ print( amplicon_variant )
                            #~ print("")
    
    # Write variant
    print( "[WRITE variants]\tStart" )
    writeVariants( args.output_variants, variant_by_pos )
    del variant_by_pos
    print( "[WRITE variants]\tEnd" )
    
    # Generate reads
    print( "[GENERATE reads]\tStart" )
    nb_pairs = {
        "min": 200,
        "max": 15000
    }
    reads_length = 150
    dual_lib = True
    idx_ampli = 1
    lib = {
        "A":{ 
            "R1": FastaIO(getFastaPath(args.output_folder, args.sample_name, "A", 1), "w"),
            "R2": FastaIO(getFastaPath(args.output_folder, args.sample_name, "A", 2), "w")
        },
        "B":{
            "R1": FastaIO(getFastaPath(args.output_folder, args.sample_name, "B", 1), "w"),
            "R2": FastaIO(getFastaPath(args.output_folder, args.sample_name, "B", 2), "w")
        }
    }
    
    for chr in amplicons_by_chr:
        for amplicon in amplicons_by_chr[chr]:
            # Select area ############################################## fct
            amplicon_area = None
            for area in area_by_chr[chr]:
                if amplicon["start"] >= area["start"]:
                    if amplicon["end"] <= area["end"]:
                        amplicon_area = area
            start_idx = amplicon["start"] - amplicon_area["start"]
            end_idx = amplicon["end"] - amplicon_area["start"]
            # Primers
            up_primer = amplicon["up_primer"]
            down_primer = amplicon["down_primer"]
            if amplicon["strand"] == "-":
                up_primer = revcom(amplicon["down_primer"])
                down_primer = revcom(amplicon["up_primer"])
            interest_start_idx = start_idx + len(up_primer)
            interest_end_idx = end_idx - len(down_primer)
            # Reads
            nb_selected_pairs = random.randrange( nb_pairs["min"], nb_pairs["max"], 100 )
            nb_from_alt = int(nb_selected_pairs/100)
            if amplicon_area["alt"] is None:
                amplicon_area["alt"] = [list(amplicon_area["ref"])]
                nb_from_alt = nb_selected_pairs
            for alt_seq in amplicon_area["alt"]:
                matrix = up_primer + "".join(alt_seq[interest_start_idx:interest_end_idx+1]) + down_primer # The protocol primers mask variation on matrix
                rvc_matrix = revcom(matrix)
                if amplicon["strand"] == "-":
                    tmp = rvc_matrix
                    rvc_matrix = matrix
                    matrix = tmp    
                for idx in range(nb_from_alt):
                    pair_id = "theoritical:" + chr + "_" + str(amplicon["start"]) + "-" + str(amplicon["end"]) + ":" + str(idx_ampli).zfill(10)
                    writeReadsPair( pair_id, matrix, rvc_matrix, reads_length, lib["A"]["R1"], lib["A"]["R2"] )
                    idx_ampli += 1
                    if dual_lib:
                        pair_id = "theoritical:" + chr + "_" + str(amplicon["start"]) + "-" + str(amplicon["end"]) + ":" + str(idx_ampli).zfill(10)
                        writeReadsPair( pair_id, rvc_matrix, matrix, reads_length, lib["B"]["R1"], lib["B"]["R2"] )
                        idx_ampli += 1
    for lib_name in lib:
        for read in lib[lib_name]:
            if lib[lib_name][read] is not None:
                lib[lib_name][read].close()
    print( "[GENERATE reads]\tEnd" )
