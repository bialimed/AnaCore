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
__version__ = '1.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import os
import sys
import time
import random
import warnings
import logging
import argparse

CURRENT_DIR = os.path.dirname(__file__)
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
class Amplicon(object):
    def __init__(self, reference, strand, up_primer, down_primer, start, end, name=None):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.reference = reference
        self.up_primer = up_primer.upper()
        self.down_primer = down_primer.upper()

    def __str__(self):
        return str({
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "region": self.reference,
            "up_primer": self.up_primer,
            "down_primer": self.down_primer
        })

    def getInterestStart(self):
        interest_start = self.start + len(self.up_primer)
        if self.strand == "-":
            interest_start = self.start + len(self.down_primer)
        return( interest_start )

    def getInterestEnd(self):
        interest_end = self.end - len(self.down_primer)
        if self.strand == "-":
            interest_end = self.end - len(self.up_primer)
        return( interest_end )

def getOverlappedComplexRegions( complex_regions, variant_start, variant_end ):
    overlapped = list()
    for curr_region in complex_regions:
        if variant_start < curr_region.end and variant_end > curr_region.start:
            overlapped.append( curr_region )
    return overlapped

def hasComplexConfiguration( complex_regions, variant_start, variant_end ):
    is_complex = False
    overlapped = getOverlappedComplexRegions( complex_regions, variant_start, variant_end )
    if len(overlapped) > 0:
        for curr_region in overlapped:
            if variant_start <= curr_region.getInterestEnd() and variant_end > curr_region.getInterestEnd():
                is_complex = True
            elif variant_start < curr_region.getInterestStart() and variant_end >= curr_region.getInterestStart():
                is_complex = True
    return is_complex

def getComplexRegions( amplicons_by_chr ):
    complex_by_chr = dict()
    for chrom in amplicons_by_chr:
        complex_by_chr[chrom] = list()
        amplicons = sorted(amplicons_by_chr[chrom], key=lambda x: (x.start, x.end))
        prev_ampli = None
        for curr_ampli in amplicons:
            if prev_ampli is not None:
                if prev_ampli.getInterestEnd() + 1 >= curr_ampli.getInterestStart(): # It does not exist a gap in coverage region
                    complex_by_chr[chrom].append(
                        Amplicon(
                            chrom,
                            "+",
                            (curr_ampli.up_primer if curr_ampli.strand == "+" else revcom(curr_ampli.down_primer)),
                            (prev_ampli.down_primer if prev_ampli.strand == "+" else revcom(prev_ampli.up_primer)),
                            curr_ampli.start,
                            prev_ampli.end
                        )
                    )
            prev_ampli = curr_ampli
    return complex_by_chr

def getAmplicons( reference_path, manifest_path ):
    """
    @summary: Returns the list of amplicons by chromosome from a Illumina's manifest.
    @param reference_path: [str] Path to the genome assembly where the amplicon have been defined.
    @param manifest_path: [str] Path to the manifest.
    @return: [dict] By chr the amplicons objects.
    """
    # Get amplicons by chr
    amplicons = getAmpliconsFromManifest( manifest_path )
    amplicons_by_chr = dict()
    for ampli_idx, ampli in enumerate(amplicons):
        chr = ampli["chromosome"]
        if chr not in amplicons_by_chr:
            amplicons_by_chr[chr] = list()
        amplicons_by_chr[chr].append(
            Amplicon(
                chr,
                ampli["probe_strand"],
                ampli["ulso_sequence"],
                ampli["dlso_sequence"],
                None,
                None,
                "ampl" + str(ampli_idx) + "_" + ampli["target_region_name"].replace(" ", "_")
            )
        )

    # Find amplicons coord
    FH_ref = SequenceFileReader.factory( reference_path )
    try:
        for record in FH_ref:
            if record.id in amplicons_by_chr:
                chr_str = record.string.upper()
                for ampli in amplicons_by_chr[record.id]:
                    up_primer = ampli.up_primer
                    down_primer = ampli.down_primer
                    if ampli.strand == "-":
                        up_primer = revcom(ampli.down_primer)
                        down_primer = revcom(ampli.up_primer)
                    # Find positions on chr
                    upstream_matches = list()
                    up_pattern = re.compile(up_primer)
                    for curr_match in up_pattern.finditer(chr_str):
                        upstream_matches.append({"start": curr_match.start()+1, "end": curr_match.end()})
                    downstream_matches = list()
                    down_pattern = re.compile(down_primer)
                    for curr_match in down_pattern.finditer(chr_str):
                        downstream_matches.append({"start": curr_match.start()+1, "end": curr_match.end()})
                    if len(upstream_matches) == 0 or len(downstream_matches) == 0:
                        raise Exception("The primers '" + up_primer + "' and '" + down_primer + "' cannot be found in " + record.id)
                    # Check multiple target in chr
                    if len(upstream_matches) > 1:
                        match_list = ", ".join(["{}:{}-{}".format(record.id, curr_match["start"], curr_match["end"]) for curr_match in upstream_matches])
                        warnings.warn("The primer '" + up_primer + "' is found multiple twice in " + record.id + " (" + match_list + ")")
                    if len(downstream_matches) > 1:
                        match_list = ", ".join(["{}:{}-{}".format(record.id, curr_match["start"], curr_match["end"]) for curr_match in downstream_matches])
                        warnings.warn("The primer '" + down_primer + "' is found multiple twice in " + record.id + " (" + match_list + ")")
                    # Select smaller amplified fragment
                    prev_length = None
                    for curr_up in upstream_matches:
                        for curr_down in downstream_matches:
                            curr_length = curr_down["start"] - curr_up["end"]
                            if curr_length >= 0:
                                if prev_length is None or prev_length > curr_length:
                                    prev_length = curr_length
                                    ampli.start = curr_up["start"]
                                    ampli.end = curr_down["end"]
                amplicons_by_chr[record.id] = sorted( amplicons_by_chr[record.id], key=lambda ampl: (ampl.start, ampl.end) )
    finally:
        FH_ref.close()

    return( amplicons_by_chr )

def revcom( seq ):
    """
    @summary: Returns the reverse complement the sequence.
    @param seq: [str] The sequence.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}

    return( "".join([complement_rules[base] for base in seq[::-1]]) )

def getAmpliconsFromManifest( manifest_path ):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons information.
    """
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

def getCoveredLength( areas_by_chr ):
    covered_len = 0
    for chr in areas_by_chr:
        for area in areas_by_chr[chr]:
            covered_len += area["end"] - area["start"] + 1
    return( covered_len )

def setAmpliconSeq( ref_sequences, areas_by_chr ):
    FH_ref = SequenceFileReader.factory( ref_sequences )
    for record in FH_ref:
        if record.id in areas_by_chr:
            for curr_area in areas_by_chr[record.id]:
                curr_area["ref"] = record.string[(curr_area["start"]-1):curr_area["end"]]
    FH_ref.close()

def getCoveredPos( areas_by_chr ):
    covered_pos = list()
    for chrom in sorted(areas_by_chr):
        for area in areas_by_chr[chrom]:
            for pos in range( area["start"], (area["end"] + 1) ):
                covered_pos.append( chrom + ":" + str(pos) )
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
    areas_by_chr = dict()
    for chrom in sorted(amplicons_by_chr):
        areas_by_chr[chrom] = list()
        amplicons = sorted(amplicons_by_chr[chrom], key=lambda x: (x.start, x.end))
        prev_start = amplicons[0].getInterestStart() if remove_primers else amplicons[0].start
        prev_end = amplicons[0].getInterestEnd() if remove_primers else amplicons[0].end
        for curr_ampl in amplicons[1:]:
            curr_start = curr_ampl.getInterestStart() if remove_primers else curr_ampl.start
            curr_end = curr_ampl.getInterestEnd() if remove_primers else curr_ampl.end
            if curr_start > prev_end + 1: # no overlap
                areas_by_chr[chrom].append({"start": prev_start, "end": prev_end, "ref":None, "alt":None})
                prev_start = curr_start
                prev_end = curr_end
            else: # overlap or inclusion
                prev_end = max(prev_end , curr_end)
        areas_by_chr[chrom].append({"start": prev_start, "end": prev_end, "ref":None, "alt":None})
    return areas_by_chr

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

def minPairsType(value):
    ivalue = int(value)
    if ivalue < 100:
        raise argparse.ArgumentTypeError( "The minimum number of reads pairs by amplicon must be 100." )
    return ivalue

def maxPairsType(value):
    ivalue = int(value)
    if ivalue % 100 != 0:
        raise argparse.ArgumentTypeError( "The maximum number of reads pairs by amplicon must be a multiple of 100." )
    return ivalue


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Creates simulated reads corresponding Illumina sequencing of amplicon double strand experiment. ***variants *********without error *************************' )
    parser.add_argument( '-s', '--random-seed', type=int, help="The seed used for the random generator. If you want reproduce results of one execution: use the same parameters AND the same random-seed. [Default: auto]" )
    parser.add_argument( '-d', '--min-distance', type=int, default=2, help="The minimum distance between two variants. [Default: %(default)s]" )
    parser.add_argument( '-n', '--sample-name', required=True, help="The name of the simulated sample." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_reads = parser.add_argument_group( 'Reads' ) # Reads
    group_reads.add_argument( '-l', '--reads-length', type=int, default=150, help='The reads length for sequencer. [Default: %(default)s]' )
    group_reads.add_argument( '-mip', '--min-pairs', type=minPairsType, default=200, help='The minimum number of reads pairs for one amplicon in one library. [Default: %(default)s]' )
    group_reads.add_argument( '-map', '--max-pairs', type=maxPairsType, default=15000, help='The maximum number of reads pairs for one amplicon in one library. [Default: %(default)s]' )
    group_reads.add_argument( '-ds', '--is-double-strand', action='store_true', help='With this argument the sequence are created for the specified design and for his reverse complement design (library A and library B in double strand experiment).' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--input-manifests', required=True, help='Path to the definition of the amplicons (format: Illumina manifest).' )
    group_input.add_argument( '-g', '--input-genome', required=True, help='Path to the genome sequences file (format: fasta).' )
    group_input.add_argument( '-p', '--input-profile', required=True, help='Path to the variants profile file (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--output-variants', required=True, help='Path to the file describing variants introduced in reads and their AF (format: VCF).' )
    group_output.add_argument( '--output-folder', default=os.getcwd(), help='Path to the output folder. [Default: %(default)s]' )
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s -- [%(name)s][pid:%(process)d][%(levelname)s] -- %(message)s' )
    logger = logging.getLogger( os.path.basename(__file__) )
    logging.info( " ".join(sys.argv) )

    # Load variant profile
    models = getVariantsProfile( args.input_profile )

    # Get reference sequences
    logging.info( "\t[GET amplicons] START" )
    amplicons_by_chr = getAmplicons( args.input_genome, args.input_manifests )
    logging.info( "\t[GET amplicons] END" ) ############################### Sortir ce code avec la recup des sequences dans un autre script

    # Get covered positions
    logging.info( "\t[GET covered areas] START" )
    areas_by_chr = getCoveredFromAmplicons( amplicons_by_chr, False )
    areas_wout_primers_by_chr = getCoveredFromAmplicons( amplicons_by_chr, True )
    logging.info( "\t[GET covered areas] END" )

    # Find warn positions
    logging.info( "\t[GET warn positions] START" )
    complex_areas = getComplexRegions( amplicons_by_chr ) # Positions where the variant can overlap primer and interest in multiple amplicons
    logging.info( "\t[GET warn positions] END" )


    # Find variants positions
    logging.info( "\t[GET variants pos] START" )
    random_seed = args.random_seed if args.random_seed is not None else int(time.time())
    random.seed(random_seed)
    logging.info( "\t[GET variants pos] Random seed used: " + str(random_seed)  )
    mutable_pos = getCoveredPos( areas_wout_primers_by_chr )
    del( areas_wout_primers_by_chr )
    covered_len = len(mutable_pos)
    variant_by_pos = dict()
    for model in models:
        nb_variants = int(model["occurence"] * covered_len)
        for idx in range(nb_variants):
            mutable_pos_len = len(mutable_pos)
            variant_len = random.randint(model["length"]["min"], model["length"]["max"])
            # Find mutated positions
            ########################################### fct
            pos_start_idx = None
            pos_end_idx = None
            variant_start_chr = None
            variant_end_chr = None
            valid_position = False
            while not valid_position:
                pos_start_idx = random.randint(0, mutable_pos_len -1)
                pos_end_idx = pos_start_idx + variant_len -1 if model["type"] != "insertion" else pos_start_idx
                variant_start_chr, variant_start_pos = mutable_pos[pos_start_idx].split(":")
                variant_start_pos = int(variant_start_pos)
                if pos_end_idx < mutable_pos_len: # not out of index
                    variant_end_chr, variant_end_pos = mutable_pos[pos_end_idx].split(":")
                    variant_end_pos = int(variant_end_pos)
                    if variant_len == 1 or model["type"] == "insertion":
                        valid_position = True
                    elif variant_end_chr == variant_start_chr: # start and end are on the same chromosome
                        if variant_end_pos == (variant_start_pos + variant_len -1): # The positions are continuous on chromosome
                            if not hasComplexConfiguration( complex_areas[variant_start_chr], variant_start_pos, variant_end_pos ):
                                valid_position = True
                            else:
                                logging.debug( "\t[GET variants pos] Simulation impossible for " + variant_start_chr + ":" + str(variant_start_pos) + "-" + str(variant_end_pos) + " the variant will be moved to another location." )
            variant_chr = variant_start_chr
            ###########################################
            # Remove mutated pos and margin from mutable positions
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
            masked_indexes = sorted(masked_indexes, reverse=True)
            clean_variant_len = variant_len if model["type"] != "insertion" else 1
            if( len(masked_indexes) < (clean_variant_len + 2*args.min_distance) ):
                logging.debug( "[GET variants pos] Limited mask on " + ", ".join([mutable_pos[idx] for idx in masked_indexes]) + " for " + model["type"] + " in position " + variant_chr + ":" + str(variant_start_pos) + "-" + str(variant_end_pos)  )
            for idx in masked_indexes:
                del mutable_pos[idx]
            # Store variant information
            variant_freq = random.randint(int(1000*model["freq"]["min"]), int(1000*model["freq"]["max"]))/float(1000)
            if variant_chr not in variant_by_pos:
                variant_by_pos[variant_chr] = dict()
            variant_by_pos[variant_chr][variant_start_pos] = {"type":model["type"], "length":variant_len, "ref":None, "alt":None, "freq":variant_freq}
    logging.info( "\t[GET variants pos] END" )

    # Create variants set origin
    logging.info( "\t[INIT covered sequences] START" )
    setAmpliconSeq( args.input_genome, areas_by_chr ) ############################################### Ne faire la recherche qu une fois pour les amplicon puis fusionner les region a ce moment du code
    for chrom in areas_by_chr:
        if chrom in variant_by_pos:
            for area in areas_by_chr[chrom]:
                area["alt"] = [list(area["ref"]) for idx in range(100)]
    logging.info( "\t[INIT covered sequences] END" )

    # Create variant set
    for chrom in sorted(areas_by_chr):
        if chrom in variant_by_pos:
            for area in areas_by_chr[chrom]:
                for pos in range(area["start"], (area["end"] + 1)):
                    if pos in variant_by_pos[chrom]:
                        # Set alternative variant
                        setAlt( variant_by_pos[chrom][pos], area["ref"], pos - area["start"] )
                        ref_allele = variant_by_pos[chrom][pos]["ref"]
                        alt_allele = variant_by_pos[chrom][pos]["alt"]
                        # Random selection of alternative sequence
                        nb_alterated_amplicons = int( (variant_by_pos[chrom][pos]["freq"]*100) + 0.5 )
                        alterated_amplicons = random.sample( area["alt"], nb_alterated_amplicons )
                        # Apply variant in selected alternatives
                        for amplicon_variant in alterated_amplicons:
                            # Apply variant
                            if "-" in alt_allele or len(alt_allele) < len(ref_allele): # Deletion
                                del_length = len(ref_allele)
                                for idx in range(del_length):
                                    amplicon_variant[(pos + idx) - area["start"]] = ""
                            elif "-" in ref_allele or len(alt_allele) > len(ref_allele): # Insertion
                                amplicon_variant[pos - area["start"]] = alt_allele
                            else: # Substitution
                                subst_length = len(alt_allele)
                                for idx in range(subst_length):
                                    amplicon_variant[(pos + idx) - area["start"]] = alt_allele[idx]
    # Write variant
    writeVariants( args.output_variants, variant_by_pos )
    del variant_by_pos

    # Generate reads
    logging.info( "\t[GENERATE reads] START" )
    lib = {
        "A": {
            "R1": FastaIO(getFastaPath(args.output_folder, args.sample_name, "A", 1), "w"),
            "R2": FastaIO(getFastaPath(args.output_folder, args.sample_name, "A", 2), "w")
        }
    }
    if args.is_double_strand:
        lib["B"] = {
            "R1": FastaIO(getFastaPath(args.output_folder, args.sample_name, "B", 1), "w"),
            "R2": FastaIO(getFastaPath(args.output_folder, args.sample_name, "B", 2), "w")
        }

    idx_ampli = 1
    for chrom in sorted(amplicons_by_chr):
        for amplicon in amplicons_by_chr[chrom]:
            # Select area ############################################## fct
            amplicon_area = None
            for area in areas_by_chr[chrom]:
                if amplicon.start >= area["start"]:
                    if amplicon.end <= area["end"]:
                        amplicon_area = area
            start_idx = amplicon.start - amplicon_area["start"]
            end_idx = amplicon.end - amplicon_area["start"]
            # Primers
            up_primer = amplicon.up_primer
            down_primer = amplicon.down_primer
            if amplicon.strand == "-":
                up_primer = revcom(amplicon.down_primer)
                down_primer = revcom(amplicon.up_primer)
            interest_start_idx = start_idx + len(up_primer)
            interest_end_idx = end_idx - len(down_primer)
            # Reads
            nb_selected_pairs = random.randrange( args.min_pairs, args.max_pairs, 100 )
            nb_from_alt = int(nb_selected_pairs/100)
            if amplicon_area["alt"] is None:
                amplicon_area["alt"] = [list(amplicon_area["ref"])]
                nb_from_alt = nb_selected_pairs
            for alt_seq in amplicon_area["alt"]:
                matrix = up_primer + "".join(alt_seq[interest_start_idx:interest_end_idx+1]) + down_primer # The protocol primers mask variation on matrix
                rvc_matrix = revcom(matrix)
                if amplicon.strand == "-":
                    tmp = rvc_matrix
                    rvc_matrix = matrix
                    matrix = tmp
                for idx in range(nb_from_alt):
                    pair_id = "theoritical:" + chrom + "_" + str(amplicon.start) + "-" + str(amplicon.end) + ":" + str(idx_ampli).zfill(10)
                    writeReadsPair( pair_id, matrix, rvc_matrix, args.reads_length, lib["A"]["R1"], lib["A"]["R2"] )
                    idx_ampli += 1
                    if args.is_double_strand:
                        pair_id = "theoritical:" + chrom + "_" + str(amplicon.start) + "-" + str(amplicon.end) + ":" + str(idx_ampli).zfill(10)
                        writeReadsPair( pair_id, rvc_matrix, matrix, args.reads_length, lib["B"]["R1"], lib["B"]["R2"] )
                        idx_ampli += 1
    for lib_name in lib:
        for read in lib[lib_name]:
            lib[lib_name][read].close()
    logging.info( "\t[GENERATE reads] END" )
