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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import sys
import pysam
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import *


########################################################################
#
# FUNCTIONS
#
########################################################################
def getADP( chrom, pos, ref, alt, aln_file, selected_RG=None ):
    """
    @summary: Returns the allele depth (AD) and the depth (DP) for the specified variant.
    @param chrom: [str] The variant region name.
    @param pos: [int] The variant position.
    @param ref: [str] The reference allele at the position.
    @param alt: [str] The variant allele at the position.
    @param aln_file: [str] The path to the alignment file (format: BAM). The AD and DP are retrieved from reads of this file. This file must be indexed.
    @param selected_RG: [list] The ID of RG used in AD and DP. Default: all read groups.
    @return: [list] The first element is the AD, the second is the DP.
    @warning: Reads ID must be unique in SAM.
    """
    ref_start = pos
    ref_end = pos + len(ref.replace(".", "")) - 1
    selected_RG = {RG:1 for RG in selected_RG}
    # Retrieve reads
    inspect_start = ref_start - 1
    inspect_end = ref_end
    reads = dict()
    with pysam.AlignmentFile( aln_file, "rb" ) as FH_sam:
        for pileupcolumn in FH_sam.pileup( chrom, inspect_start, inspect_end, max_depth=100000 ):
            for pileupread in pileupcolumn.pileups:
                if selected_RG is None or (pileupread.alignment.get_tag("RG") in selected_RG):
                    if pileupcolumn.pos >= inspect_start and pileupcolumn.pos < inspect_end:
                        # Get unique id
                        read_id = pileupread.alignment.query_name
                        if not pileupread.alignment.is_read2:
                            read_id += "_R1"
                        else:
                            read_id += "_R2"
                        # Store new reads
                        if read_id not in reads:
                            reads[read_id] = [None for pos in range(inspect_start, pileupcolumn.pos)]
                        # Store comparison with ref for current position
                        if pileupread.is_del: # Deletion
                            reads[read_id].append( "" )
                        elif pileupread.indel > 0: # Insertion
                            insert = ""
                            for insert_idx in range(pileupread.indel + 1):
                                insert += pileupread.alignment.query_sequence[pileupread.query_position + insert_idx].upper()
                            reads[read_id].append( insert )
                        elif not pileupread.is_refskip: # Substitution
                            reads[read_id].append(
                                pileupread.alignment.query_sequence[pileupread.query_position].upper()
                            )
    # Process AD and DP
    alt = alt if alt != "." else ""
    AD = 0
    DP = 0
    for read_id in reads:
        if None not in reads[read_id]: # Skip partial reads
            DP += 1
            if "".join(reads[read_id]) == alt:
                AD += 1
    # Return
    return( AD, DP )

def getNonOverlappingThread( selected_areas ):
    """
    @summary: Splits the list of regions in sub-lists of non-overlapping regions and returns them.
    @param selected_areas: [list] The selected areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    @return: [list] Each element of the list is a list of non-overlapping regions.
    """

    for curr_area in sorted_selected_areas[1:]:
        thread_idx = 0
        is_valid_thread = False
        while thread_idx < len(non_overlapping_threads) and not is_valid_thread:
            is_valid_thread = hasNoOverlap( curr_area, non_overlapping_threads[thread_idx][-1] )
            thread_idx += 1
        if is_valid_thread:
            non_overlapping_threads[thread_idx - 1].append(
                curr_area
            )
        else:
            non_overlapping_threads.append([
                curr_area
            ])
    return( non_overlapping_threads )

def getSelectedArea( input_panel ):
    """
    @summary: Returns the list of selected areas from a BED file.
    @param input_panel: [str] The path to the selected areas description (format: BED).
    @return: [list] The list of BED's areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608}.
    """
    selected_areas = list()
    with open(input_panel) as FH_panel:
        for line in FH_panel:
            if not line.startswith("browser ") and not line.startswith("track ") and not line.startswith("#"):
                fields = [elt.strip() for elt in line.split("\t")]
                selected_areas.append({
                    "region": fields[0],
                    "start": int(fields[1]) +1, # Start in BED is 0-based
                    "end": int(fields[2])
                })
    return( selected_areas )

def getOverlaps( selected_areas ):
    # Sort selected regions
    sorted_selected_areas = sorted(selected_areas, key=lambda x: (x["region"], x["start"], x["end"]))
    # Retrieve overlapping positions between area
    cover_pos_by_chr = dict() # overlapping positions between area
    prev_end = -1
    prev_region = None
    for curr_area in sorted_selected_areas:
        chrom = curr_area["region"]
        if curr_area["region"] != prev_region:
            cover_pos_by_chr[chrom] = dict()
            prev_end = -1
            prev_region = curr_area["region"]
        elif curr_area["start"] <= prev_end:
            overlap_start = curr_area["start"]
            overlap_end = prev_end
            if curr_area["end"] <= prev_end: # Inclusion
                overlap_end = curr_area["end"]
            for curr_pos in range(overlap_start, overlap_end):
                if curr_pos not in cover_pos_by_chr[chrom]:
                    cover_pos_by_chr[chrom][curr_pos] = 1
        if curr_area["end"] > prev_end:
            prev_end = curr_area["end"]
    # Convert overlapping positions in area
    cover_area_by_chr = dict() # overlapping positions area between area
    for chrom in cover_pos_by_chr:
        areas = list()
        prev_pos = None
        for curr_pos in sorted(cover_pos_by_chr[chrom]):
            if prev_pos is None or curr_pos != prev_pos + 1:
                if prev_pos is not None:
                    areas[-1]["end"] = prev_pos
                areas.append({
                    "region": chrom,
                    "start": curr_pos,
                    "end": None,
                })
            prev_pos = curr_pos
        if prev_pos is not None:
            areas[-1]["end"] = prev_pos
        cover_area_by_chr[chrom] = sorted(areas, key=lambda x: (x["start"], x["end"]))
    return cover_area_by_chr

def isOverlapping( regions_by_chr, chrom, start, end=None ):
    end = start if end is None else end
    is_overlapping = False
    if chrom in regions_by_chr:
        for region in regions_by_chr[chrom]:
            if region["start"] <= start and region["end"] >= end:
                is_overlapping = True
    return is_overlapping


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merges variants from several samples. If one variant is missing from a sample his AD, AF and DP are retrieved from the alignment file of this sample. The VCFs must come from the same process with same references.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--AF-precision', type=float, default=5, help="The AF's decimal precision. [Default: %(default)s]" )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-p', '--input-panel', required=True, help='The path to the *********************** without primers (format: BED).' )
    group_input.add_argument( '-i', '--input-variants', nargs='+', required=True, help='The path to the variants files (format: VCF).' )
    group_input.add_argument( '-a', '--input-aln', nargs='+', required=True, help='The path to the alignments files (format: BAM). Each alignment file correspond to a VCF.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()

    # Get overlapping area
    overlap_by_chr = getOverlaps( getSelectedArea(args.input_panel) )

    # Get identified variants from VCF
    variants = dict()
    sample_name = None
    out_vcf_info = None
    out_vcf_format = None
    for vcf_idx, current_vcf in enumerate(args.input_variants):
        with VCFIO(current_vcf) as FH_vcf:
            out_vcf_info = FH_vcf.info
            out_vcf_format = FH_vcf.format
            for record in FH_vcf: # For each variant
                for curr_spl in FH_vcf.samples: # For each sample in VCF
                    sample_name = curr_spl
                    vcaller_AF = record.get_AF( curr_spl )
                    vcaller_DP = record.get_DP( curr_spl )
                    for alt_idx, curr_alt in enumerate(record.alt): # For each alternative allele in in variant
                        record_allele = getAlleleRecord( FH_vcf, record, alt_idx )
                        # Get allele frequency from the variant caller
                        vcaller_curr_AF = vcaller_AF[alt_idx]
                        if len(vcaller_AF) == len(record.alt) + 1: # The AF cointains reference AF
                            vcaller_curr_AF = vcaller_AF[alt_idx + 1]
                        vcaller_curr_AD = int(vcaller_curr_AF*vcaller_DP)
                        record_allele.samples[curr_spl]["AF"] = [round(vcaller_curr_AF, args.AF_precision)]
                        record_allele.samples[curr_spl]["AD"] = [vcaller_curr_AD]
                        record_allele.samples[curr_spl]["DP"] = vcaller_DP
                        # Store allele
                        allele_id = record_allele.chrom + ":" + str(record_allele.pos) + "=" + record_allele.ref + "/" + record_allele.alt[0]
                        if allele_id not in variants:
                            # Add tags in format
                            if "AF" not in record_allele.format: record_allele.format.append("AF")
                            if "AD" not in record_allele.format: record_allele.format.append("AD")
                            if "DP" not in record_allele.format: record_allele.format.append("DP")
                            record_allele.format.append("OAAF")
                            record_allele.format.append("OADP")
                            # Add allele informations
                            variants[allele_id] = record_allele
                            variants[allele_id].samples[curr_spl]["OAID"] = [vcf_idx] # Temporary variable used to store amplicons groups where the variant has been see
                            variants[allele_id].samples[curr_spl]["OAAF"] = [round(vcaller_curr_AF, args.AF_precision)]
                            variants[allele_id].samples[curr_spl]["OADP"] = [vcaller_DP]
                        else:
                            variants[allele_id].samples[curr_spl]["AD"][0] += vcaller_curr_AD
                            variants[allele_id].samples[curr_spl]["DP"] += vcaller_DP
                            AF = variants[allele_id].samples[curr_spl]["AD"][0]/float(variants[allele_id].samples[curr_spl]["DP"])
                            variants[allele_id].samples[curr_spl]["AF"] = [round(AF, args.AF_precision)]
                            variants[allele_id].samples[curr_spl]["OAID"].append( vcf_idx ) # Temporary variable used to store amplicons groups where the variant has been see
                            OAAF = vcaller_curr_AD/float(vcaller_DP)
                            variants[allele_id].samples[curr_spl]["OAAF"].append( round(OAAF, args.AF_precision) )
                            variants[allele_id].samples[curr_spl]["OADP"].append( vcaller_DP )

    # Completes and writes variants
    with VCFIO(args.output_variants, "w") as FH_out:
        # Header
        FH_out.info = FH_vcf.info
        if "AF" in FH_out.info:
            del(FH_out.info["AF"])
        if "AD" in FH_out.info:
            del(FH_out.info["AD"])
        if "DP" in FH_out.info:
            del(FH_out.info["DP"])
        FH_out.format = FH_vcf.format
        FH_out.format["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies."}
        FH_out.format["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths."}
        FH_out.format["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Depth."}
        FH_out.format["OAAF"] = {"type": float, "type_tag": "Float", "number": None, "description": "Overlap amplicons alleles frequency."}
        FH_out.format["OADP"] = {"type": int, "type_tag": "Integer", "number": None, "description": "Overlap amplicons alleles depths."}
        FH_out.samples = [sample_name]
        FH_out._writeHeader()

        # Records
        for allele_id in variants:
            curr_var = variants[allele_id]
            ref_start = curr_var.pos
            ref_end = curr_var.pos + len(curr_var.ref) - 1
            if isOverlapping( overlap_by_chr, curr_var.chrom, ref_start, ref_end ): # Check only overlapping zone of interest otherwise the primers mask variants
                for idx_gp, aln_file in enumerate(args.input_aln):
                    if idx_gp not in curr_var.samples[curr_spl]["OAID"]:
                        # Retrieve AD, AF and DP from aln file
                        chrom_pos, ref_alt = allele_id.split("=")
                        chrom, pos = chrom_pos.split(":")
                        ref, alt = ref_alt.split("/")
                        OAAD, OADP = getADP( chrom, int(pos), ref, alt, aln_file )
                        curr_var.samples[curr_spl]["AD"][0] += OAAD
                        curr_var.samples[curr_spl]["DP"] += OADP
                        AF = curr_var.samples[curr_spl]["AD"][0]/float(curr_var.samples[curr_spl]["DP"])
                        curr_var.samples[curr_spl]["AF"] = [round(AF, args.AF_precision)]
                        OAAF = 0 if OADP == 0 else OAAD/float(OADP)
                        curr_var.samples[curr_spl]["OAAF"].append( round(OAAF, args.AF_precision) )
                        curr_var.samples[curr_spl]["OADP"].append( OADP )
            del(curr_var.samples[curr_spl]["OAID"])
            # Write variant
            if "AF" in curr_var.info:
                del(curr_var.info["AF"])
            if "AD" in curr_var.info:
                del(curr_var.info["AD"])
            if "DP" in curr_var.info:
                del(curr_var.info["DP"])
            FH_out.write( curr_var )

            #~ """
            #~ Simple case: overlap is identified, the variant does not overlap one primer and there are no repeat
                #~ ................*.PPPP
                           #~ PPPP.*................

                #~ ...............*PPPP
                           #~ PPPP*................

                #~ ..............I.PPPP
                          #~ PPPPI.................

                #~ ..............-.PPPP
                          #~ PPPP-.................

            #~ Complex case: limited by primer fixation and elongation
                #~ ..............**PPPP
                           #~ PPPP*................

            #~ Unmanaged case:
                #~ The overlap is not identified because it does not exist in the panel
                #~ ..............IPPPP
                          #~ PPPPI................
            #~ """
