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
def getADP( chrom, pos, ref, alt, aln_file ):
    """
    @summary: Returns the allele depth (AD) and the depth (DP) for the specified variant.
    @param chrom: [str] The variant region name.
    @param pos: [int] The variant position.
    @param ref: [str] The reference allele at the position.
    @param alt: [str] The variant allele at the position.
    @param aln_file: [str] The path to the alignment file (format: BAM). The AD and DP are retrieved from reads of this file. This file must be indexed.
    @return: [list] The first element is the AD, the second is the DP.
    @warning: Reads ID must be unique in SAM.
    """
    ref_start = pos
    ref_end = pos + len(ref.replace(".", "")) - 1

    FH_sam = pysam.AlignmentFile( aln_file, "rb" )
    inspect_start = ref_start - 1
    inspect_end = ref_end
    reads = dict()
        
    for pileupcolumn in FH_sam.pileup( chrom, inspect_start, inspect_end, max_depth=100000 ):
        for pileupread in pileupcolumn.pileups:
            if pileupcolumn.pos >= inspect_start and pileupcolumn.pos < inspect_end:
                # Get unique id
                read_id = pileupread.alignment.query_name
                if pileupread.is_read2:
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

def getAreaFromBed( bed_path ):
    area_by_chr = dict()
    # Group by chromosome
    with open(bed_path) as FH_in:
        for line in FH_in:
            row_data = [elt.strip() for elt in line.split("\t")]
            chrom = row_data[0]
            start = min( int(row_data[1]) + 1, int(row_data[2]) )
            end = max( int(row_data[1]) + 1, int(row_data[2]) )
            name = row_data[3]
            if chrom in kept_by_chr:
                kept_by_chr[chrom] = list()
            kept_by_chr[chrom].append( {"start": start, "end": end, "name":name} )
    # In each chromosome sort by position
    for chrom in area_by_chr:
        area_by_chr[chrom] = sorted(area_by_chr[chrom], key=lambda x: (x["start"], x["end"]))
    return area_by_chr

def getOverlappingPos( area_by_chr ): ##################################################### test
    cover_area_by_chr = dict() # overlapping position between area
    for chrom in area_by_chr:
        cover_area_by_chr[chrom] = dict()
        prev_end = -1
        prev_name = None
        for curr_area in area_by_chr[chrom]:
            if curr_area["start"] <= prev_end:
                overlap_start = curr_area["start"]
                overlap_end = prev_end
                if curr_area["end"] <= prev_end: # Inclusion
                    overlap_end = curr_area["end"]
                for curr_pos in range(overlap_start, overlap_end):
                    if curr_pos not in cover_area_by_chr[chrom]:
                        cover_area_by_chr[chrom][curr_pos] = [prev_name]
                    cover_area_by_chr[chrom][curr_pos].append( curr_area["name"] )
            if curr_area["end"] > prev_end:
                prev_end = curr_area["end"]
                prev_name = curr_area["name"]
    return cover_area_by_chr

def isOverlapping( regions_by_chr, chrom, pos ):
    is_overlapping = False
    if chrom in regions_by_chr:
        for region in regions_by_chr[chrom]:
            if region["start"] <= pos and region["end"] >= pos:
                is_overlapping = True
    return is_overlapping


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='****************************************************************.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    #~ group_input.add_argument( '-i', '--input-variants', nargs='+', required=True, help='********************************************************** (format: VCF).' )
    group_input.add_argument( '-p', '--panel', required=True, help='********************************************* (format: BED).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    #~ group_output.add_argument( '-o', '--output-variants', default="filtered.vcf", help='************************************************************* (format: VCF). [Default: %(default)s]' )
    args = parser.parse_args()

    # Process
    #~ area_by_chr = getAreaFromBed( args.panel )
    #~ overlapping_pos_by_chr = getOverlappingPos( area_by_chr )
    #~ variants_in_overlap = dict()
    #~ for curr_vcf in args.input_variants:
        #~ with VCFIO(curr_vcf) as FH_vcf:
            #~ for record in FH_vcf:
                #~ if record.chrom in overlapping_pos_by_chr:
                    #~ is_on_overlapping_pos = False
                    #~ for idx_alt in range(len(record.alt)):
                        #~ curr_allele = getAlleleRecord( FH_vcf, record, idx_alt )
                        #~ modified_pos = list()
                        #~ if len(curr_allele.ref) == 1:
                            #~ modified_pos = [curr_allele.pos]
                        #~ else:
                            #~ for curr_add_pos in range(len(curr_allele.ref)): # pb si une pos modif dans var (attention identique) est dans overlap #############################################
                                #~ modified_pos.append( curr_allele.pos + curr_add_pos )
                        #~ for curr_pos in modified_pos:
                            #~ if curr_pos in overlapping_pos_by_chr[chrom]:
                                #~ is_on_overlapping_pos = True
                    #~ if is_on_overlapping_pos: # Merge only strictly identical variants
                        #~ if chrom not in variants_in_overlap:
                            #~ variants_in_overlap[chrom] = {
                                #~ record.pos: [record]
                            #~ }
                        #~ elif record.pos not in variants_in_overlap[chrom]:
                            #~ variants_in_overlap[chrom][record.pos].append( record )
    area_by_chr = getAreaFromBed( args.panel )
    overlapping_pos_by_chr = getOverlappingPos( area_by_chr )
    variants_by_chr = dict()
    for curr_vcf in args.input_variants:
        with VCFIO(curr_vcf) as FH_vcf:
            for record in FH_vcf:
                if chrom not in variants_by_chr:
                    variants_by_chr[chrom] = {
                        record.pos: [record]
                    }
                elif record.pos not in variants_by_chr[chrom]:
                    variants_by_chr[chrom][record.pos].append( record )
    
    for curr_chrom in sorted(variants_by_chr):
        for curr_pos in sorted(variants_by_chr[chrom]):
			# Merge identical variants
			

                    """
                    Simple case: overlap is identified, the variant does not overlap one primer and there are no repeat
                        ................*.PPPP
                                   PPPP.*................
                    
                        ...............*PPPP
                                   PPPP*................
                               
                        ..............I.PPPP
                                  PPPPI.................

                        ..............-.PPPP
                                  PPPP-.................
                              
                    Complex case: limited by primer fixation and elongation
                        ..............**PPPP
                                   PPPP*................                    

                    Unmanaged case:
                        The overlap is not identified because it does not exist in the panel
                        ..............IPPPP
                                  PPPPI................
                    """
