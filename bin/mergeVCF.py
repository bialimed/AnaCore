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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFIO, getAlleleRecord



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

#~ def getADP( chrom, pos, ref, alt, aln_file ):
    #~ """
    #~ @summary: Returns the allele depth (AD) and the depth (DP) for the specified variant.
    #~ @param chrom: [str] The variant region name.
    #~ @param pos: [int] The variant position.
    #~ @param ref: [str] The reference allele at the position.
    #~ @param alt: [str] The variant allele at the position.
    #~ @param aln_file: [str] The path to the alignment file. The AD and DP are retrieved from reads of this file.
    #~ @return: [list] The first element is the AD, the second is the DP.
    #~ @warning: The sequence is returned as it is stored in the BAM file. Some mappers might have stored a reverse complement of the original read sequence.
    #~ """
    #~ ref_start = pos
    #~ ref_end = pos + len(ref.replace(".", "")) - 1

    #~ samfile = pysam.AlignmentFile( aln_file, "rb" )
    #~ inspect_start = ref_start - 1
    #~ inspect_end = ref_end
    #~ reads = list()
    #~ for curr_read in samfile.fetch( chrom, inspect_start, inspect_end ): # For each reads overlapping inspected area
        #~ inspected = [None for pos in range(inspect_start, inspect_end)]
        #~ prev_pos = None
        #~ for read_pos, ref_pos in curr_read.get_aligned_pairs(): # For position in the read alignment
            #~ if ref_pos is not None: # The read pos exists in reference
                #~ if ref_pos >= inspect_start and ref_pos < inspect_end: # If the pos is in inspected area
                    #~ inspected_pos = ref_pos - inspect_start
                    #~ if read_pos is None: # Deletion
                        #~ inspected[inspected_pos] = ""
                    #~ else: # Identical or substitution
                        #~ inspected[inspected_pos] = curr_read.query_sequence[read_pos].upper() # The sequence is returned as it is stored in the BAM file. Some mappers might have stored a reverse complement of the original read sequence.
            #~ else: # The read pos does not exist in reference
                #~ if prev_pos is not None: # Read alignment start by clipping
                    #~ if prev_pos >= inspect_start and prev_pos < inspect_end: # Insertion
                        #~ inspected_pos = prev_pos - inspect_start
                        #~ if inspected[inspected_pos] is None:
                            #~ inspected[inspected_pos] = ""
                        #~ inspected[inspected_pos] += curr_read.query_sequence[read_pos].upper() # The sequence is returned as it is stored in the BAM file. Some mappers might have stored a reverse complement of the original read sequence.
            #~ prev_pos = ref_pos if ref_pos is not None else prev_pos
        #~ reads.append( inspected )

    #~ # Process AD and DP
    #~ alt = alt if alt != "." else ""
    #~ AD = 0
    #~ DP = 0
    #~ for curr_read in reads:
        #~ if None not in curr_read: # Skip partial reads
            #~ DP += 1
            #~ if "".join(curr_read) == alt:
                #~ AD += 1

    #~ # Return
    #~ return( AD, DP )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merges variants from several samples. If one variant is missing from a sample his AD, AF and DP are retrieved from the alignment file of this sample. The VCFs must come from the same process with same references.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-p', '--AF-precision', type=float, default=5, help="The AF's decimal precision. [Default: %(default)s]" )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', nargs='+', required=True, help='The path to the variants files (format: VCF).' )
    group_input.add_argument( '-a', '--input-aln', nargs='+', required=True, help='The path to the alignments files (format: BAM). Each alignment file correspond to a VCF.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()

    # Get identified variants from VCF
    variants = dict()
    aln_by_samples = dict()
    out_vcf_info = None
    out_vcf_format = None
    for vcf_idx, current_vcf in enumerate(args.input_variants):
        current_aln = args.input_aln[vcf_idx]
        with VCFIO(current_vcf) as FH_vcf:
            out_vcf_info = FH_vcf.info
            out_vcf_format = FH_vcf.format
            for record in FH_vcf: # For each variant
                for curr_spl in FH_vcf.samples: # For each sample in VCF
                    aln_by_samples[curr_spl] = current_aln
                    vcaller_AF = record.get_AF( curr_spl )
                    vcaller_DP = record.get_DP( curr_spl )
                    for alt_idx, curr_alt in enumerate(record.alt): # For each alternative allele in in variant
                        record_allele = getAlleleRecord( FH_vcf, record, alt_idx )
                        # Get allele frequency from the variant caller
                        vcaller_curr_AF = vcaller_AF[alt_idx]
                        if len(vcaller_AF) == len(record.alt) + 1: # The AF cointains reference AF
                            vcaller_curr_AF = vcaller_AF[alt_idx + 1]
                        #~ #######################################################################################
                        #~ test_AD, test_DP = getADP( record_allele.chrom, record_allele.pos, record_allele.ref, record_allele.alt[0], current_aln )
                        #~ test_AF = float(test_AD)/test_DP
                        #~ print( record_allele.chrom + ":" + str(record_allele.pos), record_allele.ref + "/" + record_allele.alt[0], round(test_AF, 4), "vs", vcaller_curr_AF )
                        #~ print()
                        #~ #######################################################################################
                        record_allele.samples[curr_spl]["AF"] = [round(vcaller_curr_AF, args.AF_precision)]
                        record_allele.samples[curr_spl]["AD"] = [int(vcaller_curr_AF*vcaller_DP)]
                        record_allele.samples[curr_spl]["DP"] = vcaller_DP
                        # Store allele
                        allele_id = record_allele.chrom + ":" + str(record_allele.pos) + "=" + record_allele.ref + "/" + record_allele.alt[0]
                        if allele_id not in variants:
                            variants[allele_id] = record_allele
                        else:
                            variants[allele_id].samples[curr_spl] = record_allele.samples[curr_spl]

    # Completes and writes variants
    with VCFIO(args.output_variants, "w") as FH_out:
        # Header
        FH_out.info = FH_vcf.info
        FH_out.info["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies for the group of samples."}
        FH_out.info["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths for the group of samples."}
        FH_out.info["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Combined depth across samples."}
        FH_out.format = FH_vcf.format
        FH_out.format["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies."}
        FH_out.format["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths."}
        FH_out.format["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Depth."}
        FH_out.samples = [spl for spl in sorted(aln_by_samples)]
        FH_out._writeHeader()
        
        # Records
        for allele_id in variants:
            curr_var = variants[allele_id]
            # Add tag AF, AD and DP by sample
            if "AF" not in curr_var.format: curr_var.format.append("AF")
            if "AD" not in curr_var.format: curr_var.format.append("AD")
            if "DP" not in curr_var.format: curr_var.format.append("DP")
            # Process population AF, Ad and DP
            curr_var.info["AF"] = [0]
            curr_var.info["AD"] = [0]
            curr_var.info["DP"] = 0
            for spl in aln_by_samples:
                if spl not in curr_var.samples: # If the variant has not be seen in sample
                    # Retrieve AD, AF and DP from aln file
                    chrom_pos, ref_alt = allele_id.split("=")
                    chrom, pos = chrom_pos.split(":")
                    ref, alt = ref_alt.split("/")
                    AD, DP = getADP( chrom, int(pos), ref, alt, aln_by_samples[spl] )
                    curr_var.samples[spl] = {
                        "AF": [round(float(AD)/DP, args.AF_precision)],
                        "AD": [AD],
                        "DP": DP
                    }
                curr_var.info["AD"][0] += curr_var.samples[spl]["AD"][0]
                curr_var.info["DP"] += curr_var.samples[spl]["DP"]
            curr_var.info["AF"][0] = curr_var.info["AD"][0] / curr_var.info["DP"]
            # Write variant
            FH_out.write( curr_var )
