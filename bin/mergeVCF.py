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
__version__ = '1.1.2'
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
def getADPReads( chrom, pos, ref, alt, aln_file, selected_RG=None ):
    """
    @summary: Returns the allele depth (AD) and the depth (DP) for the specified variant. These counts are expressed in number of reads: if the R1 and the R2 of a sequence has overlaps the variant, each is counted.
    @param chrom: [str] The variant region name.
    @param pos: [int] The variant position.
    @param ref: [str] The reference allele at the position.
    @param alt: [str] The variant allele at the position.
    @param aln_file: [str] The path to the alignment file (format: BAM). The AD and DP are retrieved from reads of this file. This file must be indexed.
    @param selected_RG: [list] The ID of RG used in AD and DP. Default: all read groups.
    @returns: [list] The first element is the AD, the second is the DP.
    @warning: Reads ID must be unique in SAM. These counts are expressed in number of reads: if the R1 and the R2 of a sequence has overlaps the variant, each is counted.
    """
    ref_start = pos
    ref_end = pos + len(ref.replace(".", "")) - 1
    if selected_RG is not None: selected_RG = {RG:1 for RG in selected_RG}
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
    # Completes downstream positions
    inspected_len = inspect_end - inspect_start
    for read_id in reads:
        read_len = len(reads[read_id])
        for idx in range(inspected_len - read_len):
            reads[read_id].append( None )
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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merges variants from several samples. If one variant is missing from a sample his AD, AF and DP are retrieved from the alignment file of this sample. The VCFs must come from the same process with same references. ote: for a common variant all the fields values except for AF, AD and DP are retrieved from the first VCF where it has been found.' )
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
    for vcf_idx, current_vcf in enumerate(args.input_variants):
        current_aln = args.input_aln[vcf_idx]
        with VCFIO(current_vcf) as FH_vcf:
            # Manage samples
            for curr_spl in FH_vcf.samples: # For each sample in VCF
                aln_by_samples[curr_spl] = current_aln
            # Manage records
            for record in FH_vcf: # For each variant
                for curr_spl in FH_vcf.samples: # For each sample in VCF
                    vcaller_AF = record.getAF( curr_spl )
                    vcaller_DP = record.getDP( curr_spl )
                    for alt_idx, curr_alt in enumerate(record.alt): # For each alternative allele in in variant
                        record_allele = getAlleleRecord( FH_vcf, record, alt_idx )
                        # Get allele frequency from the variant caller
                        vcaller_curr_AF = vcaller_AF[alt_idx]
                        if len(vcaller_AF) == len(record.alt) + 1: # The AF cointains reference AF
                            vcaller_curr_AF = vcaller_AF[alt_idx + 1]
                        #######################################################################################
                        #~ test_AD, test_DP = getADP( record_allele.chrom, record_allele.pos, record_allele.ref, record_allele.alt[0], current_aln )
                        #~ test_AF = 0 if test_DP == 0 else float(test_AD)/test_DP
                        #~ print( record_allele.chrom + ":" + str(record_allele.pos), record_allele.ref + "/" + record_allele.alt[0], round(test_AF, 4), "vs", vcaller_curr_AF )
                        #~ print()
                        #######################################################################################
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
        FH_out.copyHeader( FH_in )
        FH_out.info["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies for the group of samples."}
        FH_out.info["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths for the group of samples."}
        FH_out.info["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Combined depth across samples."}
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
                    AD, DP = getADPReads( chrom, int(pos), ref, alt, aln_by_samples[spl] )
                    curr_var.samples[spl] = {
                        "AF": [0 if DP == 0 else round(float(AD)/DP, args.AF_precision)],
                        "AD": [AD],
                        "DP": DP
                    }
                curr_var.info["AD"][0] += curr_var.samples[spl]["AD"][0]
                curr_var.info["DP"] += curr_var.samples[spl]["DP"]
            curr_var.info["AF"][0] = curr_var.info["AD"][0] / curr_var.info["DP"]
            # Write variant
            FH_out.write( curr_var )