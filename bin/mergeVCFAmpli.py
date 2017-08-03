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

def getSelectedAreas( in_design ):
    """
    @summary: Returns the list of selected areas from a BED file.
    @param in_design: [str] The path to the selected areas description (format: BED).
    @return: [list] The list of BED's areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    """
    selected_areas = list()
    with open(in_design) as FH_design:
        for line in FH_design:
            if not line.startswith("browser ") and not line.startswith("track ") and not line.startswith("#"):
                fields = [elt.strip() for elt in line.split("\t")]
                selected_areas.append({
                    "region": fields[0],
                    "start": int(fields[1]) +1, # Start in BED is 0-based
                    "end": int(fields[2]),
                    "id": fields[3]
                })
    return( selected_areas )

def getAreasContaining( areas, chrom, start, end ):
    """
    @summary: Returns the ID of the areas containing the specidfied segment (contiguous positions).
    @param areas: [list] List of evaluated area. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    @param chrom: [str] The chromosome of the evaluated segment.
    @param start: The start position of the evaluated segment.
    @param end: The end position of the evaluated segment.
    @returns: [list] IDs of the areas containing the segment.
    """
    containers = list()
    for curr_area in areas:
        if contains( curr_area, chrom, start, end ):
            containers.append( curr_area["id"] )
    return containers

def contains( area, chrom, start, end=None ):
    """
    @summary: Returns True if the segment (contiguous positions) is contained in area.
    @param area: [dict] The evaluated area. Format: {"region":"chr1", "start":501, "end":608}.
    @param chrom: [str] The chromosome of the evaluated segment.
    @param start: The start position of the evaluated segment.
    @param end: The end position of the evaluated segment.
    @returns: [bool] True if the segment is contained in area.
    """
    contains = False
    end = start if end is None else end
    if area["region"] == chrom and area["start"] <= start and area["end"] >= end:
        contains = True
    return contains

def getRGIdByRGTag( in_aln, tag, selected_value ):
    """
    @summary: Returns the IDs of RG with a tag value in selected values.
    @param in_aln: [str] The path to the alignment file (format: BAM).
    @param tag: [str] The RG tag used in filter.
    @param selected_value: [list] The authorized values for RG tag.
    @returns: [list] IDs of the corresponding reads groups.
    """
    RG_id = list()
    with pysam.AlignmentFile( in_aln, "rb" ) as FH_sam:
        for RG in FH_sam.header["RG"]:
            if RG[tag] in selected_value:
                RG_id.append( RG["ID"] )
    return RG_id


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merges variants from several samples. If one variant is missing from a sample his AD, AF and DP are retrieved from the alignment file of this sample. The VCFs must come from the same process with same references. Note: for a common variant all the fields values except for AF, AD and DP are retrieved from the first VCF where it has been found.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--AF-precision', type=float, default=5, help="The AF's decimal precision. [Default: %(default)s]" )
    parser.add_argument( '-t', '--RG-tag', default='LB', help='RG tag used to store the area ID. [Default: %(default)s]' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-p', '--input-designs', nargs='+', required=True, help='The path to the amplicons design. The start and end of the amplicons must be without primers (format: BED).' )
    group_input.add_argument( '-i', '--input-variants', nargs='+', required=True, help='The path to the variants files (format: VCF).' )
    group_input.add_argument( '-a', '--input-aln', nargs='+', required=True, help='The path to the alignments files (format: BAM). Each alignment file correspond to a VCF.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()

    # Get identified variants from VCF
    variants = dict()
    aln_by_samples = dict()
    design_by_samples = dict()
    for vcf_idx, current_vcf in enumerate(args.input_variants):
        current_aln = args.input_aln[vcf_idx]
        current_design = getSelectedAreas( args.input_designs[vcf_idx] )
        with VCFIO(current_vcf) as FH_vcf:
            for record in FH_vcf: # For each variant
                for curr_spl in FH_vcf.samples: # For each sample in VCF
                    aln_by_samples[curr_spl] = current_aln
                    design_by_samples[curr_spl] = current_design
                    vcaller_AF = record.get_AF( curr_spl )
                    vcaller_DP = record.get_DP( curr_spl )
                    for alt_idx, curr_alt in enumerate(record.alt): # For each alternative allele in variant
                        record_allele = getAlleleRecord( FH_vcf, record, alt_idx )
                        # Get allele frequency from the variant caller
                        vcaller_curr_AF = vcaller_AF[alt_idx]
                        if len(vcaller_AF) == len(record.alt) + 1: # The AF contains reference AF
                            vcaller_curr_AF = vcaller_AF[alt_idx + 1]
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
            # Process population AF, AD and DP
            curr_var.info["AF"] = [0]
            curr_var.info["AD"] = [0]
            curr_var.info["DP"] = 0
            for spl in aln_by_samples:
                if spl not in curr_var.samples: # If the variant has not be seen in sample
                    # Get valid RG
                    ref_end = curr_var.pos + len(curr_var.ref) - 1
                    overlapped_ampl = getAreasContaining( design_by_samples[spl], curr_var.chrom, curr_var.pos, ref_end )
                    overlapped_RG = getRGIdByRGTag( aln_by_samples[spl], args.RG_tag, overlapped_ampl )
                    # Retrieve AD, AF and DP from aln file
                    AD, DP = getADP( curr_var.chrom, curr_var.pos, curr_var.ref, curr_var.alt[0], aln_by_samples[spl], overlapped_RG )
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
