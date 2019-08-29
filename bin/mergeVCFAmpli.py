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
__version__ = '2.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.bed import getAreasByChr
from anacore.region import Region
from anacore.vcf import VCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
class ConsensusException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def getADPContig(chrom, pos, ref, alt, aln_file, selected_RG=None):
    """
    @summary: Returns the allele depth (AD) and the depth (DP) for the specified variant. These counts are expressed in number of contig: if the R1 and the R2 of a sequence has overlaps the variant they only the consensus is counted for the pair.
    @param chrom: [str] The variant region name.
    @param pos: [int] The variant position.
    @param ref: [str] The reference allele at the position.
    @param alt: [str] The variant allele at the position.
    @param aln_file: [str] The path to the alignment file (format: BAM). The AD and DP are retrieved from reads of this file. This file must be indexed.
    @param selected_RG: [list] The ID of RG used in AD and DP. Default: all read groups.
    @returns: [list] The first element is the AD, the second is the DP. These counts are expressed in number of contig: if the R1 and the R2 of a sequence has overlaps the variant they only the consensus is counted for the pair.
    """
    if ref == "-":
        raise Exception('The method getADPContig() cannot be used on insertion with dash as reference (instead of "-/AGC" prefer the standard fromat like "T/TAGC").')
    # Retrieve reads
    ref_start = pos
    ref_end = pos + len(ref) - 1
    inspect_start = ref_start - 1
    inspect_end = ref_end
    reads, quals = getAlnAndQual(aln_file, chrom, inspect_start, inspect_end, selected_RG, 100000)
    # Process AD and DP
    alt = alt if alt != "-" else ""
    AD = 0
    DP = 0
    for read_id in reads:
        consensus = ""
        try:
            consensus = getSimplePairConsensus(reads[read_id], quals[read_id])
        except ConsensusException:
            # Calculate alignment cost and choose the lowest
            cost_aln_R1 = getAlnCost(ref, reads[read_id]["R1"])
            cost_aln_R2 = getAlnCost(ref, reads[read_id]["R2"])
            consensus = reads[read_id]["R1"] if cost_aln_R1 <= cost_aln_R2 else reads[read_id]["R2"]
        if None not in consensus:  # Skip partial reads
            DP += 1
            if "".join(consensus) == alt:
                AD += 1
    # Return
    return AD, DP

def getADPReads(chrom, pos, ref, alt, aln_file, selected_RG=None):
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
    ref_end = pos + len(ref.replace("-", "")) - 1
    if selected_RG is not None: selected_RG = {RG:1 for RG in selected_RG}
    # Retrieve reads
    inspect_start = ref_start - 1
    inspect_end = ref_end
    reads = dict()
    with pysam.AlignmentFile(aln_file, "rb") as FH_sam:
        for pileupcolumn in FH_sam.pileup(chrom, inspect_start, inspect_end, max_depth=100000):
            for pileupread in pileupcolumn.pileups:
                if selected_RG is None or (pileupread.alignment.get_tag("RG") in selected_RG):
                    if not pileupread.alignment.is_secondary:
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
                            if pileupread.is_del:  # Deletion
                                reads[read_id].append("")
                            elif pileupread.indel > 0:  # Insertion
                                insert = ""
                                for insert_idx in range(pileupread.indel + 1):
                                    insert += pileupread.alignment.query_sequence[pileupread.query_position + insert_idx].upper()
                                reads[read_id].append(insert)
                            elif not pileupread.is_refskip:  # Substitution
                                reads[read_id].append(
                                    pileupread.alignment.query_sequence[pileupread.query_position].upper()
                                )
    # Completes downstream positions
    inspected_len = inspect_end - inspect_start
    for read_id in reads:
        read_len = len(reads[read_id])
        for idx in range(inspected_len - read_len):
            reads[read_id].append(None)
    # Process AD and DP
    alt = alt if alt != "-" else ""
    AD = 0
    DP = 0
    for read_id in reads:
        if None not in reads[read_id]:  # Skip partial reads
            DP += 1
            if "".join(reads[read_id]) == alt:
                AD += 1
    # Return
    return AD, DP

def getAlnCost(ref, aln_seq, weights=None):
    """
    @summary: Returns cost of the alignment for alignment comparison. For two alignments on the
              same area the alignment with the lower cost is the nearest of the reference.
    @param ref: [str] The reference sequenece on the interest area. It will be used for
                differenciate matches and substitutions.
                Example:
                  Ref:   NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                  Read:  NNNNNNNNNAGAT--GGCTTACC
                  ref =  "AGATAAGGCCCA"
    @param aln_seq: [list] The read alignment on the interest area.
                    Example:
                      Ref:   NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                      Read:  NNNNNNNNNAGAT--GGCTTACC
                      aln_seq = ["A",  "G", "A", "T", "", "", "G", "G", "CTTA", "C", "C", None]
    @param weights: [dict] The weights of each evenments. Each evenment must be setted.
                    Default:
                      { "del": 3, "del_ext": 0.5, "ins": 3, "ins_ext": 0.5, "miss": 1, "subs": 1, "match": -0.5 }
    @return: [float] The cost of the alignment.
    """
    if weights is None:
        weights = {"del": 3, "del_ext": 0.5, "ins": 3, "ins_ext": 0.5, "miss": 1, "subs": 1, "match": -0.5}
    nb = {"del": 0, "del_ext": 0, "ins": 0, "ins_ext": 0, "miss": 0, "subs": 0, "match": 0}
    prev = ""
    for idx_nt, nt in enumerate(aln_seq):
        if nt is None:
            nb["miss"] += 1
            prev = "miss"
        elif nt == "":
            if prev == "del":
                nb["del_ext"] += 1
            else:
                nb["del"] += 1
            prev = "del"
        elif len(nt) > 1:
            nb["ins"] += 1
            nb["ins_ext"] = len(nt) - 2
            prev = "ins"
        elif ref[idx_nt] != nt:
            nb["subs"] += 1
            prev = "subs"
        else:
            nb["match"] += 1
            prev = "match"
    cost = sum([(nb[category] * weights[category]) for category in nb])
    return cost

def getAlnAndQual(aln_file, chrom, inspect_start, inspect_end, selected_RG, max_depth=100000):
    """
    @summary: Returns for each reads in inspected area the fragment of the alignment corresponding.
              Two elements are returned: by read ID the sequence alignment fragment for R1 and R2
              and by read ID the bases qualities in alignment fragment for R1 and R2.
              The sequence alignment is a list where each element corresponds to one position on
              reference. The deletions in reads are marked with en empty string (AT/A becomes "" at
              T position) ; the insertions are marked by a string with a length superior than 1 on
              the previous position (A/ATG becomes "ATG" at A position) ; Matches and substitutions
              are marked by a string with length equals to 1 (A/G becomes "G" at A position) ;
              Missing positions are marked by a None.
              The bases qualities in alignment is a list where each element corresponds to quality
              of the base(s) corresponding to one position in reference. The deletions and missing
              positions are marked by None ; the insertions are marked by a list of qualities.
    @param aln_file: [str] The path to the alignment file (format: BAM). This file must be indexed.
    @param chrom: [str] The reference region name.
    @param inspect_start: [int] The start position of the inspected area (0-based).
    @param inspect_end: [int] The end position of the inspected area (1-based).
    @param inspect_end: [int] The end position of the inspected area (1-based).
    @param selected_RG: [list] The retained RG. Default: all read groups.
    @param max_depth: [int] Maximum read depth processed.
    @returns: [list] First element is by read ID the sequence alignment fragment for R1 and R2.
              Second element is by read ID the base quality in alignment fragment for R1 and R2.
              Input:
                  Pos:                10           20
                                      |            |
                  chr1:      NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                  seq_A R1:  NNNNNNNNNAGAT--GACTTACC
                  seq_A R2:            GAT--GACTTACCANNNNNNNN
                  seq_B R1:  NNNNNNNNNAGATAAGGCTTAC
              Parameters:
                  chrom = "chr1"
                  inspect_start = 10
                  inspect_end = 21
                  selected_RG = None
              Returns:
                  reads = {
                    "seq_A":{
                        "R1": ["A",  "G", "A", "T", "", "", "G", "G", "CTTA", "C", "C", None],
                        "R2": [None, "G", "A", "T", "", "", "G", "A", "CTTA", "C", "C", "A"]
                    },
                    "seq_B":{
                        "R1": ["A", "G", "A", "T", "A", "A", "G", "G", "CTTA", "C", None, None]
                    }
                  }
                  quals = {
                    "seq_A":{
                        "R1": [35,   35, 35, 34, None, None, 34, 34, [33, 33, 33, 32], 32, 31, None],
                        "R2": [None, 31, 32, 32, None, None, 33, 33, [33, 34, 34, 34], 35, 35, 35]
                    },
                    "seq_B":{
                         "R1": [31, 31, 30, 30, 30, 29, 29, 29, [28, 28, 29, 28], 27, None, None]
                    }
                  }
    @warning: Insertions are stored as an alternative of the previous base (example: A/ATC becomes
              "ATC" at A position). If you want to extract insertion from alignment take care of
              provide inspected_start equal to the previous position of the insertion (example:
              position of A in A/ATC).
    """
    if selected_RG is not None: selected_RG = {RG:1 for RG in selected_RG}
    reads = dict()
    quals = dict()
    with pysam.AlignmentFile(aln_file, "rb") as FH_sam:
        for pileupcolumn in FH_sam.pileup(chrom, inspect_start, inspect_end, max_depth=max_depth):
            for pileupread in pileupcolumn.pileups:
                if selected_RG is None or (pileupread.alignment.get_tag("RG") in selected_RG):
                    if not pileupread.alignment.is_secondary:
                        if pileupcolumn.pos >= inspect_start and pileupcolumn.pos < inspect_end:
                            # Get id
                            read_id = pileupread.alignment.query_name
                            pair_id = "R1" if not pileupread.alignment.is_read2 else "R2"
                            # Store new reads
                            if read_id not in reads:
                                reads[read_id] = dict()
                                quals[read_id] = dict()
                            if pair_id not in reads[read_id]:
                                reads[read_id][pair_id] = [None for pos in range(inspect_start, pileupcolumn.pos)]
                                quals[read_id][pair_id] = [None for pos in range(inspect_start, pileupcolumn.pos)]
                            curr_read = reads[read_id][pair_id]
                            curr_qual = quals[read_id][pair_id]
                            # Store comparison with ref for current position
                            if pileupread.is_del:  # Deletion
                                curr_read.append("")
                                curr_qual.append("")
                            elif pileupread.indel > 0:  # Insertion
                                insert = ""
                                insert_qual = list()
                                for insert_idx in range(pileupread.indel + 1):
                                    insert += pileupread.alignment.query_sequence[pileupread.query_position + insert_idx].upper()
                                    insert_qual.append(pileupread.alignment.query_qualities[pileupread.query_position + insert_idx])
                                curr_read.append(insert)
                                curr_qual.append(insert_qual)
                            elif not pileupread.is_refskip:  # Substitution
                                curr_read.append(
                                    pileupread.alignment.query_sequence[pileupread.query_position].upper()
                                )
                                curr_qual.append(
                                    pileupread.alignment.query_qualities[pileupread.query_position]
                                )
    # Completes downstream positions
    inspected_len = inspect_end - inspect_start
    for read_id in reads:
        for pair_id in reads[read_id]:
            read_len = len(reads[read_id][pair_id])
            for idx in range(inspected_len - read_len):
                reads[read_id][pair_id].append(None)
                quals[read_id][pair_id].append(None)
    return reads, quals

def getRGIdByRGTag(in_aln, tag, selected_value):
    """
    @summary: Returns the IDs of RG with a tag value in selected values.
    @param in_aln: [str] The path to the alignment file (format: BAM).
    @param tag: [str] The RG tag used in filter.
    @param selected_value: [list] The authorized values for RG tag.
    @returns: [list] IDs of the corresponding reads groups.
    """
    RG_id = list()
    with pysam.AlignmentFile(in_aln, "rb") as FH_sam:
        for RG in FH_sam.header["RG"]:
            if RG[tag] in selected_value:
                RG_id.append(RG["ID"])
    return RG_id

def getSimplePairConsensus(seq, qual):
    """
    @summary: Returns pair consensus between alignment fragments from R1 and R2 of the same matrix. The consensus is based on the best base quality at position.
    @param seq: [dict] The sequence alignment fragment for R1 and R2 (see getAlnAndQual()).
                Example:
                  Ref: NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                  R1:  NNNNNNNNNAGAT--GGCTTACC
                                 |||  |*||||||
                  R2:            GAT--GACTTACCANNNNNNNN
                  seq = {
                     "R1": ["A",  "G", "A", "T", "", "", "G", "G", "CTTA", "C", "C", None],
                     "R2": [None, "G", "A", "T", "", "", "G", "A", "CTTA", "C", "C", "A"]
                  }
    @param qual: [dict] The quality of bases in alignment fragment for R1 and R2 (see getAlnAndQual()).
                 Example:
                   Ref: NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                   R1:  NNNNNNNNNAGAT--GGCTTACC
                                  |||  |*||||||
                   R2:            GAT--GACTTACCANNNNNNNN
                   Q1:  FFFFFFFFFFEDC  BA@:::?>
                   Q2:            <=>  ?0A:::BCDDDDDDDDD
                   qual = {
                      "R1": [37, 36, 35, 34, None, None, 33, 32, [31, 25, 25, 25], 30, 29, None],
                      "R2": [None, 27, 28, 29, None, None, 30, 15, [32, 25, 25, 25], 33, 34, 35]
                   }
    @returns: [list] The consensus sequence alignment fragment.
              Example:
                Ref: NNNNNNNNNAGATAAGGC---CCANNNNNNNN
                R1:  NNNNNNNNNAGAT--GGCTTACC
                               |||  |*||||||
                R2:            GAT--GACTTACCANNNNNNNN
                Q1:  FFFFFFFFFFEDC  BA@:::?>
                Q2:            <=>  ?0A:::BCDDDDDDDDD
                consensus = ["A", "G", "A", "T", "", "", "G", "G", "CTTA", "C", "C", "A"],
    """
    consensus = list()
    if "R1" in seq and "R2" in seq:  # R1 and R2 overlap inspected region
        inspect_len = len(seq["R1"])
        for idx in range(inspect_len):
            if seq["R1"][idx] is not None and seq["R2"][idx] is not None:  # R1 and R2 overlap position
                if seq["R1"][idx] == seq["R2"][idx]:  # The sequence at the current position is the same between two reads
                    consensus.append(seq["R1"][idx])
                else:  # R1 and R2 differ on position
                    if len(seq["R1"][idx]) == len(seq["R2"][idx]):  # The sequences are different but with the same frame
                        if len(seq["R1"][idx]) < 2:  # Substitution or deletion
                            if qual["R1"][idx] is None or qual["R1"][idx] >= qual["R2"][idx]:  # deletion or qual R1 superior
                                consensus.append(seq["R1"][idx])
                            else:
                                consensus.append(seq["R2"][idx])
                        else:  # Insertion
                            insert = ""
                            for idx_ins, qual_R1 in enumerate(qual["R1"][idx]):
                                qual_R2 = qual["R1"][idx][idx_ins]
                                if qual_R1 >= qual_R2:
                                    insert += seq["R1"][idx][idx_ins]
                                else:
                                    insert += seq["R2"][idx][idx_ins]
                            consensus.append(insert)
                    else:
                        raise ConsensusException("Simple consensus cannot be found between {} and {}.".format(seq["R1"], seq["R2"]))
            elif seq["R1"][idx] is not None:  # Only the R1 overlaps position
                consensus.append(seq["R1"][idx])
            elif seq["R2"][idx] is not None:  # Only the R2 overlaps position
                consensus.append(seq["R2"][idx])
            else:
                consensus.append(None)
    elif "R1" in seq:  # Only R1 overlaps inspected region
        consensus = seq["R1"]
    else:  # Only R2 overlaps inspected region
        consensus = seq["R2"]
    return consensus



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merges variants from several samples. If one variant is missing from a sample his AD, AF and DP are retrieved from the alignment file of this sample. The VCFs must come from the same process with same references. Note: for a common variant all the fields values except for AF, AD and DP are retrieved from the first VCF where it has been found.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-f', '--AF-precision', type=float, default=5, help="The AF's decimal precision. [Default: %(default)s]")
    parser.add_argument('-t', '--RG-tag', default='LB', help='RG tag used to store the area ID. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-p', '--input-designs', nargs='+', required=True, help='The path to the amplicons design (format: BED). The start and end of the amplicons must be without primers.')
    group_input.add_argument('-i', '--input-variants', nargs='+', required=True, help='The path to the variants files (format: VCF).')
    group_input.add_argument('-a', '--input-aln', nargs='+', required=True, help='The path to the alignments files (format: BAM). Each alignment file correspond to a variants file.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))
    log.info("Version: " + str(__version__))

    # Get identified variants from VCF
    variants = dict()
    aln_by_samples = dict()
    design_by_samples = dict()
    for vcf_idx, current_vcf in enumerate(args.input_variants):
        current_aln = args.input_aln[vcf_idx]
        current_design = getAreasByChr(args.input_designs[vcf_idx])
        with VCFIO(current_vcf) as FH_vcf:
            # Manage samples
            for curr_spl in FH_vcf.samples:  # For each sample in VCF
                aln_by_samples[curr_spl] = current_aln
                design_by_samples[curr_spl] = current_design
            # Manage records
            for record in FH_vcf:  # For each variant
                for curr_spl in FH_vcf.samples:  # For each sample in VCF
                    vcaller_AF = record.getAltAF(curr_spl)
                    vcaller_DP = record.getDP(curr_spl)
                    for alt_idx, curr_alt in enumerate(record.alt):  # For each alternative allele in variant
                        record_allele = getAlleleRecord(FH_vcf, record, alt_idx)
                        # Get allele frequency from the variant caller
                        vcaller_curr_AF = vcaller_AF[alt_idx]
                        if len(vcaller_AF) == len(record.alt) + 1:  # The AF contains reference AF
                            vcaller_curr_AF = vcaller_AF[alt_idx + 1]
                        record_allele.samples[curr_spl]["AF"] = [round(vcaller_curr_AF, args.AF_precision)]
                        record_allele.samples[curr_spl]["AD"] = [int(vcaller_curr_AF * vcaller_DP)]
                        record_allele.samples[curr_spl]["DP"] = vcaller_DP
                        # Store allele
                        allele_id = record_allele.getName()
                        if allele_id not in variants:
                            variants[allele_id] = record_allele
                        else:
                            variants[allele_id].samples[curr_spl] = record_allele.samples[curr_spl]

    # Completes and writes variants
    with VCFIO(args.output_variants, "w") as FH_out:
        # Header
        FH_out.copyHeader(FH_vcf)
        FH_out.info["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies for the group of samples."}
        FH_out.info["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths for the group of samples."}
        FH_out.info["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Combined depth across samples."}
        FH_out.format["AF"] = {"type": float, "type_tag": "Float", "number": None, "number_tag": "A", "description": "The alleles frequencies."}
        FH_out.format["AD"] = {"type": int, "type_tag": "Integer", "number": None, "number_tag": "A", "description": "The alleles depths."}
        FH_out.format["DP"] = {"type": int, "type_tag": "Integer", "number": 1, "description": "Depth."}
        FH_out.samples = [spl for spl in sorted(aln_by_samples)]
        FH_out.writeHeader()

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
                if spl not in curr_var.samples:  # If the variant has not be seen in sample
                    AD = 0
                    DP = 0
                    # Get valid RG
                    ref_end = curr_var.pos + len(curr_var.ref) - 1
                    curr_var_region = Region(curr_var.pos, ref_end, None, curr_var.chrom)
                    overlapped_ampl = list()
                    if curr_var.chrom in design_by_samples[spl]:
                        overlapped_ampl = design_by_samples[spl][curr_var.chrom].getContainers(curr_var_region)
                    if len(overlapped_ampl) > 0:
                        # Retrieve AD, AF and DP from aln file
                        overlapped_ampl_name = [ampl.name for ampl in overlapped_ampl]
                        overlapped_RG = getRGIdByRGTag(aln_by_samples[spl], args.RG_tag, overlapped_ampl_name)
                        AD, DP = getADPContig(curr_var.chrom, curr_var.pos, curr_var.ref, curr_var.alt[0], aln_by_samples[spl], overlapped_RG)
                    # Store AD, AF and DP for sample
                    curr_var.samples[spl] = {
                        "AF": [0 if DP == 0 else round(AD / DP, args.AF_precision)],
                        "AD": [AD],
                        "DP": DP
                    }
                curr_var.info["AD"][0] += curr_var.samples[spl]["AD"][0]
                curr_var.info["DP"] += curr_var.samples[spl]["DP"]
            curr_var.info["AF"][0] = round(curr_var.info["AD"][0] / curr_var.info["DP"], args.AF_precision)
            # Write variant
            FH_out.write(curr_var)
    logging.info("End of process")
