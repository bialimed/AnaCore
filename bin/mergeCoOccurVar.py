#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse
from statistics import mean
from anacore.sequenceIO import FastaIO
from anacore.vcf import VCFIO, VCFRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAlnCmp(read, ref_seq):
    """
    Return the dense representation of the alignment between reference sequence and read.

    :param read: The alignment.
    :type read: pysam.AlignedSegment
    :param ref_seq: The reference sequence included in the alignement.
    :type ref_seq: str
    :return: Dense representation of the alignment. First is reference alignment and second is read alignment.
    :rtype: (list, list)
    """
    ref_aln = []
    read_aln = []
    ref_seq = [elt for elt in ref_seq]  # reference sequence on alignment (read.get_reference_sequence() can be incorrect)
    read_seq = [elt for elt in read.query_alignment_sequence]  # query sequence without clipped
    for operation_id, operation_lg in read.cigartuples:
        for operation_pos in range(operation_lg):
            if operation_id != 5 and operation_id != 4:  # Is not clipped
                if operation_id in [0, 7, 8]:  # Match or mismatch
                    read_aln.append(read_seq.pop(0))
                    ref_aln.append(ref_seq.pop(0))
                elif operation_id == 1:  # Insertion
                    read_aln[-1] += read_seq.pop(0)
                elif operation_id == 2:  # Deletion
                    read_aln.append("")
                    ref_aln.append(ref_seq.pop(0))
                elif operation_id == 3:  # Refskip ########################### Are the introns in the variant ?
                    # next nt in ref but current read
                    ref_aln.pop(0)
                elif operation_id == 9:  # Back (www.seqanswers.com/forums/showthread.php?t=34440)
                    raise Exception("Parsing error on read {}. The management for the CIGAR operator B is not implemented.".format(read.query_name))
                # elif operation_id == 6:  # Padding
                #     pass
    return ref_aln, read_aln


def setRefPos(variant):
    """
    Add start and end attributes in VCFRecord. For insertions the start is defined on the first position before the insertion and the end on the last position affected by the insertion.

    :param variant: The variant to update.
    :type variant: anacore.vcf.VCFRecord
    """
    variant.start = int(variant.refStart())
    if variant.isInsertion():
        if variant.start == int(variant.refStart() + 0.5):
            variant.start -= 1
    variant.end = int(variant.refEnd())


def setSupportingReads(prev, curr, chrom_seq, FH_aln, log):
    """
    Add supporting reads for each variant. Only the reads overlapping the positions of the two variants are take into account.

    :param prev: The upstream variant to update.
    :type prev: anacore.vcf.VCFRecord updated with setRefPos() and isIns
    :param curr: The donstream variant to update.
    :type curr: anacore.vcf.VCFRecord updated with setRefPos() and isIns
    :param chrom_seq: The sequence of the chromosome.
    :type chrom_seq: str
    :param FH_aln: The file handle to the alignments file. The variants must have been defined from this alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param log: The logger object.
    :type log: logging.Logger
    """
    prev.supporting_reads = set()
    curr.supporting_reads = set()
    for read in FH_aln.fetch(curr.chrom, prev.start, curr.end):
        if not read.is_duplicate:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                overlap_the_two = (ref_start <= prev.start and ref_end >= prev.end) and (ref_start <= curr.start and ref_end >= curr.end)
                if overlap_the_two:
                    ref_aln, read_aln = getAlnCmp(read, chrom_seq[ref_start - 1:ref_end])
                    for var in [prev, curr]:
                        alt = read_aln[var.start - ref_start:var.end - ref_start + 1]
                        ref = ref_aln[var.start - ref_start:var.end - ref_start + 1]
                        if var.isIns:
                            while len(ref[0]) > 0 and len(alt[0]) > 0 and alt[0][0] == ref[0][0]:
                                alt[0] = alt[0][1:]
                                ref[0] = ref[0][1:]
                        var_alt = var.alt[0].upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                        var_ref = var.ref.upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                        if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # Is present
                            log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                            var.supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping


def mergedRecord(first, second, ref_seq):
    """
    Return the VCFRecord corresponding to the merge of first and second.

    :param first: The upstream variant to merge.
    :type first: anacore.vcf.VCFRecord
    :param second: The downstream variant to merge.
    :type second: anacore.vcf.VCFRecord
    :param ref_seq: The sequence of the chromosome.
    :type ref_seq: str
    :return: The variant corresponding to the merge of first and second.
    :rtype: anacore.vcf.VCFRecord
    :todo: Keep INFO and format on strand from FreeBayes, VarDict, ...
    """
    merged = VCFRecord(
        first.chrom,  # chrom
        first.pos,  # pos
        pFormat=first.format
    )
    # Ref and Alt
    first_end = int(round(first.refEnd() - 0.49, 0))
    second_start = int(round(second.refStart() + 0.49, 0))
    ref_add = ""
    if second_start - first_end > 0:
        ref_add = ref_seq[first_end:second_start - 1]
    merged.ref = first.ref + ref_add + second.ref
    merged.ref = merged.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
    merged.alt = [first.alt[0] + ref_add + second.alt[0]]
    merged.alt[0] = merged.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
    # Filter
    first_filters = [] if first.filter is None else first.filter
    second_filters = [] if second.filter is None else second.filter
    merged.filter = list(set(first_filters + second_filters))
    if len(merged.filter) > 1 and "PASS" in merged.filter:
        merged.filter.remove("PASS")
    # Samples
    for spl in first.samples:
        merged.samples[spl] = {}
        if "DP" in first.format:
            merged.samples[spl]["DP"] = min(
                first.getDP(spl), second.getDP(spl)
            )
        if "AD" in first.format:
            fisrt_AD = first.getAD(spl)[0]
            second_AD = second.getAD(spl)[0]
            alt_AD = min(fisrt_AD, second_AD)
            if isinstance(first.samples[spl]["AD"], (int, float)) or len(first.samples[spl]["AD"]) == 1:  # Contains only alt allele
                merged.samples[spl]["AD"] = [alt_AD]
            else:  # Contains ref allele
                # Manage ref_AD also when several variants exist on position
                ref_AD = None
                if fisrt_AD < second_AD:
                    ref_AD = first.samples[spl]["AD"][0]
                elif fisrt_AD > second_AD:
                    ref_AD = second.samples[spl]["AD"][0]
                else:
                    ref_AD = max(first.samples[spl]["AD"][0], second.samples[spl]["AD"][0])
                merged.samples[spl]["AD"] = [ref_AD, alt_AD]
        if "AF" in first.format:
            first_AF = first.getAF(spl)[0]
            second_AF = second.getAF(spl)[0]
            alt_AF = min(first_AF, second_AF)
            if isinstance(first.samples[spl]["AF"], (float, int)) or len(first.samples[spl]["AF"]) == 1:  # Contains only alt allele
                merged.samples[spl]["AF"] = [alt_AF]
            else:  # Contains ref allele
                # Manage ref_AF also when several variants exist on position
                ref_AF = None
                if first_AF < second_AF:
                    ref_AF = first.samples[spl]["AF"][0]
                elif first_AF > second_AD:
                    ref_AF = second.samples[spl]["AF"][0]
                else:
                    ref_AF = max(first.samples[spl]["AF"][0], second.samples[spl]["AF"][0])
                merged.samples[spl]["AF"] = [ref_AF, alt_AF]
    # INFO metrics
    if "AD" in first.info:
        if isinstance(first.info["AD"], (int, float)) or len(first.info["AD"]) == 1:  # Contains only alt allele
            merged.info["AD"] = merged.getPopAD()
        else:
            pop_AD = merged.getPopAD()[0]
            ref_pop_AD = 1 - pop_AD  # Erroneous when several variants exist on position
            if "AD" in merged.format and len(merged.samples.values[0]["AD"]) == 2:
                ref_pop_AD = 0
                for spl_name, spl_info in merged.samples.items():
                    ref_pop_AD += spl_info["AD"][0]
            elif "AF" in merged.format and "DP" in merged.format and len(merged.samples.values[0]["AF"]) == 2:
                ref_pop_AD = 0
                for spl_name, spl_info in merged.samples.items():
                    ref_pop_AD += spl_info["AF"][0] * spl_info["DP"]
            merged.info["AD"] = [int(ref_pop_AD), pop_AD]
    if "DP" in first.info:
        merged.info["DP"] = merged.getPopDP()
    if "AF" in first.info:
        if isinstance(first.info["AF"], (float, int)) or len(first.info["AF"]) == 1:  # Contains only alt allele
            merged.info["AF"] = merged.getPopAF()
        else:
            pop_AF = merged.getPopAF()[0]
            ref_pop_AF = 1 - pop_AF  # Erroneous when several variants exist on position
            if ("AD" in merged.info or "AD" in merged.format) and ("DP" in merged.info or "DP" in merged.format) and len(merged.getPopAD()) == 2:
                ref_pop_AF = round(merged.getPopAD()[0] / merged.getPopDP(), 6)
            merged.info["AF"] = [ref_pop_AF, pop_AF]
    # INFO Parents
    merged.info["MCO_VAR"] = []
    if "MCO_VAR" in first.info:
        for parent in first.info["MCO_VAR"]:
            merged.info["MCO_VAR"].append(parent)
    else:
        merged.info["MCO_VAR"].append(first.getName())
    if "MCO_VAR" in second.info:
        for parent in second.info["MCO_VAR"]:
            merged.info["MCO_VAR"].append(parent)
    else:
        merged.info["MCO_VAR"].append(second.getName())
    # Quality
    merged.info["MCO_QUAL"] = []
    if "MCO_QUAL" in first.info:
        for qual in first.info["MCO_QUAL"]:
            merged.info["MCO_QUAL"].append(qual)
    else:
        merged.info["MCO_QUAL"].append(first.qual)
    if "MCO_QUAL" in second.info:
        for qual in second.info["MCO_QUAL"]:
            merged.info["MCO_QUAL"].append(qual)
    else:
        merged.info["MCO_QUAL"].append(second.qual)
    if None not in merged.info["MCO_QUAL"]:
        merged.qual = mean(merged.info["MCO_QUAL"])
    # Return
    return merged


class LoggerAction(argparse.Action):
    """Manages logger level parameters (The value "INFO" becomes logging.info and so on)."""
    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    """
    :todo:
      - merge by group not upstream to downstream
      - count fragment overlapping the two and not the reads (increase length)
    """
    # Manage parameters
    parser = argparse.ArgumentParser(description='Groups variants occuring in same reads.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-r', '--intersection-rate', default=0.9, type=float, help='Minimum ratio of co-occurancy (nb_reads_containing_the_two_variants / nb_reads_overlapping_the_two_variants_but_containing_only_one). [Default: %(default)s]')
    parser.add_argument('-n', '--intersection-count', default=3, type=int, help='Minimum number of reads containing co-occurancy. [Default: %(default)s]')
    parser.add_argument('-f', '--AF-diff-rate', default=0.2, type=float, help='Maximum difference rate between AF of two merged variants. [Default: %(default)s]')
    parser.add_argument('-d', '--max-distance', default=10, type=int, help='Maximum distance between two merged variants. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='Path to the alignment file (format: BAM).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF). Variants must be ordered by position. They should not be move to upstream.')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to the reference sequences file (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the variant file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Load reference sequences
    log.info("Load reference sequences")
    seq_by_chrom = {}
    with FastaIO(args.input_sequences) as FH_seq:
        for record in FH_seq:
            seq_by_chrom[record.id] = record.string

    # Merge variants
    log.info("Process merge")
    with VCFIO(args.output_variants, "w") as FH_out:
        with pysam.AlignmentFile(args.input_aln, "rb") as FH_aln:
            with VCFIO(args.input_variants) as FH_vcf:
                # Header
                FH_out.copyHeader(FH_vcf)
                FH_out.info["MCO_VAR"] = {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Name of the variants merged because their occur on same reads."}
                FH_out.info["MCO_QUAL"] = {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Qualities of the variants merged because their occur on same reads."}
                FH_out.info["MCO_IR"] = {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Co-occurancy rate between pairs of variants."}
                FH_out.info["MCO_IC"] = {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Co-occurancy count between pairs of variants."}
                FH_out._writeHeader()
                # Records
                prev = None
                for curr in FH_vcf:
                    curr.standardizeSingleAllele()
                    prev_is_merged = False
                    if prev is not None:
                        if prev.chrom == curr.chrom:
                            setRefPos(prev)
                            setRefPos(curr)
                            if prev.start > curr.start:  # Fix order problem after standardization
                                aux = prev
                                prev = curr
                                curr = aux
                            if curr.start - prev.end < args.max_distance:  # The two records are close together
                                prev_AF = prev.getPopAF()[0]
                                curr_AF = curr.getPopAF()[0]
                                AF_diff = 1 - (min(prev_AF, curr_AF) / max(prev_AF, curr_AF))
                                if AF_diff <= args.AF_diff_rate:  # The two records have similar frequencies
                                    prev.isIns = prev.isInsertion()
                                    curr.isIns = curr.isInsertion()
                                    # Set supporting reads
                                    setSupportingReads(prev, curr, seq_by_chrom[prev.chrom], FH_aln, log)
                                    # Check co-occurence
                                    if len(prev.supporting_reads) == 0 and len(curr.supporting_reads) == 0:
                                        log.warning("Nothing read overlapp the two evaluated variants: {} and {}. In this condition the merge cannot be evaluated.".format(prev.getName(), curr.getName()))
                                    else:
                                        intersection_count = len(prev.supporting_reads & curr.supporting_reads)
                                        intersection_rate = intersection_count / len(prev.supporting_reads | curr.supporting_reads)
                                        log.info("{} and {} intersection rate: {:.5} ; number: {}.".format(prev.getName(), curr.getName(), intersection_rate, intersection_count))
                                        if intersection_rate >= args.intersection_rate and intersection_count >= args.intersection_count:
                                            # Merge variants
                                            prev_is_merged = True
                                            merged = mergedRecord(prev, curr, seq_by_chrom[prev.chrom])
                                            log.info("Merge {} and {} in {}.".format(prev.getName(), curr.getName(), merged.getName()))
                                            merged.info["MCO_IR"] = []
                                            if "MCO_IR" in prev.info:
                                                for curr_IR in prev.info["MCO_IR"]:
                                                    merged.info["MCO_IR"].append(curr_IR)
                                            merged.info["MCO_IR"].append(round(intersection_rate, 4))
                                            merged.info["MCO_IC"] = []
                                            if "MCO_IC" in prev.info:
                                                for curr_IC in prev.info["MCO_IC"]:
                                                    merged.info["MCO_IC"].append(curr_IC)
                                            merged.info["MCO_IC"].append(intersection_count)
                                            curr = merged
                        if not prev_is_merged:
                            FH_out.write(prev)
                    prev = curr
                if prev is not None:
                    FH_out.write(prev)
    log.info("End of process.")
