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
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse
from copy import deepcopy
from statistics import mean

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import IdxFastaIO
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr


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


def setRefPos(variant, seq_handler, padding=200):
    """
    Add start and end attributes in VCFRecord. For insertions the start is defined on the first position before the insertion and the end on the last position affected by the insertion.

    :param variant: The variant to update.
    :type variant: anacore.vcf.VCFRecord
    """
    if variant.ref == VCFRecord.getEmptyAlleleMarker() or variant.alt[0] == VCFRecord.getEmptyAlleleMarker():  # Normalized indel
        # Most upstream
        variant.upstream_start, variant.upstream_end = getStartEnd(variant)
        # Most downstream
        sub_region = seq_handler.getSub(variant.chrom, variant.pos - 2, variant.pos + len(variant.ref) + padding)
        chrom_pos = variant.pos
        variant.pos = 3  # Switch position from chromosome to position from subregion
        downstream_var = variant.getMostDownstream(sub_region)
        variant.pos = chrom_pos + variant.pos - 3  # Switch position from subregion to position from chromosome
        downstream_var.pos = variant.pos
        variant.downstream_start, variant.downstream_end = getStartEnd(downstream_var)
    else:
        variant.upstream_start, variant.upstream_end = getStartEnd(variant)
        variant.downstream_start = variant.upstream_start
        variant.downstream_end = variant.upstream_end


def getStartEnd(variant):
    """
    Return reference start and end for a variant. For insertions the start is defined on the first position before the insertion and the end on the last position affected by the insertion.

    :param variant: The variant.
    :type variant: anacore.vcf.VCFRecord
    :return: start and end.
    :rtype: tuple
    """
    start = int(variant.refStart())
    if variant.isInsertion():
        if start == int(variant.refStart() + 0.5):
            start -= 1
    end = int(variant.refEnd())
    return start, end


def getReadRefAlt(ref_aln, read_aln, ref_start, target_is_ins, target_start, target_end):
    """
    Return reference sequence and read sequence for the selected target.

    :param ref_aln: Reference sequence in alignment.
    :type ref_aln: str
    :param read_aln: Read sequence in alignment.
    :type read_aln: str
    :param ref_start: Start position for the alignment on the reference (1-based).
    :type ref_start: int
    :param target_is_ins: The anlyzed target correspond to an exon.
    :type target_is_ins: bool
    :param target_start: Start reference position for the target (1-based).
    :type target_start: int
    :param target_end: End reference position for the target (1-based).
    :type target_end: int
    :return: Reference and alternative sequence for the read.
    :rtype: tuple
    """
    alt = read_aln[target_start - ref_start:target_end - ref_start + 1]
    ref = ref_aln[target_start - ref_start:target_end - ref_start + 1]
    if target_is_ins:
        while len(ref[0]) > 0 and len(alt[0]) > 0 and alt[0][0] == ref[0][0]:
            alt[0] = alt[0][1:]
            ref[0] = ref[0][1:]
    return ref, alt


def getIncludingReads(FH_aln, chrom_id, target_start, target_end):
    """
    Return read ID of reads including the target.

    :param FH_aln: The file handle to the alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param chrom_id: Chromosome ID.
    :type chrom_id: str
    :param target_start: Start position for target.
    :type target_start: int
    :param target_end: End position for target.
    :type target_end: int
    :return: Reads IDs of reads including the target.
    :rtype: set
    """
    including_reads = set()
    for read in FH_aln.fetch(chrom_id, target_start - 1, target_end):
        if not read.is_duplicate:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                includes = (ref_start <= target_start and ref_end >= target_end)
                if includes:
                    including_reads.add(read.query_name)
    return including_reads


def getSupportingReads(var, chrom_seq, FH_aln, log):
    """
    Return read ID of reads supporting the altenative variant.

    :param var: The variant.
    :type var: anacore.vcf.VCFRecord updated with setRefPos() and isIns
    :param chrom_seq: The sequence of the chromosome.
    :type chrom_seq: str
    :param FH_aln: The file handle to the alignments file. The variants must have been defined from this alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param log: The logger object.
    :type log: logging.Logger
    :return: The list of supporting reads IDs.
    :rtype: set
    """
    supporting_reads = set()
    is_insertion = var.isInsertion()
    for read in FH_aln.fetch(var.chrom, var.upstream_start - 1, var.downstream_end):
        if not read.is_duplicate:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                overlap_var = (ref_start <= var.upstream_start and ref_end >= var.downstream_end)
                if overlap_var:
                    ref_aln, read_aln = getAlnCmp(read, chrom_seq[ref_start - 1:ref_end])
                    var_alt = var.alt[0].upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                    var_ref = var.ref.upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                    # Test with upstream coordinates
                    ref, alt = getReadRefAlt(ref_aln, read_aln, ref_start, is_insertion, var.upstream_start, var.upstream_end)
                    if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # The alternative is present on most upstream coordinates
                        log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                        supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping
                    # Test with downstream coordinates
                    elif var.upstream_start != var.downstream_start:
                        ref, alt = getReadRefAlt(ref_aln, read_aln, ref_start, is_insertion, var.downstream_start, var.downstream_end)
                        if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # The alternative is present on most downstream coordinates
                            log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                            supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping
    return supporting_reads


def areColocated(first, second):
    """
    Return True if one of the two variants is included in the location of the other.

    :param first: The variant.
    :type first: anacore.vcf.VCFRecord updated with setRefPos()
    :param second: The variant.
    :type second: anacore.vcf.VCFRecord updated with setRefPos()
    :return: True if one of the two variants is included in the location of the other.
    :rtype: bool
    """
    first_up = set([elt for elt in range(first.upstream_start, first.upstream_end + 1)])
    second_up = set([elt for elt in range(second.upstream_start, second.upstream_end + 1)])
    first_down = set([elt for elt in range(first.downstream_start, first.downstream_end + 1)])
    second_down = set([elt for elt in range(second.downstream_start, second.downstream_end + 1)])
    included = False
    if len(first_up - second_up) == 0:
        included = True
    elif len(first_up - second_down) == 0:
        included = True
    elif len(first_down - second_up) == 0:
        included = True
    elif len(first_down - second_down) == 0:
        included = True
    elif len(second_up - first_up) == 0:
        included = True
    elif len(second_down - first_up) == 0:
        included = True
    elif len(second_up - first_down) == 0:
        included = True
    elif len(second_down - first_down) == 0:
        included = True
    return included


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
            merged.info["AD"] = merged.getPopAltAD()
        else:
            merged.info["AD"] = [merged.getPopRefAD()] + merged.getPopAltAD()
    if "DP" in first.info:
        merged.info["DP"] = merged.getPopDP()
    if "AF" in first.info:
        if isinstance(first.info["AF"], (float, int)) or len(first.info["AF"]) == 1:  # Contains only alt allele
            merged.info["AF"] = merged.getPopAltAF()
        else:
            merged.info["AF"] = [merged.getPopRefAF()] + merged.getPopAltAF()
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
    :todo: count fragment overlapping the two and not the reads (increase length)
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
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF). Variants must be ordered by position and should be move to upstream.')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to the reference sequences file (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the variant file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Merge variants
    with IdxFastaIO(args.input_sequences, use_cache=True) as FH_seq:
        with VCFIO(args.output_variants, "w") as FH_out:
            with pysam.AlignmentFile(args.input_aln, "rb") as FH_aln:
                with VCFIO(args.input_variants) as FH_vcf:
                    # Header
                    FH_out.copyHeader(FH_vcf)
                    FH_out.info["MCO_VAR"] = HeaderInfoAttr("MCO_VAR", "Name of the variants merged because their occur on same reads.", type="String", number=".")
                    FH_out.info["MCO_QUAL"] = HeaderInfoAttr("MCO_QUAL", "Qualities of the variants merged because their occur on same reads.", type="String", number=".")
                    FH_out.info["MCO_IR"] = HeaderInfoAttr("MCO_IR", "Co-occurancy rate between pairs of variants.", type="String", number=".")
                    FH_out.info["MCO_IC"] = HeaderInfoAttr("MCO_IC", "Co-occurancy count between pairs of variants.", type="String", number=".")
                    FH_out.writeHeader()
                    # Records
                    prev_list = list()
                    for curr in FH_vcf:
                        chrom_seq = FH_seq.get(curr.chrom).string
                        std_curr = deepcopy(curr)
                        curr.normalizeSingleAllele()
                        curr.supporting_reads = None
                        setRefPos(curr, FH_seq)
                        merged_idx = set()
                        removed_idx = None
                        for idx, (prev, std_prev) in enumerate(prev_list[::-1]):
                            if removed_idx is None:  # Distance between current variant and previous most upstream variants is ok
                                if prev.chrom != curr.chrom:
                                    removed_idx = len(prev_list) - 1 - idx
                                else:
                                    variants_distance = max(
                                        max(0, curr.upstream_start - prev.downstream_end),  # prev is before
                                        max(0, prev.upstream_start - curr.downstream_end)  # prev is after
                                    )
                                    if variants_distance > args.max_distance:  # The two records are close together
                                        removed_idx = len(prev_list) - 1 - idx
                                    elif areColocated(curr, prev):
                                        log.debug("Skip colocated variants {} and {}.".format(prev.getName(), curr.getName()))
                                    else:  # The two records are close together
                                        prev_AF = prev.getPopAltAF()[0]
                                        curr_AF = curr.getPopAltAF()[0]
                                        AF_diff = 1 - (min(prev_AF, curr_AF) / max(prev_AF, curr_AF))
                                        log.debug("Allelels frequencies for {} and {}: {:.1%} and {:.1%} (diff rate: {:.2}).".format(prev.getName(), curr.getName(), prev_AF, curr_AF, AF_diff))
                                        if AF_diff <= args.AF_diff_rate:  # The two records have similar frequencies
                                            # Set supporting reads
                                            if prev.supporting_reads is None:
                                                prev.supporting_reads = getSupportingReads(prev, chrom_seq, FH_aln, log)
                                            if curr.supporting_reads is None:
                                                curr.supporting_reads = getSupportingReads(curr, chrom_seq, FH_aln, log)
                                            shared_reads = getIncludingReads(
                                                FH_aln,
                                                curr.chrom,
                                                min(prev.upstream_start, curr.upstream_start),
                                                max(prev.downstream_end, curr.downstream_end)
                                            )
                                            # Check co-occurence
                                            if len(shared_reads) == 0:
                                                log.warning("Nothing read overlapp the two evaluated variants: {} and {}. In this condition the merge cannot be evaluated.".format(prev.getName(), curr.getName()))
                                            else:
                                                prev_support_shared = (prev.supporting_reads & shared_reads)
                                                curr_support_shared = (curr.supporting_reads & shared_reads)
                                                intersection_count = len(prev_support_shared & curr_support_shared)
                                                intersection_rate = intersection_count / len(prev_support_shared | curr_support_shared)
                                                log.debug("{} and {} intersection rate: {:.5} ; number: {}.".format(prev.getName(), curr.getName(), intersection_rate, intersection_count))
                                                if intersection_rate >= args.intersection_rate and intersection_count >= args.intersection_count:
                                                    # Merge variants
                                                    first = prev
                                                    second = curr
                                                    if first.upstream_start > second.upstream_start:
                                                        first = curr
                                                        second = prev
                                                    merged = mergedRecord(first, second, chrom_seq)
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
                                                    merged.normalizeSingleAllele()
                                                    merged.supporting_reads = None
                                                    setRefPos(merged, chrom_seq)
                                                    curr = merged
                                                    std_curr = deepcopy(curr)
                                                    std_curr.fastStandardize(FH_seq, 200)
                                                    merged_idx.add(len(prev_list) - 1 - idx)
                        # Write the far records and remove them from the previous ones
                        if removed_idx is None:
                            for idx in merged_idx:
                                del(prev_list[idx])
                        else:
                            for idx in merged_idx:
                                if idx > removed_idx:
                                    del(prev_list[idx])
                            for idx in range(removed_idx + 1):
                                record, std_record = prev_list.pop()
                                if idx not in merged_idx:
                                    FH_out.write(std_record)
                        prev_list.append((curr, std_curr))
                    for record, std_record in prev_list:
                        FH_out.write(std_record)
    log.info("End of job")
