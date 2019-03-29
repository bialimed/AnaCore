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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse
from anacore.vcf import VCFIO, VCFRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAlnCmp(read):
    ref_aln = []
    read_aln = []
    ref_seq = [elt for elt in read.get_reference_sequence()]  # reference sequence on alignment
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


def setSupportingReads(prev, curr, FH_aln):
    """
    set supporting reads for each of two variants with fragment overlapping the two.
    """
    prev.start = int(prev.refStart())
    curr.end = int(curr.refEnd())
    prev.supporting_reads = set()
    curr.supporting_reads = set()
    for read in FH_aln.fetch(curr.chrom, prev.start, curr.end):
        if not read.is_duplicate:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                overlap_the_two = (ref_start <= prev.end and ref_end >= prev.start) and (ref_start <= curr.end and ref_end >= curr.start)
                if overlap_the_two:
                    ref_aln, read_aln = getAlnCmp(read)
                    for var in [prev, curr]:
                        alt = read_aln[var.start - ref_start:var.end - ref_start + 1]
                        ref = ref_aln[var.start - ref_start:var.end - ref_start + 1]
                        if var.isIns:
                            if len(ref[0]) != len(alt[0]):  # Read contains an insertion
                                while len(ref[0]) > 0 and len(alt[0]) > 0 and alt[0][0] == ref[0][0]:
                                    alt[0] = alt[0][1:]
                                    ref[0] = ref[0][1:]
                        var_alt = var.alt[0].upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                        var_ref = var.ref.upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                        if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # Is present
                            log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                            var.supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping


def mergeRecords(first, second, ref_seq):
    # ref
    # alt
    # known
    # pos
    # AF
    # AD
    # DP
    # Filter ?
    # INFO add intersection rate and interstion nb
    pass


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
    # Manage parameters
    parser = argparse.ArgumentParser(description='Groups variants occuring in same reads.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-r', '--intersection-rate', default=0.98, type=float, help='****************. [Default: %(default)s]')
    parser.add_argument('-n', '--intersection-nb', default=3, type=int, help='****************. [Default: %(default)s]')
    parser.add_argument('-f', '--af-diff-threshold', default=0.01, type=float, help='**********************. [Default: %(default)s]')
    parser.add_argument('-d', '--max-distance', default=15, type=int, help='**********************. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='Path to the alignment file (format: BAM).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the varuiants file (format: VCF).')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to the reference sequences file (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the variant file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with VCFIO(args.output_variants, "w") as FH_out:
        with pysam.AlignmentFile(args.input_aln, "rb") as FH_aln:
            with VCFIO(args.input_variants) as FH_vcf:
                # Header
                FH_out.copyHeader(FH_vcf)
                FH_out._writeHeader()
                # Records
                prev = None  #################################### pb quand plusieurs
                for curr in FH_vcf:
                    prev_is_merged = False
                    if prev is not None:
                        prev.end = int(prev.refEnd())
                        curr.start = int(curr.refStart())
                        if prev.chrom == curr.chrom and prev.end - curr.start < args.max_distance:  # The two records are close together
                            prev_AF = prev.getPopAF()[0]
                            curr_AF = curr.getPopAF()[0]
                            if min(prev_AF, curr_AF) / max(prev_AF, curr_AF) >= args.af_diff_threshold:  # The two records have similar frequencies
                                prev.isIns = prev.isInsertion()
                                curr.isIns = curr.isInsertion()
                                # Set supporting reads
                                setSupportingReads(prev, curr, FH_aln)
                                # Check co-occurence
                                if len(prev.supporting_reads) == 0 and len(curr.supporting_reads) == 0:
                                    log.warnings("Nothing read overlapp the two evaluated variants: {} and {}. In this condition the merge cannot be evaluated.".format(prev.getName(), curr.getName()))
                                else:
                                    intersection_nb = len(prev.supporting_reads & curr.supporting_reads)
                                    intersection_rate = intersection_nb / len(prev.supporting_reads | curr.supporting_reads)
                                    log.info("{} and {} intersection rate: {:.5} ; number: {}.".format(prev.getName(), curr.getName(), intersection_rate, intersection_nb))
                                if intersection_rate >= args.intersection_rate and intersection_nb >= args.intersection_nb:
                                    # Merge variants
                                    log.info("Merge {} and {}.".format(prev.getName(), curr.getName()))
                                    prev_is_merged = True
                                    mergeRecords(prev, curr, chrom_seq)
                                    ############################ merge records (complete with seq ref)
                                    ## needs seq to complete
                                    ## needs to take into account break in intron
                        if not prev_is_merged:
                            FH_out.write(prev)
                    prev = curr
                if prev is not None:
                    FH_out.write(prev)
    log.info("End of process.")
