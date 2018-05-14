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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import Sequence, FastqIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def writeReport(combined, R1, out_report):
    """
    @summary: Writes report file for combination results.
    @param combined: [str] Path to the file containing the combined reads (format: fastq).
    @param R1: [str] Path to the initial R1 file (format: fastq).
    @param out_report: [str] Path to the outputted report file (format: json).
    """
    report = {
        "nb_combined_pairs": 0,
        "nb_uncombined_pairs": 0,
        "nb_by_lengths": dict()
    }
    # Get nb combined and lengths distribution
    with FastqIO(combined) as FH_comb:
        for record in FH_comb:
            report["nb_combined_pairs"] += 1
            curr_len = len(record.string)
            if curr_len not in report["nb_by_lengths"]:
                report["nb_by_lengths"][curr_len] = 1
            else:
                report["nb_by_lengths"][curr_len] += 1
    # Get nb uncombined
    nb_total_pairs = 0
    with FastqIO(R1) as FH_R1:
        for record in FH_R1:
            nb_total_pairs += 1
    report["nb_uncombined_pairs"] = nb_total_pairs - report["nb_combined_pairs"]
    # Write report
    with open(out_report, "w") as FH_report:
        json.dump(report, FH_report, sort_keys=True)


def nucRevCom(seq):
    """
    @summary: Returns the reverse complementent of the sequence.
    @param seq: [str] The sequence to process.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return "".join([complement_rules[base] for base in seq[::-1]])


def seqRevCom(seq):
    """
    @summary: Returns the reverse complement of the object sequence.
    @param seq: [Sequence] The sequence to process.
    @return: [Sequence] The reverse complement of the object sequence.
    """
    return Sequence(
        seq.id,
        nucRevCom(seq.string),  # Reverse complement sequence
        seq.description,
        seq.quality[::-1]  # Reverse quality
    )


def process(args, log):
    """
    @summary: Combines R1 and R2 by their overlapping segment.
    param args: [Namespace] The namespace extract from the script arguments.
    param log: [Logger] The logger of the script.
    """
    nb_pairs = 0
    combined = 0
    with FastqIO(args.output_combined, "w") as FH_combined:
        with FastqIO(args.input_R1) as FH_r1:
            with FastqIO(args.input_R2) as FH_r2:
                for R1 in FH_r1:
                    R2 = seqRevCom(FH_r2.next_seq())
                    nb_pairs += 1
                    best_overlap = None
                    max_nb_support = -1
                    R1_len = len(R1.string)
                    R2_len = len(R2.string)
                    R1_start = 0
                    R2_start = R2_len - args.min_overlap
                    is_valid = R1_len >= args.min_overlap and R2_len >= args.min_overlap
                    can_be_better = True
                    while is_valid and can_be_better:  # For each shift
                        nb_support = 0
                        nb_contradict = 0
                        consensus_seq = ""
                        consensus_qual = ""
                        curr_overlap_len = min(R1_len - R1_start, R2_len - R2_start)
                        if best_overlap is not None and R1_start != 0 and curr_overlap_len < best_overlap["nb_support"]:  # R1 is first and overlap become lower than nb support
                            can_be_better = False
                        else:
                            # Evaluate overlap
                            R1_ov_s = R1.string[R1_start:R1_start + curr_overlap_len]
                            R1_ov_q = R1.quality[R1_start:R1_start + curr_overlap_len]
                            R2_ov_s = R2.string[R2_start:R2_start + curr_overlap_len]
                            R2_ov_q = R2.quality[R2_start:R2_start + curr_overlap_len]
                            for nt_R1, qual_R1, nt_R2, qual_R2 in zip(R1_ov_s, R1_ov_q, R2_ov_s, R2_ov_q):  # For each nt in overlap
                                if nt_R1 == nt_R2:
                                    nb_support += 1
                                    consensus_seq += nt_R1
                                    consensus_qual += max(qual_R1, qual_R2)
                                else:
                                    nb_contradict += 1
                                    if qual_R1 >= qual_R2:
                                        consensus_seq += nt_R1
                                        consensus_qual += qual_R1
                                    else:
                                        consensus_seq += nt_R2
                                        consensus_qual += qual_R2
                            # Filter consensus and select the best
                            if nb_support >= max_nb_support:
                                if float(nb_contradict)/curr_overlap_len <= args.max_contradict_ratio:
                                    max_nb_support = nb_support
                                    best_overlap = {
                                        "nb_support": nb_support,
                                        "nb_contradict": nb_contradict,
                                        "consensus_seq": consensus_seq,
                                        "consensus_qual": consensus_qual,
                                        "R1_start": R1_start,
                                        "R2_start": R2_start,
                                        "length": curr_overlap_len
                                    }
                            # Next shift
                            if R1_start == 0:
                                R2_start -= 1
                                if R2_start == -1:
                                    R2_start = 0
                                    R1_start = 1
                            else:
                                R1_start += 1
                                if R1_len - R1_start < args.min_overlap:
                                    is_valid = False
                    if best_overlap is not None:  # Current pair has valid combination
                        combined += 1
                        complete_seq = best_overlap["consensus_seq"]
                        complete_qual = best_overlap["consensus_qual"]
                        if best_overlap["R1_start"] > 0:  # If R1 start before R2 (insert size > read length)
                            complete_seq = R1.string[0:best_overlap["R1_start"]] + complete_seq + R2.string[best_overlap["length"]:]
                            complete_qual = R1.quality[0:best_overlap["R1_start"]] + complete_qual + R2.quality[best_overlap["length"]:]
                        consensus_record = Sequence(
                            R1.id,
                            complete_seq,
                            "Support_ratio:{}/{};R1_start:{};R2_start:{}".format(
                                best_overlap["nb_support"],
                                best_overlap["length"],
                                best_overlap["R1_start"],
                                best_overlap["R2_start"]
                            ),
                            complete_qual
                        )
                        FH_combined.write(consensus_record)
    # Log
    log.info(
        "Nb pair: {} ; Nb combined: {} ({}%)".format(
            nb_pairs,
            combined,
            (0 if nb_pairs == 0 else round(float(combined*100)/nb_pairs, 2))
        )
    )
    if args.output_report is not None:
        writeReport(args.output_combined, args.input_R1, args.output_report)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Combines R1 and R2 by their overlapping segment.')
    parser.add_argument('-o', '--min-overlap', default=20, type=int, help='Minimum overlap between R1 and R2. [Default: %(default)s]')
    parser.add_argument('-m', '--max-contradict-ratio', default=0.1, type=float, help='Error ratio in overlap region between R1 and R2. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-1', '--input-R1', required=True, help='The path to the R1 file (format: fastq).')
    group_input.add_argument('-2', '--input-R2', required=True, help='The path to the R2 file (format: fastq).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-c', '--output-combined', required=True, help='The path to the file with combined pairs (format: fastq).')
    group_output.add_argument('-r', '--output-report', help='The path to the path containing combination metrics (format: JSON).')
    args = parser.parse_args()

    # Process
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger("combinePairs")
    log.setLevel(logging.INFO)
    log.info("Start")
    log.info("Command: " + " ".join(sys.argv))
    process(args, log)
    log.info("End")
