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
__status__ = 'dev'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from sequenceIO import Sequence, FastqIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def nucRevCom(seq):
    """
    @summary: Returns the reverse complementent of the sequence.
    @param seq: [str] The sequence to process.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N'}
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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Combines R1 and R2 by their overlapping segment when sum of read length in pair is superior than the segment length.')
    parser.add_argument('-o', '--min-overlap', default=20, type=int, help='Minimum overlap between R1 and R2. [Default: %(default)s]')
    parser.add_argument('-r', '--max-contradict-ratio', default=0.1, type=float, help='Error ratio in overlap region between R1 and R2. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-1', '--input-R1', required=True, help='The path to the R1 file (format: fastq).')
    group_input.add_argument('-2', '--input-R2', required=True, help='The path to the R2 file (format: fastq).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-c', '--output-combined', help='The path of the file with combined pairs (format: fastq).')
    args = parser.parse_args()

    # Process
    nb_pairs = 0
    combined = 0

    with FastqIO(args.output_combined, "w") as FH_combined:
        with FastqIO(args.input_R1) as FH_r1:
            with FastqIO(args.input_R2) as FH_r2:
                for R1 in FH_r1:
                    R2 = seqRevCom(FH_r2.next_seq())
                    nb_pairs += 1
                    consensus_record = None
                    max_nb_support = -1
                    R1_len = len(R1.string)
                    R2_len = len(R2.string)
                    R1_start = 0
                    R2_start = R2_len - args.min_overlap
                    is_valid = R1_len >= args.min_overlap and R2_len >= args.min_overlap
                    while is_valid:
                        nb_support = 0
                        nb_contradict = 0
                        consensus_seq = ""
                        consensus_qual = ""
                        overlap_len = min(R1_len - R1_start, R2_len - R2_start)
                        for idx in range(overlap_len):
                            nt_R1 = R1.string[R1_start + idx]
                            qual_R1 = R1.quality[R1_start + idx]
                            nt_R2 = R2.string[R2_start + idx]
                            qual_R2 = R2.quality[R2_start + idx]
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
                            if float(nb_contradict)/(nb_contradict+nb_support) <= args.max_contradict_ratio:
                                max_nb_support = nb_support
                                all_seq = consensus_seq
                                all_qual = consensus_qual
                                if R1_start > 0:  # If R1 start before R2 (insert size > read length)
                                    all_seq = R1.string[0:R1_start] + consensus_seq + R2.string[overlap_len:]
                                    all_qual = R1.quality[0:R1_start] + consensus_qual + R2.quality[overlap_len:]
                                consensus_record = Sequence(
                                    R1.id,
                                    all_seq,
                                    "Support_ratio:{}/{};R1_start:{};R2_start:{}".format(nb_support, nb_contradict+nb_support, R1_start, R2_start),
                                    all_qual
                                )
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
                    if consensus_record is not None:
                        combined += 1
                        FH_combined.write(consensus_record)
    print("Nb pair: {}\nNb combined: {} ({}%)".format(combined, nb_pairs, round(float(combined*100)/nb_pairs, 2)))
