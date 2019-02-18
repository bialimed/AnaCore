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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import time
import random
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.abspath(os.path.join(APP_DIR, "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastaIO
from anacore.sv import HashedSVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
# def fitDistribParamToDup(expected_dup=0.5, distribution_law=stats.exponen):
#     better_param = 0
#     better_observed_diff = math.inf
#     for curr_param in range(0, 5, 0.001):
#         curr_dup = distribution_law.sf(1, curr_param)
#         curr_observed_diff = abs(expected_dup - curr_dup)
#         if curr_observed_diff < better_observed_diff:
#             better_param = curr_param
#             better_observed_diff = curr_observed_diff
#     return better_param

def getNbOccur(in_profile, nb_distinct_reads):
    """
    Return the number of future occurrences for each distinct reads.

    :param in_profile: Path to the file containing the percentage of distinct sequences by number of duplications (format: TSV). Header line must start with "#" and must contain "duplication_level" and "%_distinct".
    :type in_profile: str
    :param nb_distinct_reads: The duplication will be apply on this number of distinct reads.
    :type nb_distinct_reads: int
    :return: The number of future occurrences for each distinct reads.
    :rtype: list
    """
    # Get profile
    profile = None
    with HashedSVIO(in_profile) as FH_in:
        profile = FH_in.read()
    # Get nb_occurences
    nb_occurences = []
    for category in profile:
        category_nb_occur = int(round(category["%_distinct"] * nb_distinct_reads, 0))
        for idx in range(category_nb_occur):
            nb_occurences.append(int(category["duplication_level"]))
    nb_missing = len(nb_occurences) - nb_reads
    for idx in range(nb_missing):
        nb_occurences.append(1)
    # Shuffle nb_occurences
    random.shuffle(nb_occurences)
    return nb_occurences


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', '--random-seed', type=int, default=int(time.time()), help="The seed used for the random generator. If you want reproduce results of one execution: use the same parameters AND the same random-seed. [Default: auto]")
    # parser.add_argument('-d', '--duplication-rate', type=float, default=0.50, help="**************************. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-d', '--duplication-profile', required=True, help='Path to the file containing the percentage of distinct sequences by number of duplications (format: TSV). Header line must start with "#" and must contain "duplication_level" and "%_distinct".')
    group_input.add_argument('-1', '--input-R1', required=True, help='Path to the file containing R1 reads (format: fasta).')
    group_input.add_argument('-2', '--input-R2', required=True, help='Path to the file containing R2 reads (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('--output-R1', required=True, help='Path to the output containing R1 reads and their duplicates (format: fasta).')
    group_output.add_argument('--output-R2', required=True, help='Path to the output containing R2 reads and their duplicates (format: fasta).')
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s -- [%(name)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    logger = logging.getLogger(os.path.basename(__file__))
    logging.info(" ".join(sys.argv))
    logging.info("Random seed used: {}".format(args.random_seed))

    # # Get distribution law parameter to fit % of duplication
    # distribution_law = stats.exponen
    # distrib_shape_param = fitDistribParamToDup(args.duplication_rate, distribution_law)
    # logging.info("Value of shape parameter used to fit distribution law to the duplication rate: {}".format(distrib_shape_param))
    #
    # # Get number of duplications by reads
    # logging.info("Get duplication count for each read")
    # nb_reads = getNbSeq(args.input_R1)
    # nb_occurences = [int(round(elt, 0)) for elt in distribution_law.rvs(distrib_shape_param, size=nb_reads, random_state=args.random_seed)]
    # 0 will become 1

    # Get number of duplications by reads
    logging.info("Get duplication count for each read")
    random.seed(args.random_seed)
    nb_reads = FastaIO.getNbSeq(args.input_R1)
    if nb_reads < 10000:
        logging.error("The number of reads in {} is unsufficient to simulate duplication (found: {} ; expected: {}).".format(args.input_R1, nb_reads, 10000))
    nb_occurences = getNbOccur(args.duplication_profiles, nb_reads)

    # Witre reads
    logging.info("Write reads")
    with FastaIO(args.output_R1, "w") as FH_out_R1:
        with FastaIO(args.output_R2, "w") as FH_out_R2:
            with FastaIO(args.input_R1) as FH_in_R1:
                with FastaIO(args.input_R1) as FH_in_R2:
                    for curr_nb_occur, R1, R2 in zip(nb_occurences, FH_in_R1, FH_in_R2):
                        R1.description += "_dupCount:{}".format(curr_nb_occur)
                        R2.description += "_dupCount:{}".format(curr_nb_occur)
                        old_R1_id = R1.id
                        for idx in range(curr_nb_occur):
                            R1.id = old_R1_id + "_dupId:{}".format(idx)
                            R2.id = old_R1_id + "_dupId:{}".format(idx)
                            FH_out_R1.write(R1)
                            FH_out_R2.write(R2)

    logging.info("End of process")
