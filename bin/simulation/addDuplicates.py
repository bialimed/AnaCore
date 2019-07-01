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
        nb_reads_at_dup_lvl = int(round(float(category["%_distinct"]) * nb_distinct_reads / 100, 0))
        for idx in range(nb_reads_at_dup_lvl):
            nb_occurences.append(int(category["duplication_level"]))
    nb_missing = len(nb_occurences) - nb_distinct_reads
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
    parser = argparse.ArgumentParser(description='Add duplicated sequences in reads files. The duplication model (number of duplication for each distinct sequence) follow the profile provided by user.')
    parser.add_argument('-s', '--random-seed', type=int, default=int(time.time()), help="The seed used for the random generator. If you want reproduce results of one execution: use the same parameters AND the same random-seed. [Default: auto]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-d', '--duplication-profile', required=True, help='Path to the file containing the percentage of distinct sequences by number of duplications (format: TSV). Header line must start with "#" and must contain "duplication_level" and "%%_distinct".')
    group_input.add_argument('-1', '--input-R1', required=True, help='Path to the file containing R1 reads (format: fasta).')
    group_input.add_argument('-2', '--input-R2', required=True, help='Path to the file containing R2 reads (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('--output-R1', required=True, help='Path to the output containing R1 reads and their duplicates (format: fasta).')
    group_output.add_argument('--output-R2', required=True, help='Path to the output containing R2 reads and their duplicates (format: fasta).')
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info(" ".join(sys.argv))
    log.info("Random seed used: {}".format(args.random_seed))

    # Get number of duplications by reads
    log.info("Get duplication count for each read")
    random.seed(args.random_seed)
    nb_reads = FastaIO.nbSeq(args.input_R1)
    if nb_reads < 10000:
        log.error("The number of reads in {} is unsufficient to simulate duplication (found: {} ; expected: {}).".format(args.input_R1, nb_reads, 10000))
    nb_occurences = getNbOccur(args.duplication_profile, nb_reads)

    # Witre reads
    log.info("Write reads")
    with FastaIO(args.output_R1, "w") as FH_out_R1:
        with FastaIO(args.output_R2, "w") as FH_out_R2:
            with FastaIO(args.input_R1) as FH_in_R1:
                with FastaIO(args.input_R2) as FH_in_R2:
                    for curr_nb_occur, R1, R2 in zip(nb_occurences, FH_in_R1, FH_in_R2):
                        description = "dupCount={}".format(curr_nb_occur)
                        if R1.description is not None and R1.description != "":
                            description = R1.description + "_" + description
                        R1.description = description
                        description = "dupCount={}".format(curr_nb_occur)
                        if R2.description is not None and R2.description != "":
                            description = R2.description + "_" + description
                        R2.description = description
                        old_R1_id = R1.id
                        for idx in range(curr_nb_occur):
                            R1.id = old_R1_id + "_dupId={}".format(idx)
                            R2.id = old_R1_id + "_dupId={}".format(idx)
                            FH_out_R1.write(R1)
                            FH_out_R2.write(R2)

    log.info("End of job")
