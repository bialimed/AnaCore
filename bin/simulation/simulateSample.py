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
import logging
import argparse
import subprocess


########################################################################
#
# FUNCTIONS
#
########################################################################
def getLogLevelParam(str_level):
    """
    Return the level parameter for logging.basicConfig from the string version of the level.

    :param str_level: The logging level.
    :type str_level: str
    :return: The parameter for logging.basicConfig.
    :rtype: int
    """
    log_level = None
    if str_level == "DEBUG":
        log_level = logging.DEBUG
    elif str_level == "INFO":
        log_level = logging.INFO
    elif str_level == "WARNING":
        log_level = logging.WARNING
    elif str_level == "ERROR":
        log_level = logging.ERROR
    elif str_level == "CRITICAL":
        log_level = logging.CRITICAL
    else:
        raise ValueError("The value {} is invalid for logging level.".format(str_level))
    return log_level


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Simulate enrichment sample.')
    parser.add_argument('--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help='The logger level. [Default: %(default)s]')
    parser.add_argument('--random-seed', type=int, default=int(time.time()), help="The seed used for the random generator. If you want reproduce results of one execution: use the same parameters AND the same random-seed. [Default: auto]")
    parser.add_argument('--nb-random-simulations', type=int, default=30, help='Number of simulations used to find the best reads set to optimize variants frequencies simulated compared to expected. [Default: %(default)s]')
    parser.add_argument('--min-distance', type=int, default=3, help="The minimum distance between two variants. [Default: %(default)s]")
    parser.add_argument('--targets-padding', type=int, default=150, help="Padding introduced by panel around the primary targets. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_quality = parser.add_argument_group('Quality')  # Quality
    group_quality.add_argument('--qual-offset', type=int, default=33, help="The position of the first quality encoding character in ASCII table (example: 33 for Illumina 1.8+). [Default: %(default)s]")
    group_quality.add_argument('--qual-penalty', type=int, default=0, help="The penalty applied to reduce the quality of the sequences produced. With 2, the quality of each base is decrease of 2 compared to the model. [Default: %(default)s]")
    group_reads = parser.add_argument_group('Reads')  # Reads
    group_reads.add_argument('--R1-end-adapter', default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGCGGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA", help='The sequence of the Illumina p7 adapter + flowcell anchor (stretch of 10 A). It is found at the end of the R1 when the read length is superior to the fragment length. [Default: %(default)s]')
    group_reads.add_argument('--R2-end-adapter', default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTGACAAGCGCTTGTCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAA", help='The reverse complemented sequence of the Illumina p5 adapter + flowcell anchor (stretch of 10 A). It is found at the end of the R2 when the read length is superior to the fragment length. [Default: %(default)s]')
    group_reads.add_argument('--reads-length', type=int, default=150, help='The reads length. [Default: %(default)s]')
    group_reads.add_argument('--fragments-length', type=int, default=180, help='The fragments lengths. [Default: %(default)s]')
    group_reads.add_argument('--fragments-length-sd', type=int, default=30, help='The standard deviation for fragments lengths. [Default: %(default)s]')
    group_depths = parser.add_argument_group('Depths')  # Depths
    group_depths.add_argument('--min-depth', type=int, default=200, help='The minimum depth on targets. [Default: %(default)s]')
    group_depths.add_argument('--max-depth', type=int, default=15000, help='The maximum depth on targets. [Default: %(default)s]')
    group_depths.add_argument('--min-add-depth', type=int, default=5, help='When the depth increase between two consecutive positions it increase to this minimum number of reads. [Default: %(default)s]')
    group_depths.add_argument('--max-add-depth', type=int, default=15, help="When the depth increase between two consecutive positions it can't increase more than this number of reads. [Default: %(default)s]")
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('--input-duplication-profile', required=True, help='Path to the file containing the percentage of distinct sequences by number of duplications (format: TSV). Header line must start with "#" and must contain "duplication_level" and "%%_distinct".')
    group_input.add_argument('--input-targets', required=True, help='Path to the targets (format: BED).')
    group_input.add_argument('--input-reference', required=True, help='Path to the reference sequences file (format: fasta).')
    group_input.add_argument('--input-variants-profile', required=True, help='Path to the variants profile file (format: TSV).')
    group_input.add_argument('--input-R1-quality-models', required=True, nargs='+', help='The paths of the sequences files used to retrieve the error model for R1 (format: fastq). The sequences length must be at least the same as sequences provided by --reads-length.')
    group_input.add_argument('--input-R2-quality-models', required=True, nargs='+', help='The paths of the sequences files used to retrieve the error model for R2 (format: fastq). The sequences length must be at least the same as sequences provided by --reads-length.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('--output-R1', required=True, help='Path to the output containing R1 reads (format: fastq).')
    group_output.add_argument('--output-R2', required=True, help='Path to the output containing R2 reads (format: fatsq).')
    group_output.add_argument('--output-variants', required=True, help='Path to the file describing variants introduced in reads and their AF (format: VCF).')
    group_output.add_argument('--output-errors-R1', required=True, help='The path of the file reporting the sequencing errors introduced on R1 (format: TSV).')
    group_output.add_argument('--output-errors-R2', required=True, help='The path of the file reporting the sequencing errors introduced on R2 (format: TSV).')
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(level=getLogLevelParam(args.logging_level), format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info(" ".join(sys.argv))
    log.info("Random seed used: {}".format(args.random_seed))

    # Temporary files
    tmp_sim_R1 = args.output_R1 + "_tmp.fasta.gz"
    tmp_sim_R2 = args.output_R2 + "_tmp.fasta.gz"
    tmp_dup_R1 = args.output_R1 + "_tmp_dup.fasta.gz"
    tmp_dup_R2 = args.output_R2 + "_tmp_dup.fasta.gz"

    # Simulate reads
    cmd = [
        "simuReads.py",
        "--logging-level", args.logging_level,
        "--random-seed", args.random_seed,
        "--nb-random-simulations", args.nb_random_simulations,
        "--min-distance", args.min_distance,
        "--targets-padding", args.targets_padding,
        "--R1-end-adapter", args.R1_end_adapter,
        "--R2-end-adapter", args.R2_end_adapter,
        "--reads-length", args.reads_length,
        "--fragments-length", args.fragments_length,
        "--fragments-length-sd", args.fragments_length_sd,
        "--min-depth", args.min_depth,
        "--max-depth", args.max_depth,
        "--min-add-depth", args.min_add_depth,
        "--max-add-depth", args.max_add_depth,
        "--input-targets", args.input_targets,
        "--input-reference", args.input_reference,
        "--input-profile", args.input_variants_profile,
        "--output-variants", args.output_variants,
        "--output-R1", tmp_sim_R1,
        "--output-R2", tmp_sim_R2
    ]
    cmd = [str(elt) for elt in cmd]
    log.info("Simulate reads: {}".format(" ".join(cmd)))
    subprocess.check_call(cmd)

    # Add duplicates
    cmd = [
        "addDuplicates.py",
        "--random-seed", args.random_seed,
        "--duplication-profile", args.input_duplication_profile,
        "--input-R1", tmp_sim_R1,
        "--input-R2", tmp_sim_R2,
        "--output-R1", tmp_dup_R1,
        "--output-R2", tmp_dup_R2
    ]
    cmd = [str(elt) for elt in cmd]
    log.info("Add duplicates: {}".format(" ".join(cmd)))
    subprocess.check_call(cmd)
    for curr_file in [tmp_sim_R1, tmp_sim_R2]:
        os.remove(curr_file)

    # Apply quality profile
    cmd = [
        "applyQualityProfile.py",
        "--random-seed", args.random_seed,
        "--qual-offset", args.qual_offset,
        "--qual-penalty", args.qual_penalty,
        "--reads-length", args.reads_length,
        "--input-sequences", tmp_dup_R1,
        "--output-sequences", args.output_R1,
        "--output-sequences", args.output_R1,
        "--output-trace", args.output_errors_R1,
        "--input-models"
    ]
    cmd.extend(args.input_R1_quality_models)
    cmd = [str(elt) for elt in cmd]
    log.info("Apply quality profile for R1: {}".format(" ".join(cmd)))
    subprocess.check_call(cmd)
    os.remove(tmp_dup_R1)
    cmd = [
        "applyQualityProfile.py",
        "--random-seed", args.random_seed,
        "--qual-offset", args.qual_offset,
        "--qual-penalty", args.qual_penalty,
        "--reads-length", args.reads_length,
        "--input-sequences", tmp_dup_R2,
        "--output-sequences", args.output_R2,
        "--output-sequences", args.output_R2,
        "--output-trace", args.output_errors_R2,
        "--input-models"
    ]
    cmd.extend(args.input_R2_quality_models)
    cmd = [str(elt) for elt in cmd]
    log.info("Apply quality profile for R2: {}".format(" ".join(cmd)))
    subprocess.check_call(cmd)
    os.remove(tmp_dup_R2)

    log.info("End of job")
