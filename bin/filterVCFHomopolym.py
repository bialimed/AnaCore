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
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, getAlleleRecord, HeaderFilterAttr
from anacore.sequenceIO import IdxFastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def isOnHomopolymer(FH_fasta_idx, record, homopolym_length):
    """
    Return True is the variant is adjacent to an homopolymer.

    :param FH_fasta_idx: File handle to the reference file.
    :type FH_fasta_idx: anacore.sequenceIO.IdxFastaIO
    :param record: The variant.
    :type record: anacore.vcf.VCFRecord
    :param homopolym_length: The variant is flagged as adjacent to an homopolymer if the previous or next nucleotid is repeated at least this number of times.
    :type homopolym_length: int
    """
    is_on_homopolym = False
    # Check the downstream
    after_start = int(record.refEnd()) + 1
    after_end = after_start + homopolym_length - 1
    after = FH_fasta_idx.getSub(record.chrom, after_start, after_end).upper()
    after_nt = set(after)
    if len(after_nt) == 1:
        is_on_homopolym = True
    # Check the upstream
    if not is_on_homopolym:
        before_end = int(record.refStart() + 0.5) - 1
        before_start = before_end - homopolym_length + 1
        before = FH_fasta_idx.getSub(record.chrom, before_start, before_end).upper()
        before_nt = set(before)
        if len(before_nt) == 1:
            is_on_homopolym = True
    return is_on_homopolym


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter the variants adjacents of homopolymers.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-l', '--homopolym-length', type=int, default=4, help='The variant is flagged as adjacent to an homopolymer if the previous or next nucleotid is repeated at least this number of times. [Default: %(default)s]')
    group_filter.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant is adjacent to an homopolymer a tag is added in FILTER (cf. "--tag-name"). In mode "remove" if the variant is adjaacent to an homopolymer it is removed. [Default: %(default)s]')
    group_filter.add_argument('-t', '--tag-name', default="homoP", help='The name of the tag added on variant. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_input.add_argument('-r', '--input-reference', required=True, help='The path to the reference file (format: fasta). It must be indexed with faidx.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the filtered file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_variants = 0
    nb_filtered = 0
    with IdxFastaIO(args.input_reference) as FH_ref:
        with VCFIO(args.input_variants, "r") as FH_in:
            with VCFIO(args.output_variants, "w") as FH_out:
                # Header
                FH_out.copyHeader(FH_in)
                FH_out.filter[args.tag_name] = HeaderFilterAttr(args.tag_name, "The variant is adjacent to an homopolymer (repeat size >= {}).".format(args.homopolym_length))
                FH_out.writeHeader()
                # Records
                for alt_record in FH_in:
                    nb_variants += 1
                    is_on_hpolym = isOnHomopolymer(FH_ref, alt_record, args.homopolym_length)
                    if is_on_hpolym:
                        nb_filtered += 1
                    if args.mode == "tag":  # Tag mode
                        # Init filters
                        if alt_record.filter is None or len(alt_record.filter) == 0 or alt_record.filter[0] == "PASS":
                            alt_record.filter = list()
                        # Add tag
                        if is_on_hpolym:
                            alt_record.filter.append(args.tag_name)
                        elif len(alt_record.filter) == 0:
                            alt_record.filter.append("PASS")
                        # Write record
                        FH_out.write(alt_record)
                    elif not is_on_hpolym:  # Filter mode and is not on homopolymer
                        FH_out.write(alt_record)

    # Log process
    log.info(
        "{:.2%} of variants have been {} ({}/{})".format(
            0 if nb_variants == 0 else nb_filtered / nb_variants,
            "tagged" if args.mode == "tag" else "removed",
            nb_filtered,
            nb_variants
        )
    )
    log.info("End of job")
