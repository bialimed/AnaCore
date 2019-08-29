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
__version__ = '2.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, VCFRecord
from anacore.sequenceIO import IdxFastaIO
from anacore.bed import getSortedAreasByChr


########################################################################
#
# FUNCTIONS
#
########################################################################
def isOverlapping(regions_by_chr, variant, seq_handler):
    """
    Return True if the variant overlap one region in regions_by_chr.

    :param regions_by_chr: By chromosome the list of selected regions. Each list of region is sorted firstly by start position (1-based) and secondly by end position (1-based).
    :type regions_by_chr: dict
    :param variant: The evaluated variant.
    :type variant: anacore.vcf.VCFRecord
    :param seq_handler: The file handler on reference sequences file.
    :type seq_handler: anacore.sequenceIO.IdxFastaIO
    :return: True if the variant overlap one region in regions_by_chr.
    :rtype: bool
    """
    is_overlapping = False
    if variant.chrom in regions_by_chr:
        # Standard coordinates
        variant_start = variant.refStart()
        variant_end = variant.refEnd()
        for region in regions_by_chr[variant.chrom]:
            if variant_end >= region.start and variant_start <= region.end:
                is_overlapping = True
        if not is_overlapping and variant.isIndel():
            chrom_seq = seq_handler.get(variant.chrom).string
            # Upstream coordinates
            up_variant = variant.getMostUpstream(chrom_seq)
            if up_variant.pos != variant.pos or up_variant.ref != variant.ref.replace(VCFRecord.getEmptyAlleleMarker(), ""):
                variant_start = up_variant.refStart()
                variant_end = up_variant.refEnd()
                for region in regions_by_chr[variant.chrom]:
                    if variant_end >= region.start and variant_start <= region.end:
                        is_overlapping = True
            # Downstream coordinates
            if not is_overlapping:
                down_variant = variant.getMostDownstream(chrom_seq)
                if down_variant.pos != variant.pos or down_variant.ref != variant.ref.replace(VCFRecord.getEmptyAlleleMarker(), ""):
                    variant_start = down_variant.refStart()
                    variant_end = down_variant.refEnd()
                    for region in regions_by_chr[variant.chrom]:
                        if variant_end >= region.start and variant_start <= region.end:
                            is_overlapping = True
    return is_overlapping


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filters variants by location. Each variant not located on one of the selected regions is removed.')
    parser.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" the filter tag OOT is added to the variants out of target. In mode "remove" the variant out of target are removed. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_input.add_argument('-t', '--input-targets', required=True, help='Path to the file describing selected regions (format: BED).')
    group_input.add_argument('-r', '--input-reference', required=True, help='Path to the sequences file used in alignment (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', default="filtered.vcf", help='Path to the outputted variants file (format: VCF). [Default: %(default)s]')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_variants = 0
    nb_kept = 0
    kept_by_chr = getSortedAreasByChr(args.input_targets)
    with IdxFastaIO(args.input_reference, use_cache=True) as FH_seq:
        with VCFIO(args.input_variants) as FH_in:
            with VCFIO(args.output_variants, "w") as FH_out:
                # Writes header
                FH_out.copyHeader(FH_in)
                FH_out.filter["OOT"] = 'The variant is out of targeted regions (' + args.input_targets + ').'
                FH_out.writeHeader()
                # Writes variants
                for variant in FH_in:
                    nb_variants += 1
                    overlap_targets = False
                    if isOverlapping(kept_by_chr, variant, FH_seq):
                        nb_kept += 1
                        overlap_targets = True
                    if args.mode == "remove":
                        if overlap_targets:
                            FH_out.write(variant)
                    else:
                        if variant.filter is None or len(variant.filter) == 0:  # Init filter tags
                            variant.filter = ["PASS"]
                        if not overlap_targets:  # Add OOT filter tag
                            if len(variant.filter) == 1 and variant.filter[0] == "PASS":
                                variant.filter = ["OOT"]
                            else:
                                variant.filter.append("OOT")
                        FH_out.write(variant)
    # Log process
    log.info(
        "{:.2%} of variants have been {} ({}/{})".format(
            0 if nb_variants == 0 else (nb_variants - nb_kept) / nb_variants,
            "tagged" if args.mode == "tag" else "removed",
            nb_variants - nb_kept,
            nb_variants
        )
    )
    log.info("End of job")
