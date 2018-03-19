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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getKeptFromBed(bed_path):
    """
    @summary: Returns selected regions by chromosome from BED.
    @param bed_path: [str] Path to the file describing selected regions (format: BED).
    @return: [dict] By chromosome the list of selected regions. Each list of region is sorted firstly by start position (1-based) and secondly by end position (1-based). Each region is represented by a dictionary with this format: {"start":501, "end":608}.
    """
    kept_by_chr = dict()
    # Group by chromosome
    with open(bed_path) as FH_in:
        for line in FH_in:
            row_data = [elt.strip() for elt in line.split("\t")]
            chrom = row_data[0]
            start = min(int(row_data[1]) + 1, int(row_data[2]))
            end = max(int(row_data[1]) + 1, int(row_data[2]))
            if chrom not in kept_by_chr:
                kept_by_chr[chrom] = list()
            kept_by_chr[chrom].append({"start": start, "end": end})
    # In each chromosome sort by position
    for chrom in kept_by_chr:
        kept_by_chr[chrom] = sorted(kept_by_chr[chrom], key=lambda x: (x["start"], x["end"]))
    return kept_by_chr

def isOverlapping(regions_by_chr, chrom, pos):
    """
    @summary: Returns True if the specified position overlap one region in regions_by_chr.
    @param regions_by_chr: [dict] By chromosome the list of selected regions. Each list of region is sorted firstly by start position (1-based) and secondly by end position (1-based). Each region is represented by a dictionary with this format: {"start":501, "end":608}.
    @param chrom: [str] The chromosome where is located the evaluated position.
    @param pos: [int] The evaluated position on chrom.
    @return: [bool] True if the specified position overlap one region in regions_by_chr.
    """
    is_overlapping = False
    if chrom in regions_by_chr:
        for region in regions_by_chr[chrom]:
            if region["start"] <= pos and region["end"] >= pos:
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
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_input.add_argument('-b', '--selected-regions', required=True, help='Path to the file describing selected regions (format: BED).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', default="filtered.vcf", help='Path to the outputted variants file (format: VCF). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    kept_by_chr = getKeptFromBed(args.selected_regions)
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Writes header
            FH_out.copyHeader(FH_in)
            FH_out.filter["REG"] = 'The variant is located on an excluded region (' + args.selected_regions + ').'
            FH_out._writeHeader()
            # Writes variants
            for variant in FH_in:
                if isOverlapping(kept_by_chr, variant.chrom, variant.pos):
                    FH_out.write(variant)
