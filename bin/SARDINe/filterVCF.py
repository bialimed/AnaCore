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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import *


########################################################################
#
# FUNCTIONS
#
########################################################################
def getExclusionDict( excluded_regions ):
    excluded_by_chr = dict()
    # Group by chromosome
    for region in excluded_regions:
        chrom, start_end = region.split(":")
        start, end = start_end.split("-")
        if chrom not in excluded_by_chr:
            excluded_by_chr[chrom] = list()
        excluded_by_chr[chrom].append( {"start":int(start), "end":int(end)} )
    # In each chromosome sort by position
    for chrom in excluded_by_chr:
        excluded_by_chr[chrom] = sorted(excluded_by_chr[chrom], key=lambda x: (x["start"], x["end"]))
    return excluded_by_chr

def isOverlapping( regions_by_chr, chrom, pos ):
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
    parser = argparse.ArgumentParser( description='****************************************************************.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_filter = parser.add_argument_group( 'Filters' ) # Filters
    group_filter.add_argument( '-r', '--excluded-regions', nargs='+', default=[], help='*********************************************' )
    group_filter.add_argument( '-t', '--excluded-types', nargs='+', default=[], choices=["snp", "indel", "variation"], help='****************************************************' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='********************************************************** (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', default="filtered.vcf", help='************************************************************* (format: VCF). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    excluded_by_chr = getExclusionDict( args.excluded_regions )
    with open(args.output_variants, "w") as FH_out:
        # Writes header
        with open(args.input_variants) as FH_in:
            line = FH_in.readline()
            while line.startswith("#"):
                FH_out.write( line )
                line = FH_in.readline()
		# Writes variants
        with VCFIO(args.input_variants) as FH_in:
            for variant in FH_in:
                if variant.type() not in args.excluded_types:
                    if not isOverlapping(excluded_by_chr, variant.chrom, variant.pos):
                        FH_out.write( FH_in.recToVCFLine(variant) + "\n" )
