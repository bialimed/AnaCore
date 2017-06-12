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
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFIO



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Sorts VCF by coordinates.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.info = FH_in.info
            FH_out.format = FH_in.format
            FH_out.samples = FH_in.samples
            FH_out._writeHeader()
        
            # Records
            records_by_chr = dict()
            for record in FH_in:
                if record.chrom not in records_by_chr:
                    records_by_chr[record.chrom] = list()
                records_by_chr[record.chrom].append(record)
            for chrom in sorted(records_by_chr):
                sorted_records = sorted(records_by_chr[chrom], key=lambda x: (x.chrom, x.pos))
                for record in sorted_records:
                    FH_out.write( record )
