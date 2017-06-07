#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT
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
__copyright__ = 'Copyright (C) 2017 IUCT'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from VEPvcf import VEPVCFIO, getAlleleRecord



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Converts VCF annotated with VEP in separated value format (CSV, TSV, ...). One line in output file represents an annotation.' )
    parser.add_argument( '-s', '--separator', default='\t', help="Field separator in output file. [Default: tab]")
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the file file containing variants and annotated with VEP v88+ (format: VCF).' )
    args = parser.parse_args()

    # Process
    with VEPVCFIO(args.input_variants) as FH_vcf:
        print(
            "Chromosome",
            "Position",
            "Reference_allele",
            "Alternative_alleles",
            "Variant_quality",
            args.separator.join( FH_vcf.CSQ_titles ),
            sep=args.separator
        )
        for record in FH_vcf:
            for idx_csq, csq in enumerate(record.info["CSQ"]):
                csq_values = list()
                for title in FH_vcf.CSQ_titles:
                    if csq[title] is None:
                        csq_values.append("")
                    else:
                        csq_values.append( str(csq[title]) )
                print(
                    (record.chrom if idx_csq == 0 else ""),
                    (record.pos if idx_csq == 0 else ""),
                    (record.ref if idx_csq == 0 else ""),
                    (",".join(record.alt) if idx_csq == 0 else ""),
                    (record.qual if idx_csq == 0 else ""),
                    args.separator.join( csq_values ),
                    sep=args.separator
                )
