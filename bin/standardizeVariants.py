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
__version__ = '1.2.0'
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

from vcf import VCFIO, getAlleleRecord
from sequenceIO import FastaIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def getSeqByChr( genome_path ):
    genome_by_chr = dict()
    FH_seq = FastaIO( genome_path )
    for record in FH_seq:
        genome_by_chr[record.id] = record.string.upper()
    FH_seq.close()
    return genome_by_chr

def stdAndMove( genome_path, in_variant_file, out_variant_file ):
    genome_by_chr = getSeqByChr( genome_path )
    with VCFIO(out_variant_file, "w") as FH_out:
        with VCFIO(in_variant_file) as FH_in:
            # Header
            FH_out.copyHeader( FH_in )
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                curr_chrom = genome_by_chr[record.chrom]
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord( FH_in, record, alt_idx )
                    FH_out.write( alt_record.getMostUpstream(curr_chrom) )

def stdOnly( in_variant_file, out_variant_file ):
    with VCFIO(out_variant_file, "w") as FH_out:
        with VCFIO(in_variant_file) as FH_in:
            # Header
            FH_out.copyHeader( FH_in )
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord( FH_in, record, alt_idx )
                    alt_record.standardizeSingleAllele()
                    FH_out.write( alt_record )



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Splits alternatives alleles of one variants in multi-lines and removes unecessary reference and alternative nucleotids and move indel to most upstream position.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the variant file (format: VCF).' )
    group_input.add_argument( '-g', '--input-genome', help='Genome reference used in variant calling to produced the inputed VCF (format: fasta). With this option all the indel are moved at their maximum upstream. This is used only in comparison for amplicon analysis. Example with repeat: TGATCATGATGTC with variant ATG/. on position 9 becomes ATG/. on position 6. Example with shift: ATTCCTCAAATAAGATAT with variant AATAA/. on position 8 becomes AAATA/. on position 7.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()

    # Process
    if args.input_genome is None:
        stdOnly( args.input_variants, args.output_variants )
    else:
        stdAndMove( args.input_genome, args.input_variants, args.output_variants )
