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
__status__ = 'prod'

import os
import sys
import uuid
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

def mostUpstream( genome_path, in_variant_file, out_variant_file ):
    genome_by_chr = getSeqByChr( genome_path )
    with open(out_variant_file, "w") as FH_out:
        # Header
        ############################################################### Add tag BK=10:25010|./ATG
        with open(in_variant_file) as FH_vcf:
            line = FH_vcf.readline()
            while line.startswith("#"):
                FH_out.write( line )
                line = FH_vcf.readline()
        # Variants
        ############################################################### Add tag BK=10:25010|./ATG
        with VCFIO(in_variant_file) as FH_in:
            for record in FH_in:
                if record.isIndel(): # Variant is indel
                    ref = record.ref
                    alt = record.alt[0]
                    # Deletion
                    if alt == ".":
                        if genome_by_chr[record.chrom][(record.pos - 1):(record.pos - 1 + len(ref))] != ref:
                            raise Exception('The reference on position ' + record.chrom + ':' + str(record.pos) + ' does not correspond to "' + ref + '".')
                        before_var = genome_by_chr[record.chrom][0:(record.pos - 1)]
                        while before_var[-1] == ref[-1]:
                            # shift to upstream
                            before_var = before_var[:-1]
                            ref = ref[-1] + ref[:-1]
                        record.pos = len(before_var) + 1
                        record.ref = ref
                    # Insertion
                    else:
                        before_var = genome_by_chr[record.chrom][0:(record.pos - 1)]
                        while before_var[-1] == alt[-1]:
                            # shift to upstream
                            before_var = before_var[:-1]
                            alt = alt[-1] + alt[:-1]
                        record.pos = len(before_var) + 1
                        record.alt = [alt]
                FH_out.write( FH_in.recToVCFLine(record) + "\n" )

def splitAndStd( in_variant_file, out_variant_file ):
    with open(out_variant_file, "w") as FH_out:
        # Header
        with open(in_variant_file) as FH_vcf:
            line = FH_vcf.readline()
            while line.startswith("#"):
                FH_out.write( line )
                line = FH_vcf.readline()
        # Variants
        with VCFIO(in_variant_file) as FH_in:
            for record in FH_in:
                # For each variant alternative allele
                for idx_alt, alt in enumerate(record.alt):
                    new_record = getAlleleRecord( FH_in, record, idx_alt )
                    new_record.standardizeSingleAllele( "." )
                    FH_out.write( FH_in.recToVCFLine(new_record) + "\n" )



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Splits alternatives alleles of variants in multi lines and standardizes indel representation.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the variant file (format: VCF).' )
    group_input.add_argument( '-g', '--input-genome', help='Genome reference used in variant calling to produced the inputed VCF (format: fasta). With this option all the indel are moved at their maximum upstream. This is used only in comparison for amplicon analysis. Example with repeat: TGATCATGATGTC with variant ATG/. on position 9 becomes ATG/. on position 6. Example with shift: ATTCCTCAAATAAGATAT with variant AATAA/. on position 8 becomes AAATA/. on position 7.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()
    
    # Process
    out_dir = os.path.dirname( args.output_variants )
    tmp_vcf = os.path.join( out_dir, str(uuid.uuid1()) )
    if args.input_genome is None:
        tmp_vcf = args.output_variants
    splitAndStd( args.input_variants, tmp_vcf )
    if args.input_genome is not None:
        mostUpstream( args.input_genome, tmp_vcf, args.output_variants )
        os.remove( tmp_vcf )
