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
import pysam
import argparse



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='*****************************.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-a', '--allele', required=True, help='***************.' )
    group_area = parser.add_argument_group( 'Area' ) # Area
    group_area.add_argument( '-r', '--region', required=True, help='***************.' )
    group_area.add_argument( '-s', '--start', type=int, required=True, help='*************** (1-based).' )
    group_area.add_argument( '-e', '--end', type=int, required=True, help='*************** (1-based).' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-f', '--sam-files', nargs='+', required=True, help='*************** (format: SAM/BAM/CRAM).' )
    #~ group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    #~ group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).' )
    args = parser.parse_args()

    # Process
    samples = list()
    alleles = dict()
    for curr_sam in args.sam_files:
        curr_spl = os.path.basename(curr_sam)
        samples.append( curr_spl )
        alleles[curr_spl] = dict()
        sam_FH = pysam.AlignmentFile(curr_sam, "rb")
        for pileupcolumn in sam_FH.pileup(args.region, args.start - 1, args.end - 1):
            if pileupcolumn.pos >= (args.start - 1) and pileupcolumn.pos <= (args.end - 1):
                allele = {
                    'A': 0,
                    'T': 0,
                    'G': 0,
                    'C': 0,
                    'N': 0,
                    'del': 0
                }
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del:
                        allele['del'] += 1
                    elif not pileupread.is_refskip: # if the base on the padded read is not a part of CIGAR N op.
                        nt = pileupread.alignment.query_sequence[pileupread.query_position]
                        allele[nt.upper()] += 1
                alleles[curr_spl][pileupcolumn.pos + 1] = allele
        sam_FH.close()

    # Display
    print( "Pos", "\t".join(samples), sep="\t" )
    for curr_pos in range(args.start, args.end + 1):
        pos_data = [args.region + ":" + str(curr_pos)]
        for curr_spl in samples:
            evaluated_allele = 0
            others_alleles = 0
            if curr_pos in alleles[curr_spl]:
                for curr_allele in alleles[curr_spl][curr_pos]:
                    if curr_allele == args.allele:
                        evaluated_allele += alleles[curr_spl][curr_pos][curr_allele]
                    else:
                        others_alleles += alleles[curr_spl][curr_pos][curr_allele]
            pos_data.append( str(evaluated_allele) + "/" + str(evaluated_allele + others_alleles) )
        print( "\t".join(pos_data) )
