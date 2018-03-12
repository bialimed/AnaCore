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

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from bed import getAreasByChr
from sequenceIO import Sequence, FastaIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def revcom( seq ):
    """
    @summary: Returns the reverse complement the sequence.
    @param seq: [str] The sequence.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}

    return( "".join([complement_rules[base] for base in seq[::-1]]) )



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Extract primers sequence from Illumina's amplicons manifest." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-s', '--input-sequences', required=True, help='*** (format: Fasta).' )
    group_input.add_argument( '-r', '--input-regions', required=True, help='*** (format: BED).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-fasta', default="regions.fasta", help='******************************** (format: fasta). [Default: %(default)s]' )
    args = parser.parse_args()

    # Process
    ampl_by_chr = getAreasByChr( args.input_regions )
    with FastaIO(args.input_sequences) as FH_seq:
        for record in FH_seq:
            chrom = record.id
            if chrom in ampl_by_chr:
                for area in ampl_by_chr[chrom]:
                    area.annot = { "sequence": record.string[area.start-1:area.end] } ################ stranded if ask
    with FastaIO(args.output_fasta, "w") as FH_out:
        for chrom in sorted(ampl_by_chr):
            for area in ampl_by_chr[chrom]:
                FH_out.write(
                    Sequence(
                        area.name,
                        area.annot["sequence"],
                        "pos:" + area.chrom + ":" + str(area.start) + "-" + str(area.end) + " strand:" + area.strand
                    )
                )
