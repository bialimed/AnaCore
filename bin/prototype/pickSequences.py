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

from sequenceIO import SequenceFileReader, FastaIO, FastqIO


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
    group_input.add_argument( '-i', '--input-ids', required=True, help='*** (format: TXT).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output', help='******************************** (format: fasta). [Default: %(default)s]' )
    args = parser.parse_args()

    # Process
    selected_ids = dict()
    with open(args.input_ids) as FH_ids:
        for line in FH_ids:
            selected_ids[line.strip()] = 0
    FH_seq = SequenceFileReader.factory( args.input_sequences )
    try:
        FH_out = None
        if issubclass(FastqIO, FH_seq.__class__):
            FH_out = FastqIO(args.output, "w")
        else:
            FH_out = FastaIO(args.output, "w")
        try:
            for record in FH_seq:
                if record.id in selected_ids:
                    selected_ids[record.id] = 1
                    FH_out.write( record )
        finally:
            FH_out.close()
    finally:
        FH_seq.close()
    for curr_id in selected_ids:
        if selected_ids[curr_id] != 1:
            print("Pb", curr_id)
