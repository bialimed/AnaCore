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
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastqIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Change quality offset in fastq file.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-n', '--new-offset', default=33, type=int, help="New quality offset. [Default: %(default)s]")
    parser.add_argument('-b', '--old-offset', type=int, help="The current quality offset. [Default: Auto determined]")
    parser.add_argument('-mi', '--min-qual', type=int, help="The minimun quality. [Default: No minimun]")
    parser.add_argument('-ma', '--max-qual', type=int, help="The maximun quality. [Default: No maximun]")
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-file', required=True, help='Path to the sequences file (format: fastq).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-file', required=True, help='Path to the output (format: fastq).')
    args = parser.parse_args()

    # Process
    old_offset = args.old_offset if args.old_offset is not None else FastqIO.qualOffset(args.input_file)
    if old_offset is None:
        raise Exception("The quality offset in {} cannot be determined.".format(args.input_file))
    offset_modifier = args.new_offset - old_offset
    with FastqIO(args.output_file, "w") as FH_out:
        with FastqIO(args.input_file) as FH_in:
            for record in FH_in:
                new_qual = ""
                for curr_qual in record.quality:
                    new_qual_numer = ord(curr_qual) + offset_modifier
                    if args.min_qual is not None:
                        new_qual_numer = max(args.new_offset + args.min_qual, new_qual_numer)
                    if args.max_qual is not None:
                        new_qual_numer = min(args.new_offset + args.max_qual, new_qual_numer)
                    new_qual += chr(new_qual_numer)
                record.quality = new_qual
                FH_out.write(record)
