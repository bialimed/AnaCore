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
__email__ = 'escudie.frederiic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import re
import sys
import argparse



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merge HTSeq-count outputs for each samples in one.' )
    parser.add_argument( '-p', '--name-pattern', default='^([^\.]+)\.', help='The regexp used to extract sample name from HTSeq-count outputs filenames.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-files', nargs='+', required=True, help='The HTSeq-count outputs to merge (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='The path for the outputed file (format: TSV).')
    args = parser.parse_args()

    sample_names = list()

    # Get count by sample by reference
    spl_count_by_ref = dict()
    for sample_file in args.input_files:
        match = re.search( args.name_pattern, os.path.basename(sample_file) )
        if match is None:
            raise Exception( "The sample name cannot be extracted from the file '" + sample_file + "' with the regexp '" + args.name_pattern + "'" )
        spl = match.group(1)
        sample_names.append( spl )
        with open(sample_file) as FH_spl:
            for line in FH_spl:
                ref, count = [field.strip() for field in line.split("\t")]
                if ref not in ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]:
                    if ref not in spl_count_by_ref:
                        spl_count_by_ref[ref] = dict()
                    spl_count_by_ref[ref][spl] = int(count)

    # Write output
    with open(args.output_file, "w") as FH_out:
        FH_out.write( "#Reference"  + "\t" + "\t".join(sample_names) + "\n" )
        for ref in spl_count_by_ref:
            counts = list()
            for spl in sample_names:
                if spl not in spl_count_by_ref[ref]:
                    counts.append(0)
                else:
                    counts.append(spl_count_by_ref[ref][spl])
            FH_out.write( ref + "\t" + "\t".join(map(str, counts)) + "\n" )
