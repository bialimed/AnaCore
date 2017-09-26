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
import json
import argparse


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Parse areas depths distributions to list areas with an depths under the threshold on a sufficient length.' )
    parser.add_argument( '-p', '--percentile', type=int, default=15, help='The evaluated percentile. [Default: %(default)s]' )
    parser.add_argument( '-d', '--min-depth', type=int, default=200, help='Under this value of depth on the selected percentile the area is reported. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-depths', required=True, nargs='+', help='For each evaluated sample a file containing the depths by percentile in selected areas (format: JSON outputted by coverageArea.py).')
    args = parser.parse_args()

    # Process
    print("## percentile:", args.percentile, "; min_depth:", args.min_depth)
    print("#Sample", "Area", "Depth", sep="\t")
    for curr_file in args.input_depths:
        spl_name = os.path.basename(curr_file).split("_depths.json")[0]
        with open(curr_file, 'r') as FH:
            data = json.load( FH )
            for area in data:
                for name in area["depths"]:
                    curr_depth = int(area["depths"][name][str(args.percentile) + "_percentile"])
                    if curr_depth < args.min_depth:
                        print(spl_name, area["name"], curr_depth, sep="\t")
