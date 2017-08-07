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
import pysam
import argparse



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Splits one BAM in groups based on RG. Several RG can be merge in new group. Each new group becomes represented by one alignment file after split.' )
    parser.add_argument( '-r', '--remove-RG', action='store_true', help='With this parameter the RG are removed from the outputted alignments files.' )
    parser.add_argument( '-t', '--RG-tag', default='LB', help='RG tag used in link between tag value and group (see input-design parameter). [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-a', '--input-aln', required=True, help='The path to the alignments files (format: BAM).' )
    group_input.add_argument( '-d', '--input-design', required=True, help='The path to the file describing RG in each new group (format: TSV). First column is a value for a specific tag in RG, the second is the name of the new group.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-p', '--output-pattern', default="out_{GP}.bam", help='The path pattern for the outputted alignments files (format: BAM). In this path the keyword "{GP}" is replace by the group name for each group. [Default: %(default)s]' )
    args = parser.parse_args()

    # Get panel regions
    groups_names = set()
    group_by_tag = dict()
    with open(args.input_design) as FH_design:
        for line in FH_design:
            if not line.startswith( "#" ):
                read_tag, group = [elt.strip() for elt in line.split("\t")]
                group_by_tag[read_tag] = group
                groups_names.add( group )

    # Split BAM
    with pysam.AlignmentFile( args.input_aln, "rb" ) as FH_in:
        # Get new group by read group ID
        group_by_id = dict()
        for RG in FH_in.header["RG"]:
            read_tag = RG[args.RG_tag]
            group_by_id[RG["ID"]] = group_by_tag[read_tag]
        # Open FH on outputs
        FH_by_group = dict()
        for group in groups_names:
            new_header = FH_in.header.copy()
            if args.remove_RG:
                new_header["RG"] = list()
            else:
                new_header["RG"] = [curr_RG for curr_RG in FH_in.header["RG"] if group_by_tag[curr_RG[args.RG_tag]] == group]
            FH_by_group[group] = pysam.AlignmentFile(
                args.output_pattern.replace("{GP}", group),
                "wb",
                header=new_header
            )
        # Parse reads
        for curr_read in FH_in.fetch(until_eof=True):
            if curr_read.has_tag("RG"):
                RG_id = curr_read.get_tag("RG")
                if args.remove_RG:
                    curr_read.set_tag( "RG", None )
                FH_by_group[group_by_id[RG_id]].write( curr_read )
        # Close FH
        for group in groups_names:
            FH_by_group[group].close()
