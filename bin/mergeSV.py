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

import argparse



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Merges several separated values files by a linking column. Only the entry in first file are kept and annotated by entry in others files.' )
    parser.add_argument( '-t', '--new-title', help='The new title for link column. [Default: link title of the first file merged]' )
    parser.add_argument( '-l', '--links-titles', nargs='+', default=[], help='For each file the columns used as link between files (in same order as --inputs).' )
    parser.add_argument( '-is', '--in-separator', default="\t", help='The field separator in input files. [Default: %(default)s]' )
    parser.add_argument( '-os', '--out-separator', default="\t", help='The field separator in output file. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-files', nargs='+', required=True, help='The path of the files to merge (format: separated values files like CSV or TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', default='merged.tsv', help='Path of the output file (format: separated values files like CSV or TSV). [Default: %(default)s]')
    args = parser.parse_args()

    new_title = args.new_title
    link_titles = args.links_titles
    out_titles = list()
    data_by_link = dict()

    # Parses files
    processed_items = dict()
    for idx_file, curr_file in enumerate(args.input_files):
        with open(curr_file) as FH_in:
            titles = [field.strip() for field in FH_in.readline().split(args.in_separator)]
            # Store titles
            link_idx = 0
            for idx_title, curr_title in enumerate(titles):
                if curr_title != link_titles[idx_file]: # The current column is not the linking column
                    out_titles.append( curr_title.replace(args.out_separator, hex(ord(args.out_separator))) ) # Prevent bug in out separator
                else:
                    link_idx = idx_title
                    if new_title is None: # If no new title has been define for the linking column and this is the first file
                        new_title = cur_title
            # Store data
            for link_id in processed_items: # Reset items processed in current file
                processed_items[link_id] = False
            for line in FH_in:
                fields = [field.strip() for field in line.split(args.in_separator)]
                link_id = fields[link_idx]
                if idx_file == 0:
                    data_by_link[link_id] = list()
                    processed_items[link_id] = True
                if link_id in data_by_link: # Kept only entry present in first file
                    processed_items[link_id] = True
                    for idx_field, cur_field in enumerate(fields):
                        if titles[idx_field] != link_titles[idx_file]:
                            cur_field = cur_field.replace( args.out_separator, hex(ord(args.out_separator)) ) # Prevent bug in out separator
                            data_by_link[link_id].append( cur_field )
            for link_id, is_processed in processed_items.items(): # Complete items not in current file
                if not is_processed:
                    data_by_link[link_id].extend(["" for col in range(len(titles) - 1)])

    # Writes output
    with open(args.output_file, "w") as FH_out:
        # Titles
        FH_out.write( 
            "{}{}{}\n".format(
                new_title, args.out_separator, args.out_separator.join(out_titles)
            )
        )
        # Data
        for link_val in data_by_link:
            FH_out.write(
                "{}{}{}\n".format(
                    link_val, args.out_separator, args.out_separator.join(data_by_link[link_val])
                )
            )
