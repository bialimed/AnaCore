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

import argparse


########################################################################
#
# FUNCTIONS
#
########################################################################
def getNonOverlappingThread(selected_areas):
    """
    @summary: Splits the list of regions in sub-lists of non-overlapping regions and returns them.
    @param selected_areas: [list] The selected areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    @return: [list] Each element of the list is a list of non-overlapping regions.
    """
    sorted_selected_areas = sorted(selected_areas, key=lambda x: (x["region"], x["start"], x["end"]))
    non_overlapping_threads = [[
        sorted_selected_areas[0]
    ]]
    for curr_area in sorted_selected_areas[1:]:
        thread_idx = 0
        is_valid_thread = False
        while thread_idx < len(non_overlapping_threads) and not is_valid_thread:
            is_valid_thread = hasNoOverlap(curr_area, non_overlapping_threads[thread_idx][-1])
            thread_idx += 1
        if is_valid_thread:
            non_overlapping_threads[thread_idx - 1].append(
                curr_area
            )
        else:
            non_overlapping_threads.append([
                curr_area
            ])
    return non_overlapping_threads

def getSelectedArea(input_panel, margin):
    """
    @summary: Returns the list of selected areas from a BED file.
    @param input_panel: [str] The path to the selected areas description (format: BED).
    @param margin: [int] The minimum distance between two areas in same group. With 0 the group contains non-overlapping areas. With 5 the areas in the same group are non-overlapping but also separated by at least 5 nucleotids. This option is used when the sequencing adapter cannot be totally removed to prevent overlap between primer of one area and the adapters of an other.
    @return: [list] The list of BED's areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    """
    selected_areas = list()
    with open(input_panel) as FH_panel:
        for line in FH_panel:
            if not line.startswith("browser ") and not line.startswith("track ") and not line.startswith("#"):
                fields = [elt.strip() for elt in line.split("\t")]
                selected_areas.append({
                    "region": fields[0],
                    "start": max(1, (int(fields[1]) + 1 - margin)),  # Start in BED is 0-based
                    "end": int(fields[2]) + margin,
                    "id": fields[3]
                })
    return selected_areas

def hasNoOverlap(area_A, area_B):
    """
    @summary: Returns True if the area_A does not overlap area_B.
    @param area_A: The first evaluated area.
    @param area_B: The second evaluated area.
    @return: [bool] True if the area does not overlap the last area in area_thread.
    """
    has_no_overlap = False
    if area_A["region"] != area_B["region"]:  # New chromosome
        has_no_overlap = True
    elif area_A["start"] > area_B["end"]:  # Same region but without overlap
        has_no_overlap = True
    return has_no_overlap


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Writes non-overlapping areas groups from a BED file.')
    parser.add_argument('-m', '--margin', default=0, type=int, help='The minimum distance between two areas in same group. With 0 the group contains non-overlapping areas. With 5 the areas in the same group are non-overlapping but also separated by at least 5 nucleotids. This option is used when the sequencing adapter cannot be totally removed to prevent overlap between primer of one area and the adapters of an other. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-panel', required=True, help='Path to the list of selected areas (format: BED). Each area must have an unique ID in the name field.')
    group_output = parser.add_argument_group('Inputs')  # Inputs
    group_output.add_argument('-o', '--output-design', required=True, help='Path to the list of amplicons (format: TSV). The first column is the ID of the amplicon and the second is the name of the group where the amplicon has no overlap with other amplicons of this group.')
    args = parser.parse_args()

    # Get selected area from BED
    selected_areas = getSelectedArea(args.input_panel, args.margin)

    # Split overlapping area
    non_overlapping_threads = getNonOverlappingThread(selected_areas)

    # Write split design
    with open(args.output_design, "w") as FH_design:
        FH_design.write("#Area\tNon-overlapping_group\n")
        for thread_idx, thread in enumerate(non_overlapping_threads):
            for area in non_overlapping_threads[thread_idx]:
                FH_design.write(
                    '{0}\t{1}\n'.format(area["id"], "grp" + str(thread_idx))
                )
