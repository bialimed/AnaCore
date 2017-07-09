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
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import json
import numpy
import argparse


########################################################################
#
# FUNCTIONS
#
########################################################################
def getDistributionDict( values, percentile_step=25, precision=4 ):
    """
    @summary: Returns the distribution of values (min, max and percentiles).
    @param values: [list] The values.
    @param percentile_step: [int] Only this percentile and his multiples are returned.
    @param precision: [int] The decimal precision.
    @retrun: [dict] The min, max and percentiles values. Example: {"min":0, "05_percentile":10, "10_percentile":15, ..., "95_percentile":853, "max":859}
    """
    distrib = {
        "min": round(min(values), precision),
        "max": round(max(values), precision)
    }
    for curr_percentile in range(percentile_step, 100, percentile_step):
        distrib['{:02}'.format(curr_percentile) + "_percentile"] = round( numpy.percentile(values, curr_percentile), precision )
    return distrib

def getSelectedAreas( input_panel ):
    """
    @summary: Returns the list of selected areas from a BED file.
    @param input_panel: [str] The path to the amplicons with their primers (format: BED)
    @return: [list] The list of BED's areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "name":"gene_98"}.
    """
    selected_areas = list()
    with open(input_panel) as FH_panel:
        for line in FH_panel:
            if not line.startswith("browser ") and not line.startswith("track ") and not line.startswith("#"):
                fields = [elt.strip() for elt in line.split("\t")]
                selected_areas.append({
                    "region": fields[0],
                    "start": int(fields[1]) +1, # Start in BED is 0-based
                    "end": int(fields[2]),
                    "name": None if len(fields) < 4 else fields[3]
                })
    return( selected_areas )

def setDepths( depths_file, selected_areas ):
    """
    @summary: Adds the list of depths for each area in selected_areas. These depths are stored with the key "data".
    @param depths_file: [str] The path to the depths by position (format: samtools depth output). The file must only contains one sample and it must contains every positions in selected areas (see samtools depth -a option for positions with 0 reads).
    @param selected_areas: [list] The list of area. Each area is represented by a dictionary with this format: { "region":"chr1", "start":501, "end":608, "name":"gene_98" }.
    """
    # List areas by pos ID
    area_by_pos = dict()
    for curr_area in selected_areas:
        pos_id = curr_area["region"] + ":" + str(curr_area["start"])
        if pos_id not in area_by_pos:
            area_by_pos[pos_id] = list()
        area_by_pos[pos_id].append( curr_area )
    # Parse depths file
    with open(depths_file) as FH_depths:
        opened_area = list()
        for line in FH_depths:
            fields = line.strip().split("\t")
            region = fields[0]
            position = int(fields[1]) # Start in depth is 1-based
            depth = int(fields[2])
            # Open areas with start on the current position
            pos_id = region + ":" + str(position)
            if pos_id in area_by_pos:
                for curr_area in area_by_pos[pos_id]:
                    curr_area["data"] = list()
                    opened_area.append( curr_area )
            # Manage position in opened areas
            achieved_idx = list()
            for idx_area, area in enumerate(opened_area):
                if region == area["region"] and position <= area["end"]: # Current position is in area
                    area["data"].append( depth )
                else: # Current position is after the end of the area
                    achieved_idx.append( idx_area )
            # Close achieved areas
            for idx in sorted(achieved_idx, reverse=True):
                del opened_area[idx]

def writeOutputTSV( out_path, area_depths, percentile_step ):
    """
    @summary: Writes depths distribution for each area and each sample in TSV format.
    @param out_path: [str] The path to the outputted file (format: TSV).
    @param area_depths: [list] The list of area. Each area is represented by a dictionary with this format:
                        {
                          "region":"chr1",
                          "start":501,
                          "end":608,
                          "name":"gene_98",
                          "depths":
                            {
                              "splA":{"min":8, "25_percentile":10, "50_percentile":10, "75_percentile":10, "max":30},
                              "splB":{"min":0, "25_percentile":1, "50_percentile":30, "75_percentile":35, "max":50}
                            }
                        }
    @param percentile_step: [int] The step used to restrain distribution (only the depths for this percentile and his multiples have been retained).
    """
    percentiles_titles = ['{:02}'.format(perc) + "_percentile" for perc in range(percentile_step, 100, percentile_step)]
    with open(out_path, "w") as FH_out:
        print(
            "#Region",
            "start",
            "end",
            "name",
            "sample",
            "min_depth",
            "\t".join( [curr_title + "_depth" for curr_title in percentiles_titles] ),
            "max_depth",
            sep="\t", file=FH_out
        )    
        for area in area_depths:
            for spl_name in sorted(area["depths"]):
                print(
                    area["region"],
                    area["start"],
                    area["end"],
                    area["name"],
                    spl_name,
                    area["depths"][spl_name]["min"],
                    "\t".join( [str(area["depths"][spl_name][curr_title]) for curr_title in percentiles_titles] ),
                    area["depths"][spl_name]["max"],
                    sep="\t", file=FH_out
                )

def writeOutputJSON( out_path, area_depths ):
    """
    @summary: Writes depths distribution for each area and each sample in JSON format.
    @param out_path: [str] The path to the outputted file (format: JSON).
    @param area_depths: [list] The list of area. Each area is represented by a dictionary with this format:
                        {
                          "region":"chr1",
                          "start":501,
                          "end":608,
                          "name":"gene_98",
                          "depths":
                            {
                              "splA":{"min":8, "25_percentile":10, "50_percentile":10, "75_percentile":10, "max":30},
                              "splB":{"min":0, "25_percentile":1, "50_percentile":30, "75_percentile":35, "max":50}
                            }
                        }
    """
    with open(out_path, 'w') as FH_out:
        data = json.dumps( area_depths, default=lambda o: o.__dict__, sort_keys=True, indent=2 )
        FH_out.write( data )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Writes depths distribution for specified areas.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-s', '--percentile-step', type=int, default=5, help='Only the depths for this percentile and his multiples are retained. For example, with 25 only the minimum, the 1st quartile, the 2nd quartile, the 3rd quartile and the maximum depths are retained. [Default: %(default)s]' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-c', '--inputs-depths', nargs='+', required=True, help='The path to the depths by position (format: samtools depth output). Each file represents only one sample. The file must contains every positions in selected areas (see samtools depth -a option for positions with 0 reads).' )
    group_input.add_argument( '-r', '--input-regions', required=True, help='Path to the list of regions (format: BED).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-metrics', default="depths_distrib.json", help='The path to outputted file (format: JSON or TSV according to the extension).' )
    args = parser.parse_args()

    # Get interest area
    selected_areas = getSelectedAreas( args.input_regions )

    # Get coverage by area in samples
    for spl_file in args.inputs_depths:
        spl_name = os.path.basename(spl_file)
        setDepths( spl_file, selected_areas )
        for curr_area in selected_areas:
            # Transform depths to distribution
            if "data" not in curr_area:
                curr_area["data"] = [0]
            spl_distrib = getDistributionDict( curr_area["data"], args.percentile_step, 1 )
            del(curr_area["data"])
            # Store distribution
            if "depths" not in curr_area:
                curr_area["depths"] = dict()
            curr_area["depths"][spl_name] = spl_distrib

    # Write output
    if args.output_metrics.endswith("json"):
        writeOutputJSON( args.output_metrics, selected_areas )
    else:
        writeOutputTSV( args.output_metrics, selected_areas, args.percentile_step )
