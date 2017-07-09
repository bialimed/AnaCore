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

import json
import numpy
import argparse


########################################################################
#
# FUNCTIONS
#
########################################################################
def getDistributionDictForCycles( val_by_cycle, percentile_step=25, precision=4 ):
    """
    @summary: Returns the distribution of values by cycle (min, max and percentiles).
    @param val_by_cycle: [dict] The values by cycle.
    @param percentile_step: [int] Only this percentile and his multiples are returned.
    @param precision: [int] The decimal precision.
    @retrun: [dict] The min, max and percentiles values. Example: {"min":[0, 1, 5], "05_percentile":[10, 15, 15], "10_percentile":[15, 28, 30], ..., "95_percentile":[853, 853, 853], "max":[859, 859, 900]}
    """
    distrib = {
        "min": [round(min(val_by_cycle[cycle]), precision) for cycle in sorted(val_by_cycle)],
        "max": [round(max(val_by_cycle[cycle]), precision) for cycle in sorted(val_by_cycle)]
    }
    for curr_percentile in range(percentile_step, 100, percentile_step):
        distrib['{:02}'.format(curr_percentile) + "_percentile"] = [round(numpy.percentile(val_by_cycle[cycle], curr_percentile), precision) for cycle in sorted(val_by_cycle)]
    return distrib

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

def getrunMetrics( interop_path ):
    """
    @summary: Returns summary of run metrics from InterOp file.
    @param interop_path: [str] "The interop file to process (format: TXT). This file comes from the usage of dumptext on Illumina's interOp folder (cf. package on http://illumina.github.io/interop/index.html).
    @returns: [dict] The density/tile, count/tile and error_rate/tile distribution.
    """
    tiles = {
        "density": list(),
        "density_PF": list(),
        "cluster_count": list(),
        "cluster_count_PF": list(),
        "error_rate": dict()
    }

    # Parse file
    with open(interop_path) as FH:
        is_tiles_section = False
        is_error_section = False
        tiles_header = None
        for line in FH:
            if line.startswith("#") or line.strip() == "":
                is_tiles_section = False
                is_error_section = False
                if line.strip() == "# Tile,1":
                    is_tiles_section = True
                    trash = FH.readline() # Column Count: 10
                    tiles_header = [col.strip() for col in FH.readline().split(",")] # Lane,Tile,Read,ClusterCount,ClusterCountPF,Density,DensityPF,Aligned,Prephasing,Phasing
                elif line.strip() == "# Error,1":
                    is_error_section = True
                    trash = FH.readline() # Column Count: 4
                    tiles_header = [col.strip() for col in FH.readline().split(",")] # Lane,Tile,Cycle,ErrorRate
            elif is_tiles_section:
                tile_metrics = {tiles_header[col_idx]:col.strip() for col_idx, col in enumerate(line.split(","))}
                tiles["density"].append( float(tile_metrics["Density"]) )
                tiles["density_PF"].append( float(tile_metrics["DensityPF"]) )
                tiles["cluster_count"].append( int(float(tile_metrics["ClusterCount"])) )
                tiles["cluster_count_PF"].append( int(float(tile_metrics["ClusterCountPF"])) )
            elif is_error_section:
                tile_metrics = {tiles_header[col_idx]:col.strip() for col_idx, col in enumerate(line.split(","))}
                cycle = tile_metrics["Cycle"]
                if cycle not in tiles["error_rate"]:
                    tiles["error_rate"][cycle] = list()
                tiles["error_rate"][cycle].append( float(tile_metrics["ErrorRate"]) )

    # Process summary
    tiles_summary = {
        "density": getDistributionDict( tiles["density"], 5, 0 ),
        "density_PF": getDistributionDict( tiles["density_PF"], 5, 0 ),
        "cluster_count": getDistributionDict( tiles["cluster_count"], 5, 0 ),
        "cluster_count_PF": getDistributionDict( tiles["cluster_count_PF"], 5, 0 ),
        "error_rate": getDistributionDictForCycles( tiles["error_rate"], 5, 4 )
    }

    return tiles_summary


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Summarizes an Illumina's run metrics from an InterOp file." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--interop-file', required=True, help="The interop file to process (format: TXT). This file comes from the usage of dumptext on Illumina's interOp folder (cf. package on http://illumina.github.io/interop/index.html)." )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', default="interOp_summary.json", help='Path to the output file (format: JSON). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    run_metrics = getrunMetrics( args.interop_file )
    with open(args.output_file, "w") as FH_out:
        FH_out.write( json.dumps(run_metrics, default=lambda o: o.__dict__, sort_keys=True) )
