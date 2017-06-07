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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
import statistics
import argparse


########################################################################
#
# FUNCTIONS
#
########################################################################
def getrunMetrics( interop_path ):
    """
    @summary: Returns summary of run metrics from InterOp file.
    @param interop_path: [str] "The interop file to process (format: TXT). This file comes from the usage of dumptext on Illumina's interOp folder (cf. package on http://illumina.github.io/interop/index.html).
    @returns: [dict] The density/tile and count/tile distribution.
    """
    tiles = {
        "density": list(),
        "density_PF": list(),
        "cluster_count": list(),
        "cluster_count_PF": list(),
    }
    
    # Parse file
    with open(interop_path) as FH:
        is_tiles_section = False
        tiles_header = None
        for line in FH:
            if line.startswith("#") or line.strip() == "":
                tiles_header = False
                if line.strip() == "# Tile,1":
                    is_tiles_section = True
                    trash = FH.readline() # Column Count: 10
                    tiles_header = [col.strip() for col in FH.readline().split(",")] # Lane,Tile,Read,ClusterCount,ClusterCountPF,Density,DensityPF,Aligned,Prephasing,Phasing
            elif is_tiles_section:
                tile_metrics = {tiles_header[col_idx]:col.strip() for col_idx, col in enumerate(line.split(","))}
                tiles["density"].append( float(tile_metrics["Density"]) )
                tiles["density_PF"].append( float(tile_metrics["DensityPF"]) )
                tiles["cluster_count"].append( int(float(tile_metrics["ClusterCount"])) )
                tiles["cluster_count_PF"].append( int(float(tile_metrics["ClusterCountPF"])) )

    # Process summary
    tiles = {
        "density": {
            "min": min(tiles["density"]),
            "median": statistics.median(tiles["density"]),
            "mean": statistics.mean(tiles["density"]),
            "max": max(tiles["density"])
        },
        "density_PF": {
            "min": min(tiles["density_PF"]),
            "median": statistics.median(tiles["density_PF"]),
            "mean": statistics.mean(tiles["density_PF"]),
            "max": max(tiles["density_PF"])
        },
        "cluster_count": {
            "min": min(tiles["cluster_count"]),
            "median": statistics.median(tiles["cluster_count"]),
            "mean": statistics.mean(tiles["cluster_count"]),
            "max": max(tiles["cluster_count"]),
            "sum": sum(tiles["cluster_count"])
        },
        "cluster_count_PF": {
            "min": min(tiles["cluster_count_PF"]),
            "median": statistics.median(tiles["cluster_count_PF"]),
            "mean": statistics.mean(tiles["cluster_count_PF"]),
            "max": max(tiles["cluster_count_PF"]),
            "sum": sum(tiles["cluster_count_PF"])
        },
    }
    
    return tiles


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
