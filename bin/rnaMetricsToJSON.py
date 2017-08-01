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
import json
import argparse


########################################################################
#
# FUNCTIONS
#
########################################################################
def parseMetricsFile( rnaSeqMetrics_output ):
    metrics_data = dict()
    with open(rnaSeqMetrics_output) as FH_metrics:
        section_metrics = False
        metrics_titles = None
        for line in FH_metrics:
            if line.strip() == "":
                section_metrics = False
            elif line.startswith( "## METRICS CLASS" ):
                section_metrics = True
            elif section_metrics:
                if metrics_titles is None:
                    metrics_titles = [field.strip().lower() for field in line.split("\t")]
                else:
                    # Get metrics
                    data = dict()
                    for field_idx, field in enumerate(line.split("\t")):
                        if not metrics_titles[field_idx].startswith("pct_"):
                            data[metrics_titles[field_idx]] = field.strip() 
                    # Cast numbers
                    for title in data:
                        if title not in ["sample", "library", "read_group"]:
                            if title.startswith("median"):
                                data[title] = float(data[title])
                            else:
                                data[title] = int(data[title])
                    # Add data
                    if data["sample"] == "":
                        data["sample"] = os.path.basename( rnaSeqMetrics_output ).split("_")[0]
                    metrics_data[data["sample"]] = data
    return metrics_data


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Converts picard RnaSeqMetrics output in JSON format." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-metrics', required=True, help='Path to the output of picard RnaSeqMetrics (format: TXT).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-metrics', required=True, help='path to metrics in JSON format (format: JSON).' )
    args = parser.parse_args()

    # Process
    metrics_data = parseMetricsFile( args.input_metrics )
    with open(args.output_metrics, "w") as FH_out:
        FH_out.write( json.dumps(metrics_data, default=lambda o: o.__dict__, sort_keys=True ) )
