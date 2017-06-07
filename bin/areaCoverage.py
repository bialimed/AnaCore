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
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
import numpy
import statistics
import argparse



########################################################################
#
# FUNCTIONS
#
########################################################################
def getInterestArea( interest_bed ):
    """
    #~ interest_area = {
        #~ "6":{
            #~ 26115141: [26123500],
            #~ 26114621: [26114873]
        #~ }
    #~ }
    """
    nb_area = 0
    interest_area = dict()
    with open(interest_bed) as FH_bed:
        for line in FH_bed:
            region, start, end = line.strip().split("\t")[0:3]
            name = "" if len(line.split("\t")) < 4 else line.split("\t")[3].strip()
            start = int(start) + 1 # In BED start is zero-based starting position
            end = int(end) # In BED end is one-based starting position
            if start > end:
                old_start = start
                start = end
                end = old_start
            if region not in interest_area:
                interest_area[region] = dict()
            if start not in interest_area[region]:
                interest_area[region][start] = list()
            interest_area[region][start].append({"end":end, "name":name})
    return interest_area

def getCovByArea( coverage_file, interest_area ):
    """
    @warning: All the positions must exist in coverage_file.
    """
    achieved_area = list()
    with open(coverage_file) as FH_coverage:
        opened_area = list()
        for line in FH_coverage:
            fields = line.strip().split("\t")
            region = fields[0]
            position = int(fields[1])
            coverage = int(fields[2])
            if region in interest_area:
                if position in interest_area[region]:
                    for current_area in interest_area[region][position]:
                        opened_area.append( {"region":region, "start":position, "end":current_area["end"], "name":current_area["name"], "cov":[]} )
            achieved_idx = list()
            for idx_area, area in enumerate(opened_area):
                if region == area["region"] and position <= area["end"]:
                    area["cov"].append( coverage )
                else:
                    achieved_idx.append( idx_area )
            for idx in sorted(achieved_idx, reverse=True):
                achieved_area.append( opened_area[idx] )
                del opened_area[idx]
        for area in opened_area:
            achieved_area.append( area )
    return( achieved_area )

def writeOutputTSV( out_path, area_coverage ):
    with open(out_path, "w") as FH_out:
        print(
            "#Region",
            "start",
            "end",
            "name",
            "cov_min",
            "cov_lower_quartile",
            "cov_median",
            "cov_upper_quartile",
            "cov_max",
            sep="\t", file=FH_out
        )    
        for area in area_coverage:
            print(
                area["region"],
                area["start"],
                area["end"],
                area["name"],
                min(area["cov"]),
                numpy.percentile( area["cov"], 25 ),
                numpy.percentile( area["cov"], 50 ),
                numpy.percentile( area["cov"], 75 ),
                max( area["cov"] ),
                sep="\t", file=FH_out
            )

def writeOutputJSON( out_path, area_coverage ):
    with open(out_path, "w") as FH_out:
        for area in area_coverage:
            area["coverage"] = {
                "min": min(area["cov"]),
                "1st_quart": numpy.percentile( area["cov"], 25 ),
                "median": statistics.median( area["cov"] ),
                "3rd_quart": numpy.percentile( area["cov"], 75 ),
                "max": max(area["cov"]) }
            del( area["cov"] )
        FH_out.write( json.dumps(area_coverage, default=lambda o: o.__dict__, sort_keys=True ) )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Coverage distribution for specied area in alignment." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-s', '--input-selection', required=True, help='Path to the interests area description (format: BED).' )
    group_input.add_argument( '-d', '--input-depth', required=True, help='Path to the samtools depth output (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output', default="areaDepth.json", help='Path to the output (format: JSON or TSV). [Default: %(default)s]' )
    args = parser.parse_args()

    # Get interest area
    interest_area = getInterestArea( args.input_selection )

    # Get coverage by pos in area
    coverage = getCovByArea( args.input_depth, interest_area )
    coverage = sorted( coverage, key=lambda x: (x["region"], x["start"], x["end"]) )

    # Write output
    if args.output.endswith("json"):
        writeOutputJSON( args.output, coverage )
    else:
        writeOutputTSV( args.output, coverage )
