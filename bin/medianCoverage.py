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
import os
import sys
import statistics


def getInterestArea( interest_bed ):
    nb_area = 0
    interest_area = dict()
    with open(sys.argv[1]) as FH_bed:
        for line in FH_bed:
            nb_area += 1
            region, start, end = line.strip().split("\t")
            start = int(start) + 1 ##################################################### Pb start samtools depth
            end = int(end)
            if region not in interest_area:
                interest_area[region] = dict()
            if start not in interest_area[region]:
                interest_area[region][start] = list()
            interest_area[region][start].append( end )
    return     nb_area, interest_area

def getCovByArea( coverage_file, interest_area ):
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
                    for end in interest_area[region][position]:
                        opened_area.append( {"region":region, "start":position, "end":end, "cov":[]} )
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




# Get interest area
nb_area, interest_area = getInterestArea( sys.argv[1] )
#~ interest_area = {
    #~ "6":{
        #~ 26115141: [26123500],
        #~ 26114621: [26114873]
    #~ }
#~ }
# check start < end

# Get coverage by pos in area
cov_by_spl = dict()
for coverage_file in sys.argv[2:]:
    area_list = getCovByArea( coverage_file, interest_area )
    cov_by_spl[os.path.basename(coverage_file)] = sorted(area_list, key=lambda x: (x["start"], x["end"]))
# Display
print( "Area", "\t".join(sorted(cov_by_spl)), sep="\t" )
for area_idx in range(nb_area):
    area_id = ""
    coverages = []
    for spl in sorted(cov_by_spl):
        area_data = cov_by_spl[spl][area_idx]
        area_id = area_data["region"] + ":" + str(area_data["start"]) + "-" + str(area_data["end"])
        coverages.append( statistics.median(area_data["cov"]) )
    print( area_id, "\t".join(map(str, coverages)), sep="\t" )


#~ for area in sorted(area_list, key=lambda x: (x["start"], x["end"])):
    #~ print(
        #~ area["region"] + ":" + str(area["start"]) + "-" + str(area["end"]),
        #~ statistics.median(area["cov"]),
        #~ sep="\t"
    #~ )
