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

import os
import sys
import json
import argparse
import warnings
import numpy as np
from scipy.spatial.distance import pdist, squareform

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import getFreqMatrix
from anacore.sv import HashedSVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGroupsData(groups_path, samples, sample_tag="Sample", group_tag="Group", separator="\t"):
    """
    @summary: Return group name by sample, samples by grop and samples without group from separated value file.
    @param groups_path: [str] Path to the separated value file describing links between samples and groups.
    @param samples: [lsit] The list of all samples.
    @param sample_tag: [str] The title of column used to store the samples names.
    @param group_tag: [str] The title of column used to store the groups names.
    @param separator: [str] The separator used between fields in the input file.
    @return: [list] The first element is a dictionary representing the group name by sample. the second element is a dictionary representing the list of samples by group name. The last element is the list of samples without group.
    """
    group_by_spl = {}
    spl_by_group = {}
    without_group = {}
    processed_by_spl = {spl: False for spl in samples}
    # Parse groups information
    with HashedSVIO(groups_path, separator=separator, title_starter="#") as FH_gp:
        for record in FH_gp:
            sample = record[sample_tag]
            group = record[group_tag]
            processed_by_spl[sample] = True
            if sample not in processed_by_spl:
                raise Exception('The sample "{}" found in {} does not exist in expected samples.'.format(sample, groups_path))
            group_by_spl[sample] = group
            if group in spl_by_group:
                spl_by_group[group].append(sample)
            else:
                spl_by_group[group] = [sample]
    # Store samples without group
    for spl, is_in_gp in processed_by_spl.items():
        if not is_in_gp:
            without_group[spl] = True
    # Return
    return group_by_spl, spl_by_group, without_group

def writeInvalidsJSON(intruders_by_spl, group_by_spl, out_path):
    """
    @summary: Writes the list of samples with invalid group in JSON file.
    @param intruders_by_spl: [dict] By sample the list of samples coming from others groups and more similar than others samples coming from the valid group.
    @param group_by_spl: [dict] The group name by sample.
    @param out_path: [str] The path to the output file.
    """
    # Add group information
    json_data = dict()
    for spl, intruders in intruders_by_spl.items():
        json_data[spl] = {
            "expected_group": group_by_spl[spl],
            "intruders": intruders
        }
    # Writes output
    with open(out_path, "w") as FH_invalid:
        FH_invalid.write(
            json.dumps(json_data, default=lambda o: o.__dict__, sort_keys=True)
        )

def writeInvalidsSV(intruders_by_spl, group_by_spl, out_path, separator="\t"):
    """
    @summary: Writes the list of samples with invalid group in separated values file.
    @param intruders_by_spl: [dict] By sample the list of samples coming from others groups and more similar than others samples coming from the valid group.
    @param group_by_spl: [dict] The group name by sample.
    @param out_path: [str] The path to the output file.
    @param separator: [str] The separator used in output file to separate fields.
    """
    with open(out_path, "w") as FH_invalid:
        FH_invalid.write(
            separator.join(["#Sample", "Group", "Intruders\n"])
        )
        for spl in sorted(intruders_by_spl):
            FH_invalid.write(
                ("{}" + separator + "{}" + separator + "{}\n").format(
                    spl,
                    group_by_spl[spl],
                    ", ".join([elt.replace(",", ".") for elt in intruders_by_spl[spl]])
                )
            )

def writeDistances(dist_matrix, compared_elements, out_path, separator="\t"):
    """
    @summary: Writes distance matrix in separated values file.
    @param dist_matrix: [list]
    @param compared_elements: [list] The elements compared in distance matrix.
    @param out_path: [str] The path to the output file.
    @param separator: [str] The separator used in output file to separate fields.
    """
    with open(out_path, "w") as FH_dist:
        FH_dist.write(
            "{}{}\n".format(
                separator,
                separator.join([str(curr_col) for curr_col in compared_elements])
            )
        )
        for row_idx, row_val in enumerate(dist_matrix):
            row = [compared_elements[row_idx]]
            row.extend(row_val)
            FH_dist.write(
                separator.join([str(elt) for elt in row]) + "\n"
            )

def getIntruders(first_spl, samples, distances, group_by_spl, spl_by_group, excluded_spl):
    """
    @summary: Returns for first_spl the list of samples coming from others groups and more similar than others samples coming from the valid group.
    @param first_spl: [str] The evaluated sample name.
    @param samples: [list] The list of samples names in same order as distances parameters.
    @param distances: [list] The list of distances to first_spl in same order as samples parameters.
    @param group_by_spl: [dict] The group name by sample.
    @param spl_by_group: [dict] The list of samples by group name.
    @param excluded_spl: [list] The list of samples without group.
    @return: [list] Samples coming from others groups and more similar to first_spl than others samples coming from the same group.
    """
    intruders = []
    excluded_spl_dict = {elt:1 for elt in excluded_spl}
    first_spl_group = group_by_spl[first_spl]
    # Creates the list of samples by distance
    spl_by_dist = {}
    for idx, dist in enumerate(distances):
        second_spl = samples[idx]
        if second_spl != first_spl and second_spl not in excluded_spl_dict:
            if dist in spl_by_dist:
                spl_by_dist[dist].append(second_spl)
            else:
                spl_by_dist[dist] = [second_spl]
    # Sort samples by distance to first_spl and check if it exists samples in others group more similar at first_spl than farest sample coming from the same group
    nb_gp_members_found = 1
    nb_gp_members_expected = len(spl_by_group[first_spl_group])
    ordered_dist = sorted(spl_by_dist)
    dist_idx = 0
    while nb_gp_members_found != nb_gp_members_expected:
        curr_dist = ordered_dist[dist_idx]
        for spl in spl_by_dist[curr_dist]:
            if group_by_spl[spl] == first_spl_group:  # First sample and second sample are in same group
                nb_gp_members_found += 1
            else:
                intruders.append(spl)
        dist_idx += 1
    return intruders



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Returns the list of samples more similar with samples from another group than others samples coming from the same group.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-d', '--distance-method', default="euclidean", help='Used distance (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html#module-scipy.spatial.distance). [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file containing samples to compare (format: VCF).')
    group_input.add_argument('-g', '--input-groups', required=True, help='Path to the file describing links between samples and groups (format: TSV).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-invalids', required=True, help='Path to the file listing the samples more similar with samples from another group than others samples coming from the same group. (format: TSV or JSON depends on extension).')
    group_output.add_argument('-od', '--output-distances', help='Path to the 2D distances matrix (format: TSV).')
    args = parser.parse_args()

    # Get variant matrix from samples
    samples, variants, AF_matrix = getFreqMatrix(args.input_variants)
    AF_matrix = [list(x) for x in zip(*AF_matrix)]  # Transpose matrix

    # Process distances
    dist_matrix = None
    if len(samples) == 1:
        dist_matrix = np.zeros((1, 1))
    else:
        dist_matrix = pdist(AF_matrix, args.distance_method)
        nb_nan = np.count_nonzero(np.isnan(dist_matrix))
        if nb_nan > 0:
            nb_ok = np.count_nonzero(~np.isnan(dist_matrix))
            warnings.warn(
                "The {} distance on AF generates {}/{} NaN values.".format(args.distance_method, nb_nan, nb_ok)
            )
        dist_matrix = squareform(dist_matrix)

    # Evaluate groups
    intruders_by_spl = dict()
    group_by_spl, spl_by_group, without_group = getGroupsData(args.input_groups, samples)
    for row_idx, row in enumerate(dist_matrix):
        first_spl = samples[row_idx]
        intruders = getIntruders(first_spl, samples, row, group_by_spl, spl_by_group, without_group)
        if len(intruders) > 0:
            intruders_by_spl[first_spl] = intruders

    # Write outputs
    if args.output_invalids.endswith("json"):
        writeInvalidsJSON(intruders_by_spl, group_by_spl, args.output_invalids)
    else:
        writeInvalidsSV(intruders_by_spl, group_by_spl, args.output_invalids)
    if args.output_distances:
        writeDistances(dist_matrix, samples, args.output_distances)
