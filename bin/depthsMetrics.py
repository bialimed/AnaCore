#!/usr/bin/env python3
#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import numpy
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sv import SVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def loadFromDepthFile(in_path, samples):
    """
    Load depths classes and count from samtools depth output.

    :param in_path: Path to the samtools depth output.
    :type in_path: str
    :param samples: The list of samples names in order of depths columns.
    :type samples: list
    :return: The list of depths and by sample the list of counts.
    :rtype: list, dict
    """
    encountered_depths = dict()
    count_by_spl = dict()
    with SVIO(in_path, has_title=False) as FH_depths:
        count_by_spl = {curr_spl: dict() for curr_spl in samples}
        for record in FH_depths:  # record = [chr, pos, deph_spl_1, ..., depth_spl_n]
            for spl_idx, curr_spl in enumerate(samples):
                depth = int(record[spl_idx + 2])
                encountered_depths[depth] = 1
                if depth in count_by_spl[curr_spl]:
                    count_by_spl[curr_spl][depth] += 1
                else:
                    count_by_spl[curr_spl][depth] = 1
    depths_list = sorted([key for key in encountered_depths])
    for spl in samples:
        spl_counts = list()
        for depth in depths_list:
            if depth in count_by_spl[spl]:
                spl_counts.append(count_by_spl[spl][depth])
            else:
                spl_counts.append(0)
        count_by_spl[spl] = spl_counts
    return depths_list, count_by_spl


def getDistribution(values, percentile_step=25, precision=4):
    """
    Return the distribution of values (min, max and percentiles).

    :param values: The values.
    :type values: list
    :param percentile_step: Only this percentile and his multiples are returned.
    :type percentile_step: int
    :param precision: The decimal precision.
    :type precision: int
    :retrun: The min, max and percentiles values. Example: {"min":0, "05_percentile":10, "10_percentile":15, ..., "95_percentile":853, "max":859}
    :rtype: dict
    """
    distrib = {
        "min": round(min(values), precision),
        "max": round(max(values), precision)
    }
    for curr_percentile in range(percentile_step, 100, percentile_step):
        distrib['{:02}'.format(curr_percentile) + "_percentile"] = round(numpy.percentile(values, curr_percentile, interpolation="midpoint"), precision)
    return distrib


def getDistribMetrics(depths_list, count_by_spl, percentiles_step):
    """
    Return the distribution of depths by sample.

    :param depths_list: The list of depths.
    :type depths_list: list
    :param count_by_spl: By sample the list of counts corresponding to the depths_list.
    :type count_by_spl: dict
    :param percentiles_step: Only this percentile and his multiples are returned.
    :type percentiles_step: int
    :retrun: Distribution of depths by sample.
    :rtype: dict
    """
    distrib_by_spl = {}
    for spl_name, count_by_depth in count_by_spl.items():
        val_list = []
        for dp, count in zip(depths_list, count_by_depth):
            for idx in range(count):
                val_list.append(dp)
        distrib_by_spl[spl_name] = getDistribution(val_list, percentiles_step)
    return distrib_by_spl


def getSequencingMetrics(depths_list, count_by_spl):
    """
    Return by sample the number of targeted nucleotids and the number of nucleotids sequenced on target. Insertion and deletion are not counted in sequenced on target.

    :param depths_list: The list of depths.
    :type depths_list: list
    :param count_by_spl: By sample the list of counts corresponding to the depths_list.
    :type count_by_spl: dict
    :retrun: By sample the number of targeted nucleotids and the number of nucleotids sequenced on target.
    :rtype: dict
    """
    sequencing_by_spl = {}
    for spl_name, count_by_depth in count_by_spl.items():
        sequencing_by_spl[spl_name] = {"nt_targeted": 0, "nt_sequenced_on_target": 0}
        for dp, count in zip(depths_list, count_by_depth):
            sequencing_by_spl[spl_name]["nt_targeted"] += count
            sequencing_by_spl[spl_name]["nt_sequenced_on_target"] += (count * dp)
    return sequencing_by_spl


def getThresholdMetrics(depths_list, count_by_spl, min_depth):
    """
    Return by sample the count and rate of nucleotids with depth under the threshold.

    :param depths_list: The list of depths.
    :type depths_list: list
    :param count_by_spl: By sample the list of counts corresponding to the depths_list.
    :type count_by_spl: dict
    :param min_depth: Depth threshold.
    :type min_depth: int
    :retrun: By sample the count and rate of nucleotids with depth under the threshold.
    :rtype: dict
    """
    threshold_by_spl = {}
    for spl_name, count_by_depth in count_by_spl.items():
        threshold_by_spl[spl_name] = {"count": 0, "rate": 0}
        total_dp = 0
        for dp, count in zip(depths_list, count_by_depth):
            total_dp += count
            if dp < min_depth:
                threshold_by_spl[spl_name]["count"] += count
                threshold_by_spl[spl_name]["rate"] += count
        threshold_by_spl[spl_name]["rate"] = round(threshold_by_spl[spl_name]["rate"] / total_dp, 4)
    return threshold_by_spl


def writeTSV(out_path, sequencing_by_spl, depths_list, count_by_spl, distrib_by_spl, threshold_by_spl, args):
    """
    Write depths metrics on TSV file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param sequencing_by_spl: By sample the number of targeted nucleotids and the number of nucleotids sequenced on target.
    :type sequencing_by_spl: dict
    :param depths_list: The list of depths.
    :type depths_list: list
    :param count_by_spl: By sample the list of counts corresponding to the depths_list.
    :type count_by_spl: dict
    :param distrib_by_spl: By sample the distribution of depths.
    :type distrib_by_spl: dict
    :param threshold_by_spl: By sample the count and rate of nucleotids with depth under the threshold.
    :type threshold_by_spl: dict
    :param args: The script's parameters.
    :type args: NameSpace
    """
    samples = sorted([spl for spl in count_by_spl])
    with open(out_path, "w") as FH_out:
        # Summary
        FH_out.write(
            "##[SUMMARY]\n#Sample\tNb_targeted_nt\tNb_sequenced_on_target\tUnder_{}X\t%_under_{}X\tMin_depth\t{}\tMax_depth\n".format(
                args.min_depth,
                args.min_depth,
                "\t".join(['{:02}_perc_depth'.format(perc) for perc in range(args.percentiles_step, 100, args.percentiles_step)])
            )
        )
        for spl_name in samples:
            FH_out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    spl_name,
                    sequencing_by_spl[spl_name]["nt_targeted"],
                    sequencing_by_spl[spl_name]["nt_sequenced_on_target"],
                    threshold_by_spl[spl_name]["count"],
                    "{:.2%}".format(threshold_by_spl[spl_name]["rate"]).replace("%", ""),
                    distrib_by_spl[spl_name]["min"],
                    "\t".join([
                        str(distrib_by_spl[spl_name]['{:02}_percentile'.format(perc)]) for perc in range(args.percentiles_step, 100, args.percentiles_step)
                    ]),
                    distrib_by_spl[spl_name]["max"]
                )
            )
        FH_out.write("\n")
        # Details
        FH_out.write(
            "##[DETAILS]\n#Depth\t{}\n".format(
                "\t".join([spl for spl in samples])
            )
        )
        idx_by_depth = {depth: idx for idx, depth in enumerate(depths_list)}
        max_depth = max(depths_list)
        for depth in range(0, max_depth + 1):
            if depth in idx_by_depth:
                FH_out.write(
                    "{}\t{}\n".format(
                        depth,
                        "\t".join([str(count_by_spl[spl][idx_by_depth[depth]]) for spl in samples])
                    )
                )
            else:
                FH_out.write(
                    "{}\t{}\n".format(
                        depth,
                        "\t".join(["0" for spl in samples])
                    )
                )


def writeJSON(out_path, sequencing_by_spl, depths_list, count_by_spl, distrib_by_spl, threshold_by_spl, args):
    """
    Write depths metrics on JSON file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param sequencing_by_spl: By sample the number of targeted nucleotids and the number of nucleotids sequenced on target.
    :type sequencing_by_spl: dict
    :param depths_list: The list of depths.
    :type depths_list: list
    :param count_by_spl: By sample the list of counts corresponding to the depths_list.
    :type count_by_spl: dict
    :param distrib_by_spl: By sample the distribution of depths.
    :type distrib_by_spl: dict
    :param threshold_by_spl: By sample the count and rate of nucleotids with depth under the threshold.
    :type threshold_by_spl: dict
    :param args: The script's parameters.
    :type args: NameSpace
    """
    with open(out_path, "w") as FH_out:
        data = {
            "depths_classes": {
                "depths_list": depths_list,
                "count_by_spl": count_by_spl
            },
            "depths_distrib_by_spl": distrib_by_spl,
            "under_threshold": {
                "count_by_spl": threshold_by_spl,
                "threshold": args.min_depth
            },
            "sequencing_by_spl": sequencing_by_spl
        }
        FH_out.write(json.dumps(data, sort_keys=True))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Produce depths metrics (distribution, count by depth and targets under threshold) from the samtools depth output.')
    parser.add_argument('-m', '--min-depth', type=int, default=30, help='Depth threshold used to return the count of nucleotids under this value.')
    parser.add_argument('-p', '--percentiles-step', type=int, default=5, help='The distribution of depths contains this percentile and his multiples.')
    parser.add_argument('-s', '--samples', required=True, nargs='+', help='The samples names in same order of the depth columns in "--input-depths".')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-depths', required=True, help='The path to the depths file (format: TSV). This have the same format as samtools depth output and must be processed with -aa option.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-j', '--output-json', help='The path to the outputted file (format: JSON).')
    group_output.add_argument('-t', '--output-tsv', help='The path to the outputted file (format: TSV).')
    args = parser.parse_args()

    if args.output_tsv is None and args.output_json is None:
        parser.error('At least one of the output parameters is required.')

    # Process
    depths_list, count_by_spl = loadFromDepthFile(args.input_depths, args.samples)
    sequencing_by_spl = getSequencingMetrics(depths_list, count_by_spl)
    distrib_by_spl = getDistribMetrics(depths_list, count_by_spl, args.percentiles_step)
    threshold_by_spl = getThresholdMetrics(depths_list, count_by_spl, args.min_depth)
    if args.output_tsv is not None:
        writeTSV(args.output_tsv, sequencing_by_spl, depths_list, count_by_spl, distrib_by_spl, threshold_by_spl, args)
    if args.output_json is not None:
        writeJSON(args.output_json, sequencing_by_spl, depths_list, count_by_spl, distrib_by_spl, threshold_by_spl, args)
