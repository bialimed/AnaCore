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
import sys
import time
import json
import argparse
import subprocess

BIN_DIR = os.path.dirname(__file__)
os.environ['PATH'] = os.environ['PATH'] + os.pathsep + BIN_DIR



########################################################################
#
# FUNCTIONS
#
########################################################################
class TmpFiles:
    """
    @summary: Manager for temporary files.
    @copyright: FROGS's team INRA.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)


def distToHC(in_distances, out_tree, linkage_method, tree_format):
    """
    @summary: Builts hierarchical clustering from the distance matrix.
    @param in_distances: [str] Path to the file containing the distance matrix (format: TSV).
    @param out_tree: [str] Path to the hierarchical clustering tree (format: see --output-format).
    @param linkage_method: [str] Used linkage (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage).
    @param: tree_format: [str] The tree format: "newick" or "json" or "png".
    """
    cmd = [
        "distToHC.py",
        "--linkage-method", linkage_method,
        "--input-distances", in_distances,
        "--output-tree", out_tree,
        "--output-format", tree_format
    ]
    subprocess.check_output(cmd)


def filterVCF(in_variants, out_variants, min_AF):
    """
    @summary: Filters variants on their frequencies. Only variants with at least one sample where his AF >= min_AF are kept.
    @param in_variants: [str] The path to the variants file (format: VCF).
    @param out_variants: [str] The path to the filtered variants file (format: VCF).
    @param min_AF: [float] The minimum allele frequency in at least one sample.
    """
    tmp = TmpFiles(os.path.dirname(out_variants))
    # Writes filters rules
    filters_desc = {
      "class": "Filter",
      "getter": "m:getAFBySample(0).m:values.i:.0",
      "action": "select",
      "aggregator": "nb:1",
      "operator": ">=",
      "values": min_AF
    }
    filters_file = tmp.add("filters.json")
    with open(filters_file, "w") as FH_filters:
        json.dump(filters_desc, FH_filters)
    # Filter
    cmd = [
        "filterVCF.py",
        "--input-filters", filters_file,
        "--input-variants", in_variants,
        "--output-variants", out_variants
    ]
    subprocess.check_output(cmd)
    # Delete temporary
    tmp.deleteAll()

def checkGroups(in_variants, in_groups, out_distances, out_intruders, distance_method):
    """
    @summary: Writes the list of samples more similar with samples from another group than others samples coming from the same group.
    @param in_variants: [str] Path to the variants file containing samples to compare (format: VCF).
    @param in_groups: [str] Path to the file describing links between samples and groups (format: TSV).
    @param out_distances: [str] Path to the 2D distances matrix (format: TSV).
    @param out_intruders: [str] Path to the file listing the samples more similar with samples from another group than others samples coming from the same group. (format: TSV or JSON depends on extension).
    @param distance_method: [str] Used distance (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html#module-scipy.spatial.distance).
    """
    cmd = [
        "checkGroupsVCF.py",
        "--distance-method", distance_method,
        "--input-variants", in_variants,
        "--input-groups", in_groups,
        "--output-distances", out_distances,
        "--output-invalids", out_intruders
    ]
    subprocess.check_output(cmd)


def mergeVCF(in_variants, out_variants):
    """
    @summary: Merges samples in variants files.
    @param in_variants: [list] Path to the variants files.
    @param out_variants: [str] Path to the merged file.
    """
    cmd = [
        "mergeVCF.py",
        "--deactivate-completion",
        "--output-variants", out_variants,
        "--input-variants"
    ]
    cmd.extend(in_variants)
    subprocess.check_output(cmd)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Returns the hierarchical clustering of samples coming from several VCF and the evaluation of their groups cohesion. These process are based on variants alleles frequencies.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-d', '--distance-method', default="euclidean", help='Used distance (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html#module-scipy.spatial.distance). [Default: %(default)s]')
    parser.add_argument('-l', '--linkage-method', default="ward", help='Used linkage (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage). [Default: %(default)s]')
    parser.add_argument('-m', '--min-AF', type=float, default=0.0, help='Only variants with at least one sample with an allele frequency superior or equal to this value are evaluated. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', nargs='+', required=True, help='The path to the variants files (format: VCF).')
    group_input.add_argument('-g', '--input-groups', required=True, help='Path to the file describing links between samples and groups (format: TSV).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-invalids', required=True, help='Path to the file listing the samples more similar with samples from another group than others samples coming from the same group. (format: see --output-format).')
    group_output.add_argument('-t', '--output-tree', required=True, help='Path to the hierarchical clustering tree (format: see --output-format).')
    group_output.add_argument('-f', '--output-format', default="human", choices=["human", "computer"], help='The output format. With human invalids are in TSV and tree in png ; With computer invalids and computer are in JSON. [Default: %(default)s]')
    args = parser.parse_args()

    tmp = TmpFiles(os.path.dirname(args.output_invalids))

    # Merge samples variants
    merged_vcf = args.inputs_variants[0]
    if len(args.inputs_variants) > 1:
        merged_vcf = tmp.add("merged.vcf")
        mergeVCF(args.inputs_variants, merged_vcf)

    # Filter variants
    filtered_vcf = merged_vcf
    if args.min_AF > 0.0:
        filtered_vcf = tmp.add("filtered.vcf")
        filterVCF(merged_vcf, filtered_vcf, args.min_AF)

    # List samples with intruders in group
    dist_tsv = tmp.add("distances.tsv")
    out_invalids = tmp.add("invalid.tsv") if args.output_format == "human" else tmp.add("invalid.json")
    checkGroups(filtered_vcf, args.input_groups, dist_tsv, out_invalids, args.distance_method)
    os.rename(out_invalids, args.output_invalids)

    # Built hierarchical clustering
    tree_format = "png" if args.output_format == "human" else "json"
    distToHC(dist_tsv, args.output_tree, args.linkage_method, tree_format)

    # Delete temporary
    tmp.deleteAll()
