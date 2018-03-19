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
import json
import argparse
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform, is_valid_y

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.node import Node
from anacore.matrix import DistanceMatrixIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Builts hierarchical clustering from the distance matrix.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-l', '--linkage-method', default="ward", help='Used linkage (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage). [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-distances', required=True, help='Path to the file containing the distance matrix (format: TSV).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-tree', required=True, help='Path to the hierarchical clustering tree (format: see --output-format).')
    group_output.add_argument('-f', '--output-format', default="newick", choices=["newick", "json", "png"], help='The output format. [Default: %(default)s]')
    args = parser.parse_args()

    # Load distance matrix
    dist_matrix_io = DistanceMatrixIO(args.input_distances)
    dist_matrix = dist_matrix_io.dist_matrix
    rows_names = dist_matrix_io.names

    # Process tree
    tree = None
    data_link = None
    if len(rows_names) == 1:
        tree = Node(rows_names[0])
    else:
        # Computing distance and linkage
        if not is_valid_y(dist_matrix):
            dist_matrix = squareform(dist_matrix)
        data_link = linkage(dist_matrix, args.linkage_method)
        # SciPy format to Node
        hc_tree = to_tree(data_link, rd=False)
        id_2_name = dict(zip(range(len(rows_names)), rows_names))
        tree = Node.fromClusterNode(hc_tree, id_2_name)

    # Write output
    if args.output_format != "png":  # Text outputs
        out_str = None
        if args.output_format == "newick":
            out_str = "{};".format(tree.toNewick())
        elif args.output_format == "json":
            out_str = json.dumps(tree.toDict(), default=lambda o: o.__dict__, sort_keys=False)
        with open(args.output_tree, "w") as FH_out:
            FH_out.write(out_str)
    else:  # Image output
        import matplotlib
        matplotlib.use('Agg')  # Forces matplotlib to not use any Xwindows backend
        import matplotlib.pyplot as plot
        from scipy.cluster.hierarchy import dendrogram
        if len(rows_names) >= 1:
            dendro = dendrogram(data_link, labels=rows_names, orientation="left")
        plot.tight_layout()  # Adjusts the location of axes to prevent cuts in labels
        plot.savefig(args.output_tree)
