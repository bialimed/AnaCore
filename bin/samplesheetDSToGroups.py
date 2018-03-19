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
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.illumina import SampleSheetIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Creates a groups description file from an Illumina's samplesheet for amplicon double strands protocol. This file list describes the link between samples and groups. In amplicon double strand protocol each group contains the two libraries coming from the same original sample.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-s', '--selected-manifests', nargs='+', help='The list of manifests file retained in extraction. Only the samples corresponding to these manifests are processed. [Default: all manifests]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-samplesheet', required=True, help="Path to the sheet describing run and samples (format: Illumina's samplesheet).")
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-groups', required=True, help='Path to the output file (format: TSV).')
    args = parser.parse_args()

    # Process
    samplesheet = SampleSheetIO(args.input_samplesheet)
    selected_manifests = args.selected_manifests
    if args.selected_manifests is None:
        selected_manifests = samplesheet.manifests.values()
    with open(args.output_groups, "w") as FH_out:
        FH_out.write("#Sample\tGroup\tManifest\n")
        for spl_idx, spl in enumerate(samplesheet.samples):
            manifest = samplesheet.manifests[spl["Manifest"]]
            if manifest in selected_manifests:
                spl_name = spl["Sample_Name"].replace("_", "-").replace(" ", "-")
                lib_name = spl_name + "_S" + str(spl_idx + 1)
                FH_out.write(
                    "{}\t{}\t{}\n".format(lib_name, spl_name, manifest)
                )
