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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merge variants JSON coming from different sources.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--inputs-variants', nargs='+', required=True, help='Pathes to the variants files (format: JSON).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', help='Path to the outputted variant file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load variants
    log.info("Load and compare variants.")
    variants_by_id = {}
    for curr_json in args.inputs_variants:
        with open(curr_json) as FH_in:
            variants = json.load(FH_in)
            for curr_var in variants:
                curr_var_id = "{}:{}={}/{}".format(
                    curr_var["coord"]["region"],
                    curr_var["coord"]["pos"],
                    curr_var["coord"]["ref"],
                    curr_var["coord"]["alt"]
                )
                if curr_var_id not in variants_by_id:
                    variants_by_id[curr_var_id] = curr_var
                else:
                    for curr_support in curr_var["supports"]:
                        variants_by_id[curr_var_id]["supports"].append(curr_support)

    # Write output
    log.info("Write output.")
    variants = [curr_var for curr_id, curr_var in variants_by_id.items()]
    with open(args.output_variants, "w") as FH_out:
        json.dump(variants, FH_out)
    log.info("END of process.")
