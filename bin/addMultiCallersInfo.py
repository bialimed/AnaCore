#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add sources and lowAF filter on each VCF coming from different variant caller.')
    parser.add_argument('-l', '--lowAF-tag', default="lowAF", help='FILTER tag used for variants with population AF <= than --min-AF.')
    parser.add_argument('-m', '--min-AF', type=float, help='Variants with population AF <= than this values in all calllers are tagged. [Default: No tag]')
    parser.add_argument('-c', '--calling-sources', nargs='+', help='Name of callers corresponding to the inputs-variants.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', nargs='+', help='Pathes to the variants files coming from different callers (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--outputs-variants', nargs='+', help='Pathes to the file outputted files (format: JSON).')
    args = parser.parse_args()
    if len(args.calling_sources) != len(args.inputs_variants) or len(args.inputs_variants) != len(args.outputs_variants):
        parser.error("Each input VCF must correspond to one calling source and one output.")

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get sources and AF
    info_by_var = {}
    for idx_in, curr_in in enumerate(args.inputs_variants):
        with VCFIO(curr_in) as FH_in:
            sorted_spl = sorted(FH_in.samples)
            for record in FH_in:
                if len(record.alt) > 1:
                    raise Exception("The script {} cannot be used on VCF containing multi-allelic variant.".format(sys.argv[0]))
                variant_name = record.getName()
                if variant_name not in info_by_var:
                    info_by_var[variant_name] = {
                        "src": [args.calling_sources[idx_in]],
                        "max_AF": record.getPopAltAF()[0]
                    }
                else:
                    info_by_var[variant_name]["src"].append(args.calling_sources[idx_in])
                    info_by_var[variant_name]["max_AF"] = max(info_by_var[variant_name]["max_AF"], record.getPopAltAF()[0])

    # Write
    for idx_in, curr_in in enumerate(args.inputs_variants):
        with VCFIO(args.outputs_variants[idx_in], "w") as FH_out:
            with VCFIO(curr_in) as FH_in:
                # Header
                FH_out.copyHeader(FH_in)
                FH_out.info["SRC"] = {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Variant callers where the variant is identified"}
                if args.min_AF is not None:
                    FH_out.filter[args.lowAF_tag] = "Variants with population AF <= {} in all variant callers".format(args.min_AF)
                FH_out.writeHeader()
                # Records
                for record in FH_in:
                    variant_name = record.getName()
                    record.info["SRC"] = info_by_var[variant_name]["src"]
                    if args.min_AF is not None and info_by_var[variant_name]["max_AF"] < args.min_AF:
                        record.filter.append(args.lowAF_tag)
                    FH_out.write(record)

    log.info("End of job")
