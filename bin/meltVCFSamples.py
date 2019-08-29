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
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, getAlleleRecord


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Melts all the samples contained in VCF in one sample. For the output only the fields AD, AF and DP will be writted in FORMAT.')
    parser.add_argument('-f', '--AF-precision', type=float, default=5, help="The AF's decimal precision. [Default: %(default)s]")
    parser.add_argument('-n', '--new-spl-name', default="all", help="The name of the outputted sample. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file with several samples (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.filter = FH_in.filter
            removed_from_info = list()
            FH_out.info = FH_in.info.copy()
            for field in ["AD", "AF", "DP"]:
                # Remove population AD, AF and DP
                if field in FH_out.info:
                    removed_from_info.append(field)
                    del(FH_out.info[field])
                # Create new FORMAT
                if field in FH_in.format:
                    FH_out.format[field] = FH_in.format[field]
                    if field != "DP":
                        FH_out.format[field]["number"] = None
                        FH_out.format[field]["number_tag"] = "A"
                else:
                    FH_out.format[field] = {
                        "type": int if field != "AF" else float,
                        "type_tag": "Integer" if field != "AF" else "Float",
                        "number": None,
                        "number_tag": "A",
                        "description": ""
                    }
            FH_out.samples = [args.new_spl_name]
            FH_out.writeHeader()

            # Records
            for record in FH_in:
                pop = {
                    "AD": None,
                    "DP": record.getPopDP(),
                    "AF": [round(curr_AF, args.AF_precision) for curr_AF in record.getPopAltAF()]
                }
                pop["AD"] = [int(curr_AF * pop["DP"]) for curr_AF in pop["AF"]]
                # Replace samples by the new sample
                record.samples = {args.new_spl_name: pop}
                record.format = ["AF", "AD", "DP"]
                # Remove old population AD, AF and DP
                for field in removed_from_info:
                    if field in record.info:
                        del(record.info[field])
                FH_out.write(record)
