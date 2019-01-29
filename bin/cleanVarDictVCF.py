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
    parser = argparse.ArgumentParser(description='Fix bug in INFO fields format of VarDict outputs.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Rules to clean file
    clean_info = {
        "REFBIAS": {
            "declaration": {"type": int, "type_tag": "Integer", "number_tag": "2", "number": 2, "description": "Reference depth by strand"},
            "process": lambda val: [int(elt) for elt in val.split(":")]
        },
        "VARBIAS": {
            "declaration": {"type": int, "type_tag": "Integer", "number_tag": "2", "number": 2, "description": "Variant depth by strand"},
            "process": lambda val: [int(elt) for elt in val.split(":")]
        }
    }

    # Process
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            for tag in FH_out.info:
                if tag in clean_info:
                    prev = FH_out.info[tag]
                    new = clean_info[tag]["declaration"]
                    if prev["type_tag"] != new["type_tag"] or prev["number_tag"] != new["number_tag"]:
                        FH_out.info[tag] = new
                    else:
                        del(clean_info[tag])
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                for tag, value in record.info.items():
                    if tag in clean_info:
                        record.info["tag"] = clean_info[tag]["process"](value)
                FH_out.write(record)
