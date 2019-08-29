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

from anacore.vcf import VCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='***************************************.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-s', '--samples-names', nargs='+', required=True, help='*************************.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            idx_removed = [idx_spl for idx_spl, spl_name in enumerate(FH_out.samples) if spl_name in args.samples_names]
            for idx in sorted(idx_removed, reverse=True):
                del(FH_out.samples[idx])
            for tag in ["AD", "AF"]:
                if tag in FH_out.info:
                    FH_out.info[tag]["number_tag"] = "A"
                    FH_out.info[tag]["number"] = None
            FH_out.writeHeader()

            # Records
            for record in FH_in:
                for spl_name in args.samples_names:
                    del(record.samples[spl_name])
                has_AD = False
                if "AD" in record.info:
                    has_AD = True
                    del(record.info["AD"])
                has_AF = False
                if "AF" in record.info:
                    has_AF = True
                    del(record.info["AF"])
                has_DP = False
                if "DP" in record.info:
                    has_DP = True
                    del(record.info["DP"])
                if has_DP:
                    record.info["DP"] = record.getPopDP()
                if has_AF:
                    record.info["AF"] = record.getPopAltAF()
                if has_AD:
                    record.info["AD"] = record.getPopAltAD()
                FH_out.write(record)
