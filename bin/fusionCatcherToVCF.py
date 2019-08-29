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

from anacore.annotVcf import AnnotVCFIO
from anacore.fusionCatcher import FusionCatcherIO, setFusionCatcherVCFHeader


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Converts the TSV output of fusionCatcher to VCF.")
    parser.add_argument('-s', '--sample-name', required=True, help='The sample name.')
    parser.add_argument('-a', '--annotation-field', default="FCANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-fusions', required=True, help='Path to the output of fusionCatcher (format: TSV).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-fusions', required=True, help='Path to the output file (format: VCF).')
    args = parser.parse_args()

    # Process
    with AnnotVCFIO(args.output_fusions, "w", args.annotation_field) as FH_out:
        # Header
        setFusionCatcherVCFHeader(FH_out, args.annotation_field)
        FH_out.samples = [args.sample_name]
        FH_out.writeHeader()
        # Records
        with FusionCatcherIO(args.input_fusions, "r", args.sample_name, args.annotation_field) as FH_in:
            for first_bnd, second_bnd in FH_in:
                FH_out.write(first_bnd)
                FH_out.write(second_bnd)
