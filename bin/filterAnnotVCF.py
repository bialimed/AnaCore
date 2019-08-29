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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse
from copy import deepcopy

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.filters import filtersFromDict
from anacore.annotVcf import AnnotVCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filters VCF variants and annotations on criteria described in JSON files.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-r', '--input-filters-variants', help='The path to the filters file on variants (format: JSON).')
    group_input.add_argument('-a', '--input-filters-annotations', help='The path to the filters file on annotations (format: JSON).')
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get filters
    record_filters = None
    if args.input_filters_variants is not None:
        with open(args.input_filters_variants) as data_file:
            record_filters = filtersFromDict(json.load(data_file))
    annot_filters = None
    if args.input_filters_annotations is not None:
        with open(args.input_filters_annotations) as data_file:
            annot_filters = filtersFromDict(json.load(data_file))
    if record_filters is None and annot_filters is None:
        parser.error('"--input-filters-annotations" and/or "--input-filters-variants" must be specified.')

    # Process
    with AnnotVCFIO(args.output_variants, mode="w") as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotation_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                if record_filters is None or record_filters.eval(record):
                    if annot_filters is not None:
                        if args.annotation_field in record.info and len(record.info[args.annotation_field]) != 0:
                            kept_annot = []
                            for annot in record.info[args.annotation_field]:
                                upgraded_annot = deepcopy(annot)
                                upgraded_annot["variant"] = record
                                if annot_filters.eval(upgraded_annot):
                                    kept_annot.append(annot)
                            record.info[args.annotation_field] = kept_annot
                    FH_out.write(record)
    log.info("End of job")
