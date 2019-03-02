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

from anacore.bed import BEDIO, getAreas
from anacore.region import consolidated


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add padding to each regions and merge the overlapping and optionally the contiguous regions.')
    parser.add_argument('-c', '--concatenate-names', action="store_true", help="If this option is used the name of a merged region will be the concatenation of all names of the merged regions.")
    parser.add_argument('-m', '--merge-contiguous', action="store_true", help="If this option is used the contiguous regions are merged.")
    parser.add_argument('-p', '--padding-size', type=int, default=50, help="The size adding before and after each regions. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-regions', required=True, help='Path to the file containing the regions to extend.')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-regions', required=True, help='Path to the extended regions (format: BED).')
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(level=logging.INFO, format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info(" ".join(sys.argv))

    # Get regions
    bed_nb_col = BEDIO.getMaxNbCol(args.input_regions)
    if bed_nb_col > 6:
        log.warning("The BED file {} contains {} columns all the columns after the 6th will be removed on output.".format(args.input_regions, bed_nb_col))
    regions = getAreas(args.input_regions)

    # Add padding
    for curr_region in regions:
        curr_region.start = max(1, curr_region.start - args.padding_size)
        curr_region.end += args.padding_size  # The max reference size is not managed

    # Merge overlapping
    consolidated_regions = consolidated(regions, args.merge_contiguous, True)

    # Write output
    out_nb_col = (6 if bed_nb_col >= 6 else bed_nb_col)
    with BEDIO(args.output_regions, "w", write_nb_col=out_nb_col) as FH_out:
        for curr_region in consolidated_regions:
            if args.concatenate_names and "merge_traceback" in curr_region.annot:
                curr_region.name = ";".join([str(elt) for elt in curr_region.annot["merge_traceback"]])
            FH_out.write(curr_region)

    log.info("End of process")
