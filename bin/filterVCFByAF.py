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

from anacore.vcf import VCFIO, HeaderFilterAttr


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter variants on their AF.')
    parser.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag "CSQ" and/or "popAF" is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]')
    parser.add_argument('-a', '--min-AF', default=0.02, type=float, help='Filter variants with AF <= than this values. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', help='Path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', help='Path to the file outputted file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_variants = 0
    nb_filtered = 0
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            if "ADSRC" in FH_in.format:  # Variants file come from merging of different calling
                FH_out.filter["lowAF"] = HeaderFilterAttr("lowAF", "Variants with population AF <= {} in all variant calling sources".format(args.min_AF))
            else:
                FH_out.filter["lowAF"] = HeaderFilterAttr("lowAF", "Variants with population AF <= {}".format(args.min_AF))
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                if len(record.alt) > 1:
                    raise Exception("The multi-allelic variants cannot be processed: {}.".format(record.getName()))
                nb_variants += 1
                is_filtered = False
                if "ADSRC" in FH_in.format:  # Variants file come from merging of different calling
                    pop_AD = [0 for src in record.info["SRC"]]  # Population AD for each calling source
                    pop_DP = [0 for src in record.info["SRC"]]  # Population DP for each calling source
                    for spl_name, spl_info in record.samples.items():
                        for idx, (AD, DP) in enumerate(zip(spl_info["ADSRC"], spl_info["DPSRC"])):
                            pop_AD[idx] += AD
                            pop_DP[idx] += DP
                    pop_AF = [AD / DP for AD, DP in zip(spl_info["ADSRC"], spl_info["DPSRC"])]
                    if max(pop_AF) < args.min_AF:
                        is_filtered = True
                else:
                    if record.getPopAltAF()[0] < args.min_AF:
                        is_filtered = True
                # Filter record
                if args.mode == "remove":
                    if is_filtered:
                        nb_filtered += 1
                    else:
                        FH_out.write(record)
                else:
                    if record.filter is None or len(record.filter) == 0 or record.filter[0] == "PASS":
                        record.filter = list()
                    if is_filtered:
                        record.filter.append("lowAF")
                        nb_filtered += 1
                    if len(record.filter) == 0:
                        record.filter.append("PASS")
                    FH_out.write(record)

    # Log process
    log.info(
        "{:.2%} of variants have been {} ({}/{})".format(
            0 if nb_variants == 0 else nb_filtered / nb_variants,
            "tagged" if args.mode == "tag" else "removed",
            nb_filtered,
            nb_variants
        )
    )
    log.info("End of job")
