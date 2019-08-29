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

from anacore.vcf import VCFIO, getAlleleRecord


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
    parser.add_argument('-g', '--group-name', help='*************************.')
    group_thresholds = parser.add_argument_group('Thresolds')  # Thresholds
    group_thresholds.add_argument('-ad', '--AD-threshold', type=int, default=4, help='*************************.')
    group_thresholds.add_argument('-af', '--AF-threshold', type=float, default=0.05, help='*************************.')
    group_thresholds.add_argument('-rt', '--nb-spl-ratio-threshold', type=float, default=round(2/float(3), 2), help='*************************.')
    group_thresholds.add_argument('-odp', '--others-DP-threshold', type=int, default=5, help='*************************.')
    group_thresholds.add_argument('-ovd', '--others-valid-DP', type=int, default=1, help='*************************.')
    group_thresholds.add_argument('-mc', '--max-contradict', type=int, default=0, help='*************************.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    spl_in_gp = { spl:1 for spl in args.samples_names }
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.info["SUP"] = { "type": int, "type_tag": "Integer", "number": 1, "number_tag": "1", "description": "The number of samples in group supporting variant (AD >= " + str(args.AD_threshold) + " and AF >= " + str(args.AF_threshold) + "." }
            FH_out.info["CONT"] = { "type": int, "type_tag": "Integer", "number": 1, "number_tag": "1", "description": "The number of samples out group supporting variant (AD >= " + str(args.AD_threshold) + " or AF >= " + str(args.AF_threshold) + "." }
            FH_out.filter["SPEC"] = "The variant is not specific at the group" + ("" if args.group_name is None else " " + args.group_name) + " (" + ", ".join(args.samples_names) + "). The variant is specific if: (1) the ratio number_of_support_samples_in_group/number_of_samples_in_group is >= " + str(args.nb_spl_ratio_threshold) + " ; (2) the number of support samples out of the group is <= " + str(args.max_contradict) + " ; (3) at least " + str(args.others_valid_DP) + " samples out of the group have DP >= " + str(args.others_DP_threshold) + "."
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(FH_in, record, idx)
                    nb_valid_DP_in_others = 0
                    curr_allele.info["SUP"] = 0
                    curr_allele.info["CONT"] = 0
                    for curr_spl in FH_in.samples:
                        if curr_spl not in spl_in_gp and curr_allele.getDP(curr_spl) > args.others_DP_threshold:
                            nb_valid_DP_in_others += 1
                        if curr_allele.getAltAD(curr_spl)[0] >= args.AD_threshold and curr_allele.getAltAF(curr_spl)[0] >= args.AF_threshold:
                            if curr_spl in spl_in_gp:
                                curr_allele.info["SUP"] += 1
                        if curr_allele.getAltAD(curr_spl)[0] >= args.AD_threshold or curr_allele.getAltAF(curr_spl)[0] >= args.AF_threshold:
                            if curr_spl not in spl_in_gp:
                                curr_allele.info["CONT"] += 1
                    # Evaluates specificity of variant
                    if nb_valid_DP_in_others >= args.others_valid_DP and curr_allele.info["CONT"] <= args.max_contradict:
                        if curr_allele.info["SUP"]/float(len(spl_in_gp)) > args.nb_spl_ratio_threshold:
                            FH_out.write(curr_allele)
                            print(record.chrom, record.pos)
                            print("\t", "support", curr_allele.info["SUP"])
                            print("\t", "contradict", curr_allele.info["CONT"])
                            print("\t", "other valid DP", nb_valid_DP_in_others)
