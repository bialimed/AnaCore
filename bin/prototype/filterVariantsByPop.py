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
    parser.add_argument('-g', '--group', nargs='+', help='*************************.')
    parser.add_argument('-gn', '--group-name', help='*************************.')
    group_thresholds = parser.add_argument_group('Group thresolds') # Group thresholds
    group_thresholds.add_argument('-ad', '--AD-threshold', type=int, default=3, help='*************************.')
    group_thresholds.add_argument('-dp', '--DP-threshold', type=int, default=10, help='*************************.')
    group_thresholds.add_argument('-af', '--AF-threshold', type=float, default=0.05, help='*************************.')
    group_thresholds.add_argument('-sr', '--support-ratio-threshold', type=float, default=round(1/float(3), 2), help='*************************.')
    group_thresholds.add_argument('-paf', '--popAF-threshold', type=float, default=0.01, help='*************************.')
    group_thresholds.add_argument('-sb', '--strand-ratio-threshold', type=float, default=0.5, help=' min(SAF, SAR)/max(SAF, SAR)*************************.')
    group_thresholds.add_argument('-q', '--variant-quality-threshold', type=float, default=1, help=' *************************.')
    out_thresholds = parser.add_argument_group('Out group thresolds') # Out group thresholds
    out_thresholds.add_argument('-oad', '--others-AD-threshold', type=int, default=3, help='*************************.')
    out_thresholds.add_argument('-odp', '--others-DP-threshold', type=int, default=10, help='*************************.')
    out_thresholds.add_argument('-oaf', '--others-AF-threshold', type=float, default=0.02, help='*************************.')
    out_thresholds.add_argument('-ondp', '--others-number-valid-DP', type=int, default=2, help='*************************.')
    out_thresholds.add_argument('-mc', '--max-contradict', type=int, default=0, help='*************************.')
    group_input = parser.add_argument_group('Inputs') # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs') # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.info["SUP"] = {
                "type": int,
                "type_tag": "Integer",
                "number": 1,
                "number_tag": "1",
                "description": "The number of samples in population supporting the variant. A sample supports a variant if the variant has an AD >= {}, a DP >= {}, an AF >= {} and an strand ratio >= {}.".format(
                    args.AD_threshold,
                    args.DP_threshold,
                    args.AF_threshold,
                    args.strand_ratio_threshold
                )
            }
            if args.group:
                FH_out.info["gpAF"] = {
                    "type": float,
                    "type_tag": "Float",
                    "number": 1,
                    "number_tag": "1",
                    "description": "The AF in group."
                }
                FH_out.info["CONT"] = {
                    "type": int,
                    "type_tag": "Integer",
                    "number": 1,
                    "number_tag": "1",
                    "description": "The number of samples out group supporting variant. A sample supports a variant if the variant has an AD >= {}, a DP >= {}, an AF >= {} and an strand ratio >= {}.".format(
                        args.others_AD_threshold,
                        args.others_DP_threshold,
                        args.others_AF_threshold,
                        args.strand_ratio_threshold
                    )
                }
                FH_out.filter["SPEC"] = "The variant is not specific at the group {}. The variant is specific if: (1) the ratio number_of_support_samples_in_group/number_of_samples_in_group is >= {} ; (2) the number of support samples out of the group is <= {} ; (3) at least {} samples out of the group have DP >= {}.".format(
                    ("" if args.group_name is None else " " + args.group_name),
                    args.support_ratio_threshold,
                    args.max_contradict,
                    args.others_number_valid_DP,
                    args.others_DP_threshold
                )
            FH_out._writeHeader()
            # Records
            nb_var = 0
            nb_valid = 0
            group = args.group if args.group is not None else FH_in.samples
            for record in FH_in:
                for idx in range(len(record.alt)):
                    nb_var += 1
                    if record.qual >= args.variant_quality_threshold:
                        curr_allele = getAlleleRecord(FH_in, record, idx)
                        # Count support and contradict
                        curr_allele.info["SUP"] = 0
                        group_AD = 0
                        group_DP = 0
                        nb_others_support = 0
                        nb_others_valid_DP = 0
                        test0 = []
                        test1 = []
                        test2 = []
                        for curr_spl in FH_in.samples:
                            strand_ratio = 1
                            if "SAF" in curr_allele.format and "SAR" in curr_allele.format:
                                SAF = [int(elt) for elt in curr_allele.samples[curr_spl]["SAF"]]
                                SAR = [int(elt) for elt in curr_allele.samples[curr_spl]["SAR"]]
                                strand_ratio = 0
                                if max(SAF[idx], SAR[idx]) != 0:
                                    strand_ratio = min(SAF[idx], SAR[idx])/max(SAF[idx], SAR[idx])
                            if curr_spl in group:
                                curr_AD = curr_allele.getAD(curr_spl)[0]
                                curr_DP = curr_allele.getDP(curr_spl)
                                group_AD += curr_AD
                                group_DP += curr_DP
                                if strand_ratio >= args.strand_ratio_threshold:
                                    if curr_allele.getAD(curr_spl)[0] >= args.AD_threshold and curr_allele.getAF(curr_spl)[0] >= args.AF_threshold and curr_allele.getDP(curr_spl) >= args.DP_threshold:
                                        test1.append(strand_ratio)
                                        test0.append(curr_allele.getAD(curr_spl)[0])
                                        curr_allele.info["SUP"] += 1
                            else:
                                #~ if strand_ratio >= args.strand_ratio_threshold:
                                if curr_allele.getDP(curr_spl) >= args.others_DP_threshold:
                                    nb_others_valid_DP += 1
                                    if curr_allele.getAD(curr_spl)[0] >= args.others_AD_threshold and curr_allele.getAF(curr_spl)[0] >= args.others_AF_threshold:
                                        test2.append(strand_ratio)
                                        nb_others_support += 1
                        # Evaluates specificity of variant
                        if curr_allele.info["SUP"]/float(len(group)) >= args.support_ratio_threshold: # The number of supporting samples in group is sufficient
                            if group_AD/float(group_DP) >= args.popAF_threshold: # The number of supporting samples in group is sufficient
                                    if args.group:
                                        curr_allele.info["gpAF"] = group_AD/float(group_DP)
                                        curr_allele.info["CONT"] = nb_others_support
                                        if nb_others_valid_DP >= args.others_number_valid_DP and nb_others_support <= args.max_contradict:
                                            print(
                                                curr_allele.ref, "/", curr_allele.alt[0],
                                                curr_allele.info["SUP"], "=>", curr_allele.info["SUP"]/float(len(group)) >= args.support_ratio_threshold, " ; ",
                                                group_AD/float(group_DP), "=>", group_AD/float(group_DP) >= args.popAF_threshold, " ; ",
                                                nb_others_support, "=>", nb_others_support <= args.max_contradict, " ; ",
                                                nb_others_valid_DP, "=>", nb_others_valid_DP >= args.others_number_valid_DP, " ; ",
                                                test0,
                                                test1,
                                                test2
                                            )
                                            nb_valid += 1
                                            FH_out.write(curr_allele)
                                    else:
                                        nb_valid += 1
                                        FH_out.write(curr_allele)
            print(nb_valid, nb_var)
