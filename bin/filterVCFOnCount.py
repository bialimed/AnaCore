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
from argparse import RawTextHelpFormatter

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getFilterTag(DP_is_ok, AF_is_ok):
    """
    @summary: Returns the tag corresponding to the filters results (DP and AF) for the variant.
    @param DP_is_ok: [list] Each element is a sample's result against the DP's filter. Example: [False, True] => sample 1 has an insufficient DP and sample 2 has a sufficient DP.
    @param AF_is_ok: [list] Each element is a sample's result against the AF's filter. Example: [True, False] => sample 1 has a sufficient AF and sample 2 has an insufficient AF.
    @return: [str] The filter tag.
    """
    tag = "PASS"
    if not DP_is_ok[0] and not DP_is_ok[1]:  # -- **
        tag = "lowDP"
    elif DP_is_ok[0] and DP_is_ok[1]:  # ++ **
        if not AF_is_ok[0] and not AF_is_ok[1]:  # ++ --
            tag = "lowAF"
        elif not AF_is_ok[0] or not AF_is_ok[1]:  # ++ +-/-+
            tag = "libSpe"
    else:  # +-/-+ **
        if AF_is_ok[0] and AF_is_ok[1]:  # +-/-+ ++
            tag = "incomplete"
        elif DP_is_ok[0] and not DP_is_ok[1] and AF_is_ok[0] and not AF_is_ok[1]:  # +- +-
            tag = "incomplete"
        elif not DP_is_ok[0] and DP_is_ok[1] and not AF_is_ok[0] and AF_is_ok[1]:  # -+ -+
            tag = "incomplete"
        else:  # +- -+ / -+ +-
            tag = "invalid"
    return tag


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='''Filters variants on AF and DP accross samples. The AF must be sufficient in two samples.
  Case   Sample_1   Sample_2   FILTER
         DP  AF     DP  AF
  1      +   +      +   +      PASS
  2      +   -      +   -      lowAF
  3      -   +      -   +      lowDP
  4      -   -      -   -      lowDP
  5      +   +      -   -      incomplete
  6      +   -      -   +      invalid
  7      +   +      -   +      incomplete
  8      +   -      -   -      invalid
  9      -   +      +   -      invalid
  10     -   -      +   +      incomplete
  11     -   +      +   +      incomplete
  12     -   -      +   -      invalid
  13     -   -      -   +      lowDP
  14     +   -      +   +      libSpe
  15     -   +      -   -      lowDP
  16     +   +      +   -      libSpe
  In "remove" mode only "PASS" and "incomplete" variants are kept.
''')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]')
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-a', '--AF-threshold', type=float, default=0.02, help='The minimum allele frequency to validate the variant. [Default: %(default)s]')
    group_filter.add_argument('-d', '--DP-threshold', type=int, default=120, help='The minimum depth in sample to validate the variant. [Default: %(default)s].')
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
            FH_out.filter["incomplete"] = 'The variant has a sufficient frequency in one sample (AF>=' + str(args.AF_threshold * 100) + '% and DP>=' + str(args.DP_threshold) + ') but the depth in second sample it has an insufficiant depth (the AF cannot be take into account).'
            FH_out.filter["invalid"] = 'The variant has an insufficient frequency in one sample (AF<' + str(args.AF_threshold * 100) + '% and DP>=' + str(args.DP_threshold) + ') but the depth in second sample it has an insufficiant depth (the AF cannot be take into account).'
            FH_out.filter["libSpe"] = 'The variant has a frequency lower than ' + str(args.AF_threshold * 100) + '% in one of the sample and depth is superior than ' + str(args.DP_threshold) + ' in all the samples.'
            FH_out.filter["lowAF"] = 'The variant has a frequency lower than ' + str(args.AF_threshold * 100) + '% in all the samples and depth is superior than ' + str(args.DP_threshold) + ' in all the samples.'
            FH_out.filter["lowDP"] = 'The variant has a depth lower than ' + str(args.DP_threshold) + ' in all the samples.'
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord(FH_in, record, alt_idx)
                    # Evaluates filters
                    AF_is_ok = list()
                    DP_is_ok = list()
                    for spl_idx, curr_spl in enumerate(alt_record.samples):
                        AF_is_ok.append(True)
                        if alt_record.getAltAF(curr_spl)[0] < args.AF_threshold:
                            AF_is_ok[spl_idx] = False
                        DP_is_ok.append(True)
                        if alt_record.getDP(curr_spl) < args.DP_threshold:
                            DP_is_ok[spl_idx] = False
                    # Apply filters
                    tag = getFilterTag(DP_is_ok, AF_is_ok)
                    if tag == "PASS":
                        if len(alt_record.filter) == 0:
                            alt_record.filter.append("PASS")
                    else:
                        if len(alt_record.filter) == 1 and alt_record.filter[0] == "PASS":
                            alt_record.filter[0] = tag
                        else:
                            alt_record.filter.append(tag)
                    if args.mode == "tag":
                        FH_out.write(alt_record)
                    else:
                        if tag in ["PASS", "incomplete"]:
                            FH_out.write(alt_record)
