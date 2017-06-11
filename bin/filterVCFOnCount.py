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
__version__ = '1.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFIO, getAlleleRecord



########################################################################
#
# FUNCTIONS
#
########################################################################
def getFilterTag( DP_is_ok, AF_is_ok ):
    tag = "PASS"
    if not DP_is_ok[0] and not DP_is_ok[1]: # -- **
        tag = "lowDP"
    elif DP_is_ok[0] and DP_is_ok[1]: # ++ **
        if not AF_is_ok[0] and not AF_is_ok[1]: # ++ --
            tag = "lowAF"
        elif not AF_is_ok[0] or not AF_is_ok[1]: # ++ +-/-+
            tag = "libSpe"
    else: # +-/-+ **
        if AF_is_ok[0] and AF_is_ok[1]: # +-/-+ ++
            tag = "incomplete"
        elif DP_is_ok[0] and not DP_is_ok[1] and AF_is_ok[0] and not AF_is_ok[1]: # +- +-
            tag = "incomplete"
        elif not DP_is_ok[0] and DP_is_ok[1] and not AF_is_ok[0] and AF_is_ok[1]: # -+ -+
            tag = "incomplete"
        else: # +- -+ / -+ +-
            tag = "invalid"
    return tag


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='************************************************' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]' )
    group_filter = parser.add_argument_group( 'Filters' ) # Filters
    group_filter.add_argument( '-a', '--AF-threshold', type=float, default=0.02, help='*****************************************.' )
    group_filter.add_argument( '-d', '--DP-threshold', type=int, default=120, help='*****************************************.' )    
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    with open(args.output_variants, "w") as FH_out:
        # Filter variant file
        with open(args.input_variants) as FH_vcf:
            line = FH_vcf.readline()
            while line.startswith("#"):
                ######################################### add filter info
                FH_out.write( line )
                line = FH_vcf.readline()
        # Find retained variants
        with VCFIO(args.input_variants) as FH_in:
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord( FH_in, record, alt_idx )
                    # Evaluates filters
                    AF_is_ok = list()
                    DP_is_ok = list()
                    for spl_idx, curr_spl in enumerate(alt_record.samples):
                        AF_is_ok.append( True )
                        if alt_record.get_AF(curr_spl)[0] < args.AF_threshold:
                            AF_is_ok[spl_idx] = False
                        DP_is_ok.append( True )
                        if alt_record.get_DP(curr_spl) < args.DP_threshold:
                            DP_is_ok[spl_idx] = False
                        print(curr_spl, alt_record.get_AF(curr_spl), alt_record.get_DP(curr_spl))
                    print(AF_is_ok)
                    print(DP_is_ok)
                    # Apply filters
                    tag = getFilterTag( DP_is_ok, AF_is_ok )
                    if tag == "PASS":
                        if len(alt_record.filter) == 0:
                            alt_record.filter.append( "PASS" )
                    else:
                        alt_record.filter.append( tag )
                    if args.mode == "tag":
                        FH_out.write( FH_in.recToVCFLine(alt_record) + "\n" )
                    else:
                        if tag in ["PASS", "incomplete"]:
                            FH_out.write( FH_in.recToVCFLine(alt_record) + "\n" )
