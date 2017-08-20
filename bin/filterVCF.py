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
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def getRecordValue( record, target ):
    target_tree = target.split(".")
    value = record
    for key in target_tree:
        if issubclass(value.__class__, dict):
            value = value[key]
        elif issubclass(value.__class__, object):
            value = value.__getattribute__(key)
    return value

def isKept( eval_fct, eval_value, curr_filter ):
    match_filter = eval_fct( record_value, curr_filter["values"] )
    is_ok = False
    if match_filter and curr_filter["action"] == "select":
        is_ok = True
    elif not match_filter and curr_filter["action"] == "exclude":
        is_ok = True
    return is_ok

def getStdEvalFct( operator ):
    """
    @summary: Returns a function used for evaluate a value against a reference. This function take two parameters and returns a boolean.
    @param operator: [str] The operator used in comparison.
    @return [fct] The evaluation function.
    """
    eval_fct = None
    if operator in ["=", "==", "eq"]:
        def eval_fct(val, ref): return( val == ref )
    elif operator in ["!=", "<>", "ne"]:
        def eval_fct(val, ref): return( val != ref )
    elif operator in ["<=", "le"]:
        def eval_fct(val, ref): return( val <= ref )
    elif operator in [">=", "ge"]:
        def eval_fct(val, ref): return( val >= ref )
    elif operator in ["<", "lt"]:
        def eval_fct(val, ref): return( val < ref )
    elif operator in [">", "gt"]:
        def eval_fct(val, ref): return( val > ref )
    elif operator == "in":
        def eval_fct(val, ref): return( val in ref )
    elif operator == "not in":
        def eval_fct(val, ref): return( val not in ref )
    elif operator == "at least one in":
        def eval_fct(val, ref):
            one_in = False
            if not issubclass(val.__class__, list) and not issubclass(val.__class__, tuple) and issubclass(val.__class__, set):
                raise Exception('"atLeastOneIn" must be applied on list.')
            for curr_elt in val:
                if curr_elt in ref:
                    one_in = True
            return( one_in )
    else:
        raise Exception("The operator '" + operator + "' is not implemented.")
    return eval_fct


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Sorts VCF by coordinates.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--filters', required=True, help='The path to the filters file (format: JSON).' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    filters = None
    with open(args.filters) as data_file:
        filters = json.load(data_file)
    eval_fct = getStdEvalFct( filters["operator"] )
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader( FH_in )
            FH_out._writeHeader()

            # Records
            for record in FH_in:
                record_value = getRecordValue( record, filters["target"] )
                if isKept(eval_fct, record_value, filters):
                    FH_out.write( record )
