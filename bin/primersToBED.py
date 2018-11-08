#!/usr/bin/env python3
#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import os
import sys
import warnings
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.bed import BEDIO, BEDRecord
from anacore.sequenceIO import FastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def revcom(seq):
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}

    return("".join([complement_rules[base] for base in seq[::-1]]))

def getAmplicons(amplicons_path):
    amplicons = []
    with open(amplicons_path) as FH_ampl:
        for line in FH_ampl:
            if not line.startswith("#"):
                fields = [elt.strip() for elt in line.split("\t")]
                amplicons.append({
                    "name": fields[0].replace(" ", "_"),
                    "f_primer": fields[1].upper(),
                    "r_primer": revcom(fields[2].upper())
                })
    return(amplicons)

def getBEDRecords(ref_path, amplicons):
    for ampl in amplicons:
        ampl["found"] = False
    bed_ampl = []
    with FastaIO(ref_path) as FH_seq:
        for record in FH_seq:
            chr_id = record.id
            chr_str = record.string.upper()
            for ampli in amplicons:
                # Primers are on strand +
                up_primer = ampli["f_primer"].upper()
                down_primer = ampli["r_primer"].upper()
                start, end = findPosOnSequence(chr_id, chr_str, up_primer, down_primer)
                if start is not None:
                    ampli["found"] = True
                    bed_ampl.append(
                        BEDRecord(
                            chr_id,
                            start,
                            end,
                            ampli["name"],
                            0,
                            "+",
                            start + len(up_primer),
                            end - len(down_primer)
                        )
                    )
                # Primers are on strand -
                up_primer = revcom(ampli["r_primer"].upper())
                down_primer = revcom(ampli["f_primer"].upper())
                start, end = findPosOnSequence(chr_id, chr_str, up_primer, down_primer)
                if start is not None:
                    ampli["found"] = True
                    bed_ampl.append(
                        BEDRecord(
                            chr_id,
                            start,
                            end,
                            ampli["name"],
                            0,
                            "-",
                            start + len(up_primer),
                            end - len(down_primer)
                        )
                    )
    for ampl in amplicons:
        if not ampl["found"]:
            warnings.warn('The amplicons {} with primers fwd:{}, rvs:{} cannot be found in {}.'.format(
                ampl["name"], ampl["f_primer"], ampl["r_primer"], ref_path
            ))
    return(bed_ampl)

def findPosOnSequence(chr_id, chr_str, up_primer, down_primer):
    ampl_start = None
    ampl_end = None
    upstream_matches = list()
    up_pattern = re.compile(up_primer)
    for curr_match in up_pattern.finditer(chr_str):
        upstream_matches.append({"start": curr_match.start()+1, "end": curr_match.end()})
    downstream_matches = list()
    down_pattern = re.compile(down_primer)
    for curr_match in down_pattern.finditer(chr_str):
        downstream_matches.append({"start": curr_match.start()+1, "end": curr_match.end()})
    if len(upstream_matches) != 0 and len(downstream_matches) != 0:
        # Check multiple target in chr
        if len(upstream_matches) > 1:
            match_list = ", ".join(["{}:{}-{}".format(chr_id, curr_match["start"], curr_match["end"]) for curr_match in upstream_matches])
            warnings.warn("The primer '" + up_primer + "' is found multiple twice in " + chr_id + " (" + match_list + ")")
        if len(downstream_matches) > 1:
            match_list = ", ".join(["{}:{}-{}".format(chr_id, curr_match["start"], curr_match["end"]) for curr_match in downstream_matches])
            warnings.warn("The primer '" + down_primer + "' is found multiple twice in " + chr_id + " (" + match_list + ")")
        # Select smaller amplified fragment
        prev_length = None
        for curr_up in upstream_matches:
            for curr_down in downstream_matches:
                curr_length = curr_down["start"] - curr_up["end"]
                if curr_length >= 0:
                    if prev_length is None or prev_length > curr_length:
                        prev_length = curr_length
                        ampl_start = curr_up["start"]
                        ampl_end = curr_down["end"]
    return(ampl_start, ampl_end)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Converts amplicons list to BED.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-amplicons', required=True, help='Path to the definition of the amplicons (format: TSV). The columns are: Name, Forward_primer and Reverse_primer.')
    group_input.add_argument('-g', '--input-reference', required=True, help='Path to the reference sequence where the amplicon will be searched (format: fasta).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-BED', default="amplicons.bed", help='The amplicons description (format: BED). [Default: %(default)s]')
    args = parser.parse_args()

    amplicons = getAmplicons(args.input_amplicons)
    bed_records = getBEDRecords(args.input_reference, amplicons)
    with BEDIO(args.output_BED, "w", 8) as FH_out:
        for record in bed_records:
            FH_out.write(record)
