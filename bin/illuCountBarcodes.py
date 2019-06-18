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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import copy
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastqIO
from anacore.illumina import getInfFromSeqDesc


########################################################################
#
# FUNCTIONS
#
########################################################################
def getCountByBarcode(in_seq):
    """
    Return the number of reads by barcode in the fastq file.

    :param in_seq: The path to the sequence file (format: fastq).
    :type in_seq: int
    :return: The number of reads by barcode.
    :rtype: dict
    """
    count_by_barcode = dict()
    for curr_seq in in_seq:
        with FastqIO(curr_seq) as FH_in:
            for record in FH_in:
                barcode = getInfFromSeqDesc(record.description)["barcode"]
                if barcode not in count_by_barcode:
                    count_by_barcode[barcode] = 1
                else:
                    count_by_barcode[barcode] += 1
    return count_by_barcode

def getTopN(count_by_key, nb_top, is_desc=True):
    """
    Return the number of reads by barcode for the N most or least abundant.

    :param count_by_key: The number of reads by barcode.
    :type count_by_key: dict
    :param nb_top: The number of selected barcodes.
    :type nb_top: int
    :param is_desc: If true the most abudant are selected. otherwise, the least abundant are selected.
    :type is_desc: bool
    :return: The number of reads by barcode for the N most or least abundant.
    :rtype: dict
    """
    restricted_count_by_key = dict()
    if len(count_by_key) <= nb_top:
        restricted_count_by_key = copy.deepcopy(count_by_key)
    else:
        restricted_count_by_key = dict()
        for key in sorted(count_by_key, key=count_by_key.get, reverse=is_desc)[:nb_top]:
            restricted_count_by_key[key] = count_by_key[key]
    return restricted_count_by_key

def writeJSON(data, out_path):
    """
    Write the results in a JSON file.

    :param data: Analysis results.
    :type data: dict
    :param out_path: Path to the output file.
    :type out_path: str
    """
    with open(out_path, "w") as FH_out:
        FH_out.write(
            json.dumps(data, default=lambda o: o.__dict__, sort_keys=True)
        )

def writeSV(data, out_path, separator="\t"):
    """
    Write the results in a JSON file.

    :param data: Analysis results.
    :type data: dict
    :param out_path: Path to the output file.
    :type out_path: str
    :param separator: The field separator used in the output file.
    :type separator: str
    """
    with open(out_path, "w") as FH_out:
        FH_out.write(
            separator.join(["#Barcode", "Nb_barcodes", "Count\n"])
        )
        for key in sorted(data["count_by_barcode"], key=data["count_by_barcode"].get, reverse=True):
            FH_out.write(
                separator.join([key, "1", str(data["count_by_barcode"][key]) + "\n"])
            )
        # Add count for others barcodes
        if data["limit_barcodes"] is not None and data["limit_barcodes"] < data["nb_barcodes"]:
            restricted_nb_barcodes = data["limit_barcodes"]
            restricted_total_count = sum([count for barcode, count in data["count_by_barcode"].items()])
            FH_out.write(
                separator.join([
                    "others_barcodes",
                    str(data["nb_barcodes"] - restricted_nb_barcodes),
                    str(data["total_count"] - restricted_total_count) + "\n"
                ])
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Count the number o sequences by barcode in fastq file(s).")
    parser.add_argument('-m', '--nb-max', type=int, help='The maximum number of barcodes to report. With this option only the top N of most represented barcodes are detailed on output the others are merged inone group.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-seq', required=True, nargs='+', help='The sequence files produced by an Illumina sequencer (format: fastq).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-oj', '--output-json', help='Path to the output file (format: JSON).')
    group_output.add_argument('-ot', '--output-tsv', help='Path to the output file (format: TSV).')
    args = parser.parse_args()

    if args.output_json is None and args.output_tsv is None:
        parser.error('"--output-json" and/or "--output-tsv" must be specified.')

    # Process
    count_by_barcode = getCountByBarcode(args.inputs_seq)
    data = {
        "count_by_barcode": count_by_barcode,
        "total_count": sum([count for barcode, count in count_by_barcode.items()]),
        "nb_barcodes": len(count_by_barcode),
        "limit_barcodes": args.nb_max
    }
    if args.nb_max is not None:
        data["count_by_barcode"] = getTopN(count_by_barcode, args.nb_max)
    if args.output_json is not None:
        writeJSON(data, args.output_json)
    if args.output_tsv is not None:
        writeSV(data, args.output_tsv)
