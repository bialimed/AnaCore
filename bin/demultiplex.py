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
import re
import sys
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastaIO, FastqIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getBarcodes(in_barcodes):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons information.
    """
    barcodes = list()
    with open(in_barcodes) as FH_in:
        for line in FH_in:
            if not line.startswith("#"):
                barcode_id, barcode_fwd, barcode_rvs = [elt.strip() for elt in line.split("\t")]
                barcodes.append({
                    "id": barcode_id,
                    "fwd": barcode_fwd,
                    "rvs": barcode_rvs
                })
    return barcodes

def get_seq_ids(fastq_path):
    ids = list()
    with FastqIO(fastq_path) as FH:
        for record in FH:
            ids.append(record.id)
    return ids

def pickSequences(in_path, out_path, retained_ids):
    dict_retained_ids = {curr_id: 1 for curr_id in retained_ids}
    with FastqIO(out_path, "w") as FH_out:
        with FastqIO(in_path) as FH_in:
            for record in FH_in:
                if record.id in dict_retained_ids:
                    FH_out.write(record)

def cutadapt(in_fastq, out_fastq, adapter_seq, error_rate=0.1):
    cmd = [
        "cutadapt",
        "--error-rate", str(error_rate),
        "-g", adapter_seq,
        "--discard-untrimmed",
        "-o", out_fastq,
        in_fastq
    ]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL)

def getLibNameFromReadPath(fastq_path):
    library_name = os.path.basename(fastq_path).split(".")[0]
    if re.search('_[rR][1-2]$', library_name):
        library_name = library_name[:-3]
    elif re.search('_[rR][1-2]_\d\d\d$', library_name):
        library_name = library_name[:-7]
    return library_name


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # design avec un des 2 primers proche d'un autre a 10% près c'est le premier dans la liste qui va prendre le pas ? Dans ce cas faire de la recup par couple de prime ?
    # au moins un primer si non ambigu ? ############################

    # Manage parameters
    parser = argparse.ArgumentParser(description='**************************************.')
    parser.add_argument('-e', '--error-rate', default=0.1, type=float, help='************************************. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs') # Inputs
    group_input.add_argument('-b', '--barcodes', required=True, help='************ (format: TSV).')
    group_input.add_argument('-R1', '--R1-path', required=True, help='************ (format: fastq).')
    group_input.add_argument('-R2', '--R2-path', required=True, help='************ (format: fastq).')
    group_output = parser.add_argument_group('Outputs') # Outputs
    group_output.add_argument('-o', '--output-filename-pattern', help='**********************.')
    group_output.add_argument('-d', '--output-dir', default=os.getcwd(), help='**********************. [Default: %(default)s]')
    args = parser.parse_args()

    # Set output variables
    if args.output_filename_pattern is None:
        args.output_filename_pattern = getLibNameFromReadPath(args.R1_path) + "_<AMPLI>_<R>.fastq.gz"

    # Load amplicons param
    barcodes = getBarcodes(args.barcodes)

    # Demultiplex
    ids_by_barcodes = dict()
    for curr_barcode in barcodes:
        ampl_out_R1 = os.path.join(args.output_dir, args.output_filename_pattern.replace("<R>", "R1").replace("<AMPLI>", curr_barcode["id"]))
        ampl_out_R2 = os.path.join(args.output_dir, args.output_filename_pattern.replace("<R>", "R2").replace("<AMPLI>", curr_barcode["id"]))

        # Find adapter
        cutadapt_out_R1 = os.path.join(args.output_dir, "tmp_" + args.output_filename_pattern.replace("<R>", "R1").replace("<AMPLI>", curr_barcode["id"]))
        cutadapt_out_R2 = os.path.join(args.output_dir, "tmp_" + args.output_filename_pattern.replace("<R>", "R2").replace("<AMPLI>", curr_barcode["id"]))
        cutadapt(args.R1_path, cutadapt_out_R1, "^" + curr_barcode["fwd"], args.error_rate)
        cutadapt(args.R2_path, cutadapt_out_R2, "^" + curr_barcode["rvs"], args.error_rate)

        # Select reads for the amplicon
        R1 = set(get_seq_ids(cutadapt_out_R1))
        R2 = set(get_seq_ids(cutadapt_out_R2))
        retained_ids = R1.intersection(R2)
        ids_by_barcodes[curr_barcode["id"]] = retained_ids

        # Filter
        pickSequences(args.R1_path, ampl_out_R1, retained_ids)
        pickSequences(args.R2_path, ampl_out_R2, retained_ids)

        # Clean
        os.remove(cutadapt_out_R1)
        os.remove(cutadapt_out_R2)

        print(curr_barcode["id"], len(retained_ids), curr_barcode["fwd"], curr_barcode["rvs"], sep="\t")

    # Check ambiguous
    ids = dict()
    non_uniq = dict()
    for barcode in ids_by_barcodes:
        for seq_id in ids_by_barcodes[barcode]:
            if seq_id in ids:
                if seq_id not in non_uniq:
                    non_uniq[seq_id] = 1
                non_uniq[seq_id] += 1
            ids[seq_id] = 1
    print("Ambiguous", len(non_uniq), sep="\t")
