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

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastaIO
from anacore.vcf import VCFIO, VCFRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getChromSeq(chrom_name, in_fasta):
    """
    Return the sequence corresponding to the chromosome.

    :param chrom_name: The name of the selected chromosome.
    :type chrom_name: str
    :param in_fasta: The path to the file sequences file (format: fasta).
    :type in_fasta: str
    :return: the sequence corresponding to the chromosome.
    :rtype: str
    """
    seq = None
    with FastaIO(in_fasta) as FH_ref:
        for record in FH_ref:
            if record.id == chrom_name:
                seq = record.string
    if seq is None:
        raise Exception('The chromosome "{}" cannot be rertrieved from "{}".'.format(chrom_name, in_fasta))
    return seq


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Replaces in/del marker by the syntax described in the specification 4.1 of VCF. Example: chr20<tab>1234568<tab>.<tab>TC<tab>-,TCT becomes chr20<tab>1234567<tab>.<tab>GTC<tab>G,GTCT.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-ir', '--input-reference', required=True, help='The path to the reference used for the variant calling (format: fasta).')
    group_input.add_argument('-iv', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).')
    args = parser.parse_args()

    # Process
    curr_chrom = {"name": "", "seq": None}
    with VCFIO(args.output_variants, "w") as FH_out_vcf:
        with VCFIO(args.input_variants) as FH_in_vcf:
            # Header
            FH_out_vcf.copyHeader(FH_in_vcf)
            FH_out_vcf._writeHeader()
            # Records
            for record in FH_in_vcf:
                if record.ref == VCFRecord.getEmptyAlleleMarker() or any([alt == VCFRecord.getEmptyAlleleMarker() for alt in record.alt]):  # record is a standardized in/del
                    # Get previous nt
                    if record.chrom != curr_chrom["name"]:
                        curr_chrom["name"] = record.chrom
                        curr_chrom["seq"] = getChromSeq(record.chrom, args.input_reference)
                    prev_nt = curr_chrom["seq"][record.pos - 2]
                    # Update record
                    record.pos -= 1
                    if record.ref == VCFRecord.getEmptyAlleleMarker():  # Insertion
                        record.ref = prev_nt
                    else:  # Deletion
                        record.ref = prev_nt + record.ref
                    for idx, alt in enumerate(record.alt):
                        alt = alt.replace(VCFRecord.getEmptyAlleleMarker(), "")  # One of alt can be a deletion
                        record.alt[idx] = prev_nt + alt
                FH_out_vcf.write(record)
