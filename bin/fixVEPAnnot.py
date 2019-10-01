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

from anacore.annotVcf import AnnotVCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getVEPAlt(ref, alt):
    """
    Return the alternative allele in same format as annotation allele in VEP.

    :param ref: The reference allele.
    :type ref: str
    :param alt: The alternative allele.
    :type alt: str
    :return: The alternative allele in same format as annotation allele in VEP.
    :rtype: str
    """
    alleles = [ref] + alt
    # Replace empty marker by empty string
    for idx, cur_allele in enumerate(alleles):
        if cur_allele == "-":
            alleles[idx] = ""
    # Select shorter allele
    shorter_allele = alleles[0]
    for current_alt in alleles[1:]:
        if len(current_alt) < len(shorter_allele):
            shorter_allele = current_alt
    # Trim alleles
    trim = True
    while len(shorter_allele) != 0 and shorter_allele != "" and trim:
        for cur_allele in alleles:
            if len(cur_allele) == 0:
                trim = False
            elif cur_allele[0] != shorter_allele[0]:
                trim = False
        if trim:
            shorter_allele = shorter_allele[1:]
            for idx, cur_allele in enumerate(alleles):
                alleles[idx] = cur_allele[1:]
    # Replace empty by empty_marker
    for idx, cur_allele in enumerate(alleles):
        if cur_allele == "":
            alleles[idx] = "-"
    return alleles[1:]


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Reverse normalisation produced by VEP in allele annotation field.')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the output file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotations_field) as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotations_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord(FH_in, record, alt_idx)
                    vep_alt = getVEPAlt(alt_record.ref, alt_record.alt)[0]
                    for idx_ann, annot in enumerate(alt_record.info[FH_in.annot_field]):
                        if annot["Allele"] == vep_alt:
                            annot["Allele"] = alt_record.alt[0]
                FH_out.write(record)
    log.info("End of job")
