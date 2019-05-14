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
__version__ = '1.3.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.annotVcf import AnnotVCFIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGeneByNM(gene_to_id_file, trim_version=False):
    """
    Return gene name by RNA_id.

    :param gene_to_id_file: Path to the file describing the link between genes and RNA_id (format: TSV). Each line has the following format: <GENE>\t<RNA_ID>.
    :type gene_to_id_file: str
    :param trim_version: With True the version number is removed from the id.
    :type trim_version: bool
    :return: The genes by RNA_id.
    :rtype: dict
    """
    gene_by_NM = dict()
    with open(gene_to_id_file) as FH_ref:
        for line in FH_ref:
            if not line.startswith("#"):
                gene, NM = [field.strip() for field in line.split("\t")]
                if trim_version:
                    NM = NM.split(".")[0]
                gene_by_NM[NM] = gene
    return gene_by_NM


def filterRecordAnnot(record, kept_id, trim_version=False, annotation_field="ANN"):
    """
    Remove annotations that does not come from a kept id.

    :param record: The annotated record.
    :type record: anacore.vcf.VCFRecord
    :param kept_id: The keys of this dictionary are the kept id.
    :type kept_id: dict
    :param trim_version: With True the version number is removed from the id.
    :type trim_version: bool
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    """
    removed_annot_idx = list()
    for annot_idx, annot in enumerate(record.info[annotation_field]):
        RNA_id = annot["Feature"]
        if RNA_id is not None and trim_version:
            RNA_id = RNA_id.split(".")[0]
        if RNA_id not in kept_id:
            removed_annot_idx.append(annot_idx)
    for curr_idx in sorted(removed_annot_idx, reverse=True):
        del(record.info[annotation_field][curr_idx])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Removes annotations that does not come from a list of RNA. The RNA source of one annotation is retrieved from the field "Feature".')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-w', '--without-version', action='store_true', help='With this option the version number of the NM is not used in filter.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-r', '--reference-RNA', required=True, help='The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.')
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Get kept RNA_ID
    kept_ID = getGeneByNM(args.reference_RNA, args.without_version)

    # Filter annotations
    with AnnotVCFIO(args.input_variants, "r", args.annotation_field) as FH_in:
        with AnnotVCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                filterRecordAnnot(record, kept_ID, args.without_version, FH_in.annot_field)
                FH_out.write(record)
