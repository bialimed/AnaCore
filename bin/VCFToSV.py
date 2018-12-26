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
__version__ = '1.5.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.annotVcf import VEPVCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAlleleCounts(FH_vcf, record, allele_separator=","):
    """
    @summary: Returns the AD/DP for each allele for each sample.
    @param FH_vcf: [VCFIO] The VCF file object.
    @param record: [VCFRecord] The variant record.
    @param allele_separator: [str] The character used to separate the alleles information.
    @return: [list] Each element in the list represents one sample. In each sample the information for each allele is separated by "," (example: ["10/280,5/8", "800/812,62/62"]).
    """
    count_by_spl = list()
    for curr_spl in FH_vcf.samples:
        if curr_spl not in record.samples:
            count_by_spl.append("")
        else:
            AD = record.getAD(curr_spl)
            DP = record.getDP(curr_spl)
            alleles_count = list()
            for idx_allele, curr_AD in enumerate(AD):
                alleles_count.append("{}/{}".format(curr_AD, DP))
            count_by_spl.append(allele_separator.join(alleles_count))
    return count_by_spl


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Converts VCF annotated with VEP in separated value format (CSV, TSV, ...). One line in output file represents an annotation.')
    parser.add_argument('-s', '--separator', default='\t', help="Field separator in output file. [Default: tab]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file containing variants annotated with VEP v88+ (format: VCF).')
    args = parser.parse_args()

    # Process
    with VEPVCFIO(args.input_variants) as FH_vcf:
        print(
            "Chromosome",
            "Position",
            "Reference_allele",
            "Alternative_alleles",
            "Variant_quality",
            "Filters",
            "Alt_alleles_frequencies",
            args.separator.join(FH_vcf.samples),
            args.separator.join(FH_vcf.ANN_titles),
            sep=args.separator
        )
        for record in FH_vcf:
            record_fields = [
                record.chrom,
                str(record.pos),
                record.ref,
                ",".join(record.alt),  # Alternatives alleles
                ("." if record.qual is None else str(record.qual)),  # Quality
                ("." if record.filter is None else ",".join(record.filter)),  # Filters
                ",".join([str(round(AF, 5)) for AF in record.getPopAF()]),  # Alternatives AF
            ]
            record_fields.extend(
                getAlleleCounts(FH_vcf, record)  # Samples AF/DP
            )
            if len(record.info["CSQ"]) == 0:  # Record without consequence
                print(
                    args.separator.join(record_fields),
                    args.separator.join(["" for col in FH_vcf.ANN_titles]),
                    sep=args.separator
                )
            else:  # Record with at least one consequence
                for idx_csq, csq in enumerate(record.info["CSQ"]):  # For each consequence
                    # Pre-process record display
                    start_fields = ["" for col in record_fields]
                    if idx_csq == 0:
                        start_fields = record_fields
                    # Pre-process consequence display
                    csq_values = list()
                    for title in FH_vcf.ANN_titles:
                        if title not in csq or csq[title] is None:
                            csq_values.append("")
                        else:
                            csq_values.append(str(csq[title]))
                    # Display
                    print(
                        args.separator.join(start_fields),
                        args.separator.join(csq_values),
                        sep=args.separator
                    )
