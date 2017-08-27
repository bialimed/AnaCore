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
__version__ = '1.3.0'
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

from VEPvcf import VEPVCFIO, getAlleleRecord



########################################################################
#
# FUNCTIONS
#
########################################################################
def getAlleleCounts( FH_vcf, record, allele_separator="," ):
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
            AD = record.get_AD(curr_spl)
            DP = record.getDP(curr_spl)
            alleles_count = list()
            for idx_allele, curr_AD in enumerate(AD):
                alleles_count.append( str(curr_AD) + "/" + str(DP) )
            count_by_spl.append( allele_separator.join(alleles_count) )
    return( count_by_spl )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Converts VCF annotated with VEP in separated value format (CSV, TSV, ...). One line in output file represents an annotation.' )
    parser.add_argument( '-s', '--separator', default='\t', help="Field separator in output file. [Default: tab]")
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-variants', required=True, help='The path to the file containing variants annotated with VEP v88+ (format: VCF).' )
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
            args.separator.join( FH_vcf.samples ),
            args.separator.join( FH_vcf.CSQ_titles ),
            sep=args.separator
        )
        for record in FH_vcf:
            alt = ",".join(record.alt)
            qual = "." if record.qual is None else record.qual
            filters = "." if record.filter is None else ",".join(record.filter)
            alt_AF = ",".join([str(round(AF, 5)) for AF in record.getPopAF()])
            spl_counts = args.separator.join( getAlleleCounts(FH_vcf, record) )
            for idx_csq, csq in enumerate(record.info["CSQ"]):
                csq_values = list()
                for title in FH_vcf.CSQ_titles:
                    if title not in csq or csq[title] is None:
                        csq_values.append("")
                    else:
                        csq_values.append( str(csq[title]) )
                print(
                    (record.chrom if idx_csq == 0 else ""),
                    (record.pos if idx_csq == 0 else ""),
                    (record.ref if idx_csq == 0 else ""),
                    (alt if idx_csq == 0 else ""),
                    (qual if idx_csq == 0 else ""),
                    (filters if idx_csq == 0 else ""),
                    (alt_AF if idx_csq == 0 else ""),
                    (spl_counts if idx_csq == 0 else args.separator*(len(FH_vcf.samples) - 1)),
                    args.separator.join( csq_values ),
                    sep=args.separator
                )
