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

from anacore.sv import HashedSVIO
from anacore.annotVcf import AnnotVCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Add in annotations the pathways related to the gene affected by the variant.")
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-g', '--gene-field', default="SYMBOL", help='The field used in annotation to store the gene HUGO symbol. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-p', '--input-pathways', required=True, help='Path to the pathways file (format: GMT). The pathway ID is strored in the 2nd column and the genes list start at the 4th column (see https://reactome.org/download/current/ReactomePathways.gmt.zip).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))
    log.info("Version: " + str(__version__))

    # Get pathways
    pathways_by_gene = {}
    with open(args.input_pathways) as FH_pathways:
        for line in FH_pathways:
            fields = [elt.strip() for elt in line.split("\t")]
            for gene in fields[3:]:
                if gene not in pathways_by_gene:
                    pathways_by_gene[gene] = set()
                pathways_by_gene[gene].add(fields[1])

    # Write output
    with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotation_field) as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotation_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.ANN_titles.append("Pathways")
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                for annot in record.info[FH_in.annot_field]:
                    if annot[args.gene_field] is not None and annot[args.gene_field] != "":
                        pathways = set()
                        for gene in annot[args.gene_field].split(","):
                            if gene in pathways_by_gene:
                                pathways = pathways | pathways_by_gene[gene]
                    annot["Pathways"] = "&".join([curr_pathway.replace(",", "%2C") for curr_pathway in sorted(pathways)])
                FH_out.write(record)

    log.info("End of job")
