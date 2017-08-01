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
__email__ = 'escudie.frederiic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from GTFI import *


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Write the number of reads by biotype from HTSeq-count gene output.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-g', '--input-gtf', required=True, help='The ensembl GTF used in HTSeq-count (format: GTF).' )
    group_input.add_argument( '-c', '--input-count', required=True, help='The HTSeq-count output (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='The path for the outputed file (format: TSV).')
    args = parser.parse_args()

    # Get biotype by gene ID
    biotype_by_id = dict()
    FH_gff = GTFI(args.input_gtf)
    FH_gff.open()
    try:
        for record in FH_gff:
            gene_id = record["attr"]["gene_id"]
            gene_biotype = "Unknown"
            if "gene_biotype" in record["attr"]:
                gene_biotype = record["attr"]["gene_biotype"]
            if gene_id in biotype_by_id and biotype_by_id[gene_id] != gene_biotype:
                raise Exception("The gene with ID '" + gene_id + "' has already be described in '" + args.input_gtf + "'.")
            biotype_by_id[gene_id] = gene_biotype
    finally:
        FH_gff.close()

    # Get count by biotype 
    ################################################# Pb gene de diff longueur donne diff de count
    count_by_biotype = dict()
    samples = list()
    count_by_spl = dict()
    with open(args.input_count) as FH_count:
        samples = [field.strip() for field in FH_count.readline().split("\t")[1:]]
        count_by_spl = {spl:0 for spl in samples}
        for line in FH_count:
            ref = line.split("\t", 1)[0].strip()
            counts = [int(count.strip()) for count in line.split("\t")[1:]]
            biotype = biotype_by_id[ref]
            if biotype not in count_by_biotype:
                count_by_biotype[biotype] = [0 for spl in samples]
            for idx_spl, count in enumerate(counts):
                count_by_biotype[biotype][idx_spl] += count
                count_by_spl[samples[idx_spl]] += count

    # Write output
    with open(args.output_file, "w") as FH_out:
        FH_out.write( "#Biotype\t" + "\t".join(samples) + "\n" )
        for biotype in count_by_biotype:
            biotype_prct = list()
            for idx_spl, spl in enumerate(samples):
                biotype_count = count_by_biotype[biotype][idx_spl]
                prct = float(biotype_count * 100)/count_by_spl[spl]
                biotype_prct.append( prct )
            FH_out.write( biotype + "\t" + "\t".join(map(str, biotype_prct)) + "\n" )
