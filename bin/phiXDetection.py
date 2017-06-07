#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT
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
__copyright__ = 'Copyright (C) 2017 IUCT'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='***********************.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--fields', default=['@observation_name', '@observation_sum', '@sample_count'], nargs='+', help="Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count' '@rdp_tax_and_bootstrap', '@seed_sequence' . The others columns must be metadata title. [Default: %(default)s]")
    parser.add_argument( '-s', '--list-separator', default=';', help='Separator for complex metadata. [Default: %(default)s]')
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file', required=True, help='Path to the abundance file (format: BIOM).' )
    group_input.add_argument( '-a', '--input-fasta', default=None, required=False, help='Path to the sequences file (format: fasta).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='Path to the output file (format: TSV).')
    args = parser.parse_args()

    cmd_aln = "bwa " + algorithm + \
              " -t " + str(nb_threads) + \
              " " + reference + \
              " " + read1_file + \
              ("" if read2 is None else " " + read2_file) + \
              " > " + tmp_BAM

"samtools view -bS -@ ".$cpu." ".$sam." > ".$unsorted_bam
"samtools sort -@ ".$cpu." ".$unsorted_bam." ".$bam
"mv ".$bam.".bam ".$unmerged_bams[$file_idx]
