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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class AmpliVariantCalling (Component):

    def define_parameters(self, in_genome_seq, in_regions_with_primers, in_regions_wout_primers, non_overlapping_design, in_aln, min_AF=0.02, min_base_qual=25):
        # Parameters
        self.add_parameter("min_AF", "Variants with an allele frequency under this value are not emitted.", type=float, default=min_AF, required=True)
        self.add_parameter("min_base_qual", "The phred score for a base to be considered a good call.", type=int, default=min_base_qual, required=True)

        # Input files
        self.add_input_file_list("in_aln", "Path to the alignment file (format: BAM).", default=in_aln, required=True)
        self.add_input_file("in_genome_seq", "Path to the reference genome (format: fasta).", default=in_genome_seq, required=True)
        self.add_input_file("in_regions_with_primers", "Path to the amplicons design with their primers (format: BED).", default=in_regions_with_primers, required=True)
        self.add_input_file("in_regions_wout_primers", "Path to the amplicons design without their primers (format: BED).", default=in_regions_wout_primers, required=True)
        self.add_input_file("non_overlapping_design", "Path to the list of amplicons (format: TSV). The first column is the ID of the amplicon and the second is the name of the group where the amplicon has no overlap with other amplicons of this group.", default=non_overlapping_design, required=True)

        # Output files
        self.add_output_file_list("out_variants", "Path to the outputted variants file (format: VCF).", pattern='{basename_woext}.vcf', items=self.in_aln)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.in_aln)
        self.add_output_file_list("stdout", "Path to the stdout file (format: txt).", pattern='{basename_woext}.stdout', items=self.in_aln)

    def process(self):
        cmd = self.get_exec_path("ampliVariantCalling.py") + \
            " --min-AF " + str(self.min_AF) + \
            " --min-base-qual " + str(self.min_base_qual) + \
            " --input-genome " + self.in_genome_seq + \
            " --input-design-with-primers " + self.in_regions_with_primers + \
            " --input-design-wout-primers " + self.in_regions_wout_primers + \
            " --input-non-overlapping-design " + self.non_overlapping_design + \
            " --input-aln $1" + \
            " --output-variants $2" + \
            " --output-log $3" + \
            " 2> $4"
        calling_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(
            calling_fct,
            inputs=[self.in_aln],
            outputs=[self.out_variants, self.stdout, self.stderr],
            includes=[self.in_genome_seq, self.in_regions_with_primers, self.in_regions_wout_primers, self.non_overlapping_design]
        )
