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

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class VarDictAmpli (Component):

    def define_parameters(self, in_genome_seq, in_regions, in_aln, min_AF=0.02):
        # Parameters
        self.add_parameter( "min_AF", "Variants with an allele frequency under this value are not emitted.", type=float, default=min_AF, required=True )

        # Input files
        self.add_input_file_list( "in_aln", "Path to the alignment file (format: BAM).", default=in_aln, required=True )
        self.add_input_file( "in_genome_seq", "The path to the reference genome (format: fasta).", default=in_genome_seq, required=True )
        self.add_input_file( "in_regions", "Path to the amplicons design with their primers if the interest area is located in fields thickStart and thickEnd (format: BED).", default=in_regions, required=True )

        # Output files
        self.add_output_file_list( "out_calling", "******************", pattern='{basename_woext}_calling.tsv', items=self.in_aln )
        self.add_output_file_list( "out_filter", "******************", pattern='{basename_woext}_filter.tsv', items=self.in_aln )
        self.add_output_file_list( "out_variants", "Path to the outputted variants file (format: VCF).", pattern='{basename_woext}.vcf', items=self.in_aln )
        self.add_output_file_list( "calling_stderr", "The stderr file of the calling step (format: tsv).", pattern='{basename_woext}_calling.stderr', items=self.in_aln )
        self.add_output_file_list( "filter_stderr", "The stderr file of the filter step (format: tsv).", pattern='{basename_woext}_filter.stderr', items=self.in_aln )
        self.add_output_file_list( "convert_stderr", "The stderr file of the conversion step (format: tsv).", pattern='{basename_woext}_convert.stderr', items=self.in_aln )

    def process(self):
        # Variant calling
        cmd = self.get_exec_path("vardict.pl") + \
            " -f " + str(self.min_AF) + \
            " -F 0" + \
            " -c 1 -S 2 -E 3 -g 4" + \
            " -b $1 " + \
            " -G " + self.in_genome_seq + \
            " " + self.in_regions + \
            " > $2" + \
            " 2> $3"
        calling_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( calling_fct, inputs=[self.in_aln], outputs=[self.out_calling, self.calling_stderr], includes=[self.in_genome_seq, self.in_regions] )

        # Filter on strand bias
        cmd = "cat $1 |" + \
            " " + self.get_exec_path("teststrandbias.R") + \
            " > $2" + \
            " 2> $3"
        filter_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( filter_fct, inputs=[self.out_calling], outputs=[self.out_filter, self.filter_stderr] )

        # Convert to VCF
        cmd = "cat $1 |" + \
            " " + self.get_exec_path("var2vcf_valid.pl") + \
            " -a" + \
            " -E" + \
            " -f " + str(self.min_AF) + \
            " > $2" + \
            " 2> $3"
        convert_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( convert_fct, inputs=[self.out_filter], outputs=[self.out_variants, self.convert_stderr] )
