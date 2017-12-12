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


class FreeBayes (Component):

    def define_parameters(self, aln_files, ref_seq, min_alt_count=5, min_alt_fraction=0.05):
        # Parameters
        self.add_parameter("min_alt_count", "Only position with this number of reads supporting an alternate allele within a single individual are evaluated.", type=int, default=min_alt_count)
        self.add_parameter("min_alt_fraction", "Only position with this fraction of reads supporting an alternate allele within a single individual are evaluated.", type=float, default=min_alt_fraction)

        # Input files
        self.add_input_file("ref_seq", "Path to the reference sequences (format: fasta).", default=ref_seq)
        self.add_input_file_list("aln_files", "Path to the alignments path (format: BAM).", default=aln_files)

        # Output files
        self.add_output_file_list("variants_files", "Path to the variants files (format: VCF).", pattern='{basename_woext}.vcf', items=self.aln_files)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}_stderr.txt', items=self.aln_files)

    def process(self):
        freebayes_cmd = self.get_exec_path("freebayes") + \
             ' --min-alternate-count ' + str(self.min_alt_count) + \
             ' --min-alternate-fraction ' + str(self.min_alt_fraction) + \
             ' --f ' + self.ref_seq + \
             ' $1' + \
             ' > $2' + \
             ' 2> $3'
        freebayes_fct = ShellFunction(freebayes_cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(freebayes_fct, inputs=[self.aln_files], outputs=[self.variants_files, self.stderr], includes=[self.ref_seq])
