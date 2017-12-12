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


class CigarlineGraph (Component):

    def define_parameters(self, aln_files):
        self.add_input_file_list("aln_files", "Path to the alignment file (format: BAM).", default=aln_files, required=True)
        self.add_output_file_list("cigar_R1", "Path to the data outputted for R1 (format: TSV).", pattern='{basename_woext}_cigarR1.tsv', items=self.aln_files)
        self.add_output_file_list("cigar_R2", "Path to the data outputted for R2 (format: TSV).", pattern='{basename_woext}_cigarR2.tsv', items=self.aln_files)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.aln_files)

    def process(self):
        cmd = self.get_exec_path("samtools") + ' view -F256 $1 | ' + self.get_exec_path("cigarlineGraph") + " --readsplit -i - -t $2 $3 2> $4"
        cigarline_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(cigarline_fct, inputs=[self.aln_files], outputs=[self.cigar_R1, self.cigar_R2, self.stderr])
