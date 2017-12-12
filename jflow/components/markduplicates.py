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


class MarkDuplicates (Component):

    def define_parameters(self, aln_files, remove_seq_dup=False, remove_dup=False, memory=5, validation_stringency="lenient"):
        # Parameters
        self.add_parameter("memory", "The memory size to allocate (in GigaBytes).", type="int", default=memory)
        self.add_parameter("remove_dup", "If True the duplicates (sequencing and optical) are removed from the alignments.", type=bool, default=remove_dup)
        self.add_parameter("remove_seq_dup", "If True the sequencing duplicates are removed from the alignments.", type=bool, default=remove_seq_dup)
        self.add_parameter("validation_stringency", "Validation stringency for alignments files.", choices=["strict", "lenient", "silent"], default=validation_stringency)

        # Input files
        self.add_input_file_list("aln_files", "Path to the alignments files (format: BAM).", default=aln_files)

        # Output files
        self.add_output_file_list("marked_files", "Path to the marked/cleaned alignments files (format: BAM).", pattern='{basename_woext}.bam', items=self.aln_files)
        self.add_output_file_list("metrics_files", "Path to the duplication metrics (format: txt).", pattern='{basename_woext}_metrics.txt', items=self.aln_files)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}_stderr.txt', items=self.aln_files)

    def process(self):
        mark_cmd = 'java -Xmx' + self.memory + 'G -jar ' + self.get_exec_path("picard.jar") + \
            ' MarkDuplicates' + \
            ' CREATE_INDEX=true' + \
            ' VALIDATION_STRINGENCY=' + self.validation_stringency.upper() + \
            ' REMOVE_SEQUENCING_DUPLICATES=' + str(self.remove_seq_dup).upper() + \
            ' REMOVE_DUPLICATES=' + str(self.remove_dup).upper() + \
            ' INPUT=$1' + \
            ' OUTPUT=$2' + \
            ' METRICS_FILE=$3' + \
            ' 2> $4'
        mark_fct = ShellFunction(mark_cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(mark_fct, inputs=[self.aln_files], outputs=[self.marked_files, self.metrics_files, self.stderr])
