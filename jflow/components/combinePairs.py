#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap
from weaver.function import ShellFunction


class CombinePairs (Component):

    def define_parameters(self, R1, R2, names=None, mismatch_ratio=0.25, min_overlap=20, min_frag_length=None, max_frag_length=None):
        # Parameters
        self.add_parameter("max_frag_length", "Maximum length for the resulting fragment. This filter is applied after best overlap selection.", default=max_frag_length, type=int)
        self.add_parameter("min_frag_length", "Minimum length for the resulting fragment. This filter is applied after best overlap selection.", default=min_frag_length, type=int)
        self.add_parameter("min_overlap", "The minimum required overlap length between two reads to provide a confident overlap.", default=min_overlap, type=int)
        self.add_parameter("mismatch_ratio", "Maximum allowed ratio between the number of mismatched base pairs and the overlap length. Two reads will not be combined with a given overlap if that overlap results in a mismatched base density higher than this value.", default=mismatch_ratio, type=float)
        self.add_parameter_list("names", "The basenames of the output fastq in order of the R1. By default the basename is automatically determined.", default=names)
        if len(self.names) == 0:
            self.prefixes = self.get_outputs('{basename_woext}', [R1, R2])
        else:
            self.prefixes = self.get_outputs('{basename_woext}', names)

        # Inputs files
        self.add_input_file_list("R1", "fastq read R1 (format: fastq).", default=R1, required=True)
        self.add_input_file_list("R2", "fastq read R2 (format: fastq).", default=R2, required=True)

        # Outputs files
        self.add_output_file_list("out_combined", "Pathes to the files containing combined pairs (format: fastq).", pattern='{basename_woext}_combined.fastq.gz', items=self.prefixes)
        self.add_output_file_list("out_report", "Pathes to the files containing combination metrics (format: JSON).", pattern='{basename_woext}_report.json', items=self.prefixes)
        self.add_output_file_list("stderr", "Pathes to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.prefixes)


    def process(self):
        cmd = self.get_exec_path("combinePairs.py") + \
            ("" if self.max_frag_length == None else " --max-frag-length " + str(self.max_frag_length)) + \
            ("" if self.min_frag_length == None else " --min-frag-length " + str(self.min_frag_length)) + \
            " --min-overlap " + str(self.min_overlap) + \
            " --max-contradict-ratio " + str(self.mismatch_ratio) + \
            " --input-R1 $1" + \
            " --input-R2 $2" + \
            " --output-combined $3" + \
            " --output-report $4" + \
            " 2> $5"
        combinePairs_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(
            combinePairs_fct,
            inputs=[self.R1, self.R2],
            outputs=[self.out_combined, self.out_report, self.stderr]
        )
