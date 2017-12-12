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


class Coverage (Component):

    def define_parameters(self, aln_files, bed_file=None, max_depth=50000):
        # Parameters
        self.add_parameter("max_depth", "Truncate reported depth at a maximum of INT reads.", type="int", default=max_depth)

        # Input files
        self.add_input_file_list("aln_files", "Computes the depth at each position for the specified BAM files (format: BAM).", default=aln_files, required=True)
        self.add_input_file("bed_file", "Compute depth at list of positions or regions in specified file (format: BED).", default=bed_file)

        # Output files
        self.add_output_file_list("depth_files", "The depth for each posiion (format: TSV)", pattern='{basename_woext}_coverage.tsv', items=self.aln_files)
        self.add_output_file_list("stderr", "The path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.aln_files)

    def process(self):
        bed_opt = "-b " + self.bed_file + " " if self.bed_file != None else ""
        depth_opt = "-m " + str(self.max_depth) + " "
        cmd = self.get_exec_path("samtools") + " depth -a " + bed_opt + depth_opt + "$1 > $2 2> $3"
        coverage_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(coverage_fct, inputs=[self.aln_files], outputs=[self.depth_files, self.stderr])
