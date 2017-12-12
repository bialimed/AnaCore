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

from weaver.function import PythonFunction


def bam_index(exec_path, in_aln, out_aln, out_idx, stderr_path):
    import os
    import subprocess
    # Create symlink
    if os.path.exists(out_aln):
        os.remove(out_aln)
    os.symlink(in_aln, out_aln)
    # Submit cmd
    cmd = exec_path + \
        " index" + \
        " " + out_aln + \
        " 2> " + stderr_path
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, shell=True)


class BAMIndex (Component):

    def define_parameters(self, in_aln):
        self.add_input_file_list("in_aln", "PATH to the alignment file (format: BAM).", default=in_aln, required=True)
        self.add_output_file_list("out_aln", "Path to the outputted alignments file (format: BAM).", pattern='{basename_woext}.bam', items=self.in_aln)
        self.add_output_file_list("out_idx", "Path to the outputted index (format: BAI).", pattern='{basename_woext}.bam.bai', items=self.in_aln)
        self.add_output_file_list("stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.in_aln)

    def process(self):
        index_fct = PythonFunction(bam_index, cmd_format='{EXE} ' + self.get_exec_path("samtools") + ' {IN} {OUT}')
        MultiMap(index_fct, inputs=[self.in_aln], outputs=[self.out_aln, self.out_idx, self.stderr])
