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

import os

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class Salmon (Component):

    def define_parameters(self, transcriptome_idx, R1, R2, libtype="A", nb_threads=1):
        # Parameters
        self.add_parameter("libtype", "The library type string consists of three parts: the relative orientation of the reads, the strandedness of the library, and the directionality of the reads (see http://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype).", default=libtype)
        self.add_parameter("nb_threads", "Number of threads used", type=int, default=nb_threads)

        # Input files
        self.add_input_file("transcriptome_idx", "Path to the transcriptome index (format: salmon index).", default=transcriptome_idx, required=True)
        self.add_input_file_list("R1", "Which R1 files should be used.", default=R1, required=True)
        self.add_input_file_list("R2", "Which R2 files should be used.", default=R2, required=True)

        # Output files
        self.add_output_file_list("quant_files", "The salmon outputed file", pattern='{basename_woext}/quant.sf', items=self.R1)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.R1)

    def process(self):
        quant_folders = [os.path.dirname(curr_file) for curr_file in self.quant_files]
        cmd = self.get_exec_path("salmon") + \
            " quant" + \
            " --threads " + str(self.nb_threads) + \
            " --libType " + self.libtype + \
            " -i " + self.transcriptome_idx + \
            " -1 $1" + \
            " -2 $2" + \
            " -o $3" + \
            " 2> $4"
        salmon_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(salmon_fct, inputs=[self.R1, self.R2], outputs=[quant_folders, self.stderr])
