#
# Copyright (C) 2015 INRA
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
from weaver.function import ShellFunction


class BWAIndex (Component):

    def define_parameters(self, input_fasta, algorithm="bwtsw"):
        # Parameters
        self.add_parameter("algorithm", "Which algorithm should be used to index the fasta file", default=algorithm, choices=["bwtsw", "div", "is"])

        # Files
        self.add_input_file("input_fasta", "Which fasta file should be indexed", file_format="fasta", default=input_fasta, required=True)
        self.add_output_file("databank", "The indexed databank", filename=os.path.basename(input_fasta))
        self.add_output_file("stdout", "The BWAIndex stdout file", filename="bwaindex.stdout")
        self.add_output_file("stderr", "The BWAIndex stderr file", filename="bwaindex.stderr")

    def process(self):
        cmd = "ln -s $1 $2; " + self.get_exec_path("bwa") + " index -a " + self.algorithm + " -p $2 $1 > $3 2> $4"
        index_fct = ShellFunction(cmd, cmd_format="{EXE} {IN} {OUT}")
        index_fct(inputs=self.input_fasta, outputs=[self.databank, self.stdout, self.stderr])
