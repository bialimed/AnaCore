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


class VariantsCtrlCheck (Component):

    def define_parameters(self, expected_file, variants_files, out_format="tsv"):
        # Parameters
        self.add_parameter("out_format", "The format of the output.", choices=["tsv", "json"], default=out_format)

        # Input files
        self.add_input_file("expected_file", "The expected variants in control samples (format: VCF).", default=expected_file)
        self.add_input_file_list("variants_files", "The variants found in control samples (format: VCF).", default=variants_files, required=True)

        # Output files
        self.add_output_file_list("eval_files", "The divergence between expected and found variants (format: see out_format).", pattern='{basename_woext}_eval.' + self.out_format, items=self.variants_files)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.variants_files)

    def process(self):
        eval_cmd = self.get_exec_path("evalVariantControl.py") + \
             ' --expected-file "' + self.expected_file + '"' + \
             ' --detected-file $1' + \
             ' --output-file $2' + \
             ' 2> $3'
        eval_fct = ShellFunction(eval_cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(eval_fct, inputs=[self.variants_files], outputs=[self.eval_files, self.stderr])
