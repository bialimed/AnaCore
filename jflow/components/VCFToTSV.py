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


class VCFToTSV (Component):

    def define_parameters(self, in_variants):
        self.add_input_file_list("in_variants", "Path to the files containing variants annotated with VEP v88+ (format: VCF).", default=in_variants, required=True)
        self.add_output_file_list("out_variants", "The path to the outputted variants files (format: TSV).", pattern='{basename_woext}.tsv', items=self.in_variants)
        self.add_output_file_list("stderr", "The path to the stderr files (format: txt).", pattern='{basename_woext}.stderr', items=self.in_variants)

    def process(self):
        cmd = self.get_exec_path("VCFToSV.py") + \
            " --input-variants $1" + \
            " > $2" + \
            " 2> $3"
        convert_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(convert_fct, inputs=[self.in_variants], outputs=[self.out_variants, self.stderr])
