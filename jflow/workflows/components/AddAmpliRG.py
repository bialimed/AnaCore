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
__copyright__ = 'Copyright (C) 2017 IUCT'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class AddAmpliRG (Component):

    def define_parameters(self, in_regions, in_aln):
        self.add_input_file( "in_regions", "*****", default=in_regions, required=True )    
        self.add_input_file_list( "in_aln", "*****", default=in_aln, required=True )
        self.add_output_file_list( "out_aln", "******************", pattern='{basename_woext}.bam', items=self.in_aln )
        self.add_output_file_list( "stdout", "******************", pattern='{basename_woext}.tsv', items=self.in_aln )
        self.add_output_file_list( "stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.in_aln )

    def process(self):
        cmd = self.get_exec_path("addPanelRG.py") + \
            " --input-panel " + self.in_regions + \
            " --input-aln $1" + \
            " --output-aln $2" + \
            " > $3" + \
            " 2> $4"
        add_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( add_fct, inputs=[self.in_aln], outputs=[self.out_aln, self.stdout, self.stderr], includes=[self.in_regions] )
