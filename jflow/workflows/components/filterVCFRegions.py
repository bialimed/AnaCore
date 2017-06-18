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

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class FilterVCFRegion (Component):

    def define_parameters(self, in_regions, in_variants):
        self.add_input_file( "in_regions", "*****", default=in_regions, required=True )    
        self.add_input_file_list( "in_variants", "*****", default=in_variants, required=True )
        self.add_output_file_list( "out_variants", "******************", pattern='{basename_woext}.vcf', items=self.in_variants )
        self.add_output_file_list( "stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.in_variants )

    def process(self):
		cmd = self.get_exec_path("addPanelRG.py") + \
            " --selected-regions " + self.in_regions + \
            " --input-variants $1" + \
            " --output-variants $2" + \
            " 2> $3"
		filter_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
		MultiMap(filter_fct, inputs=[self.in_variants], outputs=[self.out_variants, self.stderr], includes=[self.in_regions] )
