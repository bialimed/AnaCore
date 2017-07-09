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


class DepthsDistribution (Component):

    def define_parameters(self, in_regions, in_depths, out_format="json", percentile_step=5):
        # Parameters
        self.add_parameter( "percentile_step", "***************", type=int, default=percentile_step )
        self.add_parameter( "out_format", "***************", choices=["tsv", "json"], default=out_format )
        
        # Input files
        self.add_input_file( "in_regions", "********* (format: BED).", default=in_regions, required=True )
        self.add_input_file_list( "in_depths", "********** (format: TSV).", default=in_depths, required=True )

        # Output files
        self.add_output_file_list( "out_depths", "***************** (format: see out_format)", pattern='{basename_woext}.' + self.out_format, items=self.in_depths )
        self.add_output_file_list( "stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.in_depths )

    def process(self):
        cmd = self.get_exec_path("areaCoverage.py") + \
              " --percentile-step " + str(self.percentile_step) + \
              " --input-regions " + self.in_regions + \
              " --inputs-depths $1" + \
              " --output-metrics $2" + \
              " 2> $3"
        convertion = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( convertion, inputs=[self.in_depths], outputs=[self.out_depths, self.stderr], includes=[self.in_regions] )
