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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'


from jflow.workflow import Workflow
from workflows.DSVF import DSVF

class DSVFAnapath (DSVF):
    def get_description(self):
        return "Variant analysis on amplicon double strand librairies."

    def define_parameters(self, parameters_section=None):
        DSVF.define_parameters(self)
        self.add_input_directory( "output_dir", "****************************************", group="Output data" )

    def process(self):
        # store in database
        DSVF.process(self)
        # Convert BAM to CRAM
        print("process")

    def post_process(self):
        DSVF.post_process(self)
        #~ print( self.components )
		# mv files os.rename( sef., )
        print("post-process")
