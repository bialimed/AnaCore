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


class Demultiplex (Component):

    def define_parameters(self, manifest, R1, R2, error_rate=0.01):
        # Parameters
        self.add_parameter( "error_rate", "***********************", default=error_rate, type=float )

        # Files        
        self.add_input_file( "manifest", "**************************************", default=manifest, required=True )
        self.add_input_file_list( "in_R1", "Which reads files should be used. ***************************", default=R1, required=True )
        self.add_input_file_list( "in_R2", "Which reads files should be used. ***************************", default=R2, required=True )

        self.add_output_file_endswith( "out_R1", "******************", pattern='R1.fastq.gz' )
        self.add_output_file_endswith( "out_R2", "******************", pattern='R2.fastq.gz' )
        self.add_output_file_list( "stderr", "Demultiplex stderr files", pattern='{basename_woext}.stderr', items=self.in_R1 )

    def process(self):
        cmd = self.get_exec_path("demultiplex") + \
            " --error-rate " + str(self.error_rate) + \
            " --manifest-path " + self.manifest + \
            " --R1-path $1" + \
            " --R2-path $2" + \
            " --output-dir " + self.output_directory + \
            " 2> $3"
        demultiplex = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( demultiplex, inputs=[self.in_R1, self.in_R2], outputs=[self.stderr], includes=[self.manifest] )
