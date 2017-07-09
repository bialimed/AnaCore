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


class Cutadapt (Component):

    def define_parameters(self, adaptor_type, adaptor_file, R1, R2=None, error_rate=0.01, min_overlap=None, discard_untrimmed=False):
        # Parameters
        self.add_parameter( "adaptor_type", "***********************", default=adaptor_type, choices=["a", "g", "b"] )
        self.add_parameter( "discard_untrimmed", "***********************", default=discard_untrimmed, type=bool )
        self.add_parameter( "error_rate", "***********************", default=error_rate, type=float )
        self.add_parameter( "min_overlap", "***********************", default=min_overlap, type=int )

        # Files        
        self.add_input_file( "adaptor_file", "***********************", default=adaptor_file, required=True )
        self.add_input_file_list( "in_R1", "************************", default=R1, required=True )
        self.add_input_file_list( "in_R2", "************************", default=R2, required=True )
        self.add_output_file_list( "out_R1", "***********************", pattern='{basename_woext}.fastq.gz', items=self.in_R1 )
        if len(self.in_R2) != 0:
            self.add_output_file_list( "out_R2", "***********************", pattern='{basename_woext}.fastq.gz', items=self.in_R2 )
        self.add_output_file_list( "stderr", "Cutadapt stderr files", pattern='{basename_woext}.stderr', items=self.in_R1 )

    def process(self):
        cmd = self.get_exec_path("cutadapt") + \
            " --error-rate " + str(self.error_rate) + \
            (" --overlap " + str(self.min_overlap) if self.min_overlap == None else "") + \
            (" --discard-untrimmed" if self.discard_untrimmed else "") + \
            " -" + self.adaptor_type + " file:" + self.adaptor_file
        
        if len(self.in_R2) == 0: # Process single read cutadapt
            cmd += " --output $2" + \
                " $1" + \
                " 2> $3"
            cutadapt = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
            MultiMap( cutadapt, inputs=[self.in_R1], outputs=[self.out_R1, self.stderr], includes=[self.adaptor_file] )
        else: # Process paired-end cutadapt
            cmd += " --output $3" + \
                " --paired-output $4" + \
                " $1" + \
                " $2" + \
                " 2> $5"
            cutadapt = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
            MultiMap( cutadapt, inputs=[self.in_R1, self.in_R2], outputs=[self.out_R1, self.out_R2, self.stderr], includes=[self.adaptor_file] )