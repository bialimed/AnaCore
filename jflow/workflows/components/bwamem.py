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

from weaver.function import PythonFunction


def bwa( bwa_exec_path, samtools_exec_path, reference_genome, R1, R2, aln_file, stderr ):
    """
    """
    import os
    import time
    import subprocess

    tmp_aln = os.path.join( os.path.dirname(aln_file), "tmp_" + str(time.time()) + ".sam" )

    # Alignment
    cmd_aln = bwa_exec_path + \
               ' mem' + \
               ' "' + reference_genome + '"' + \
               ' "' + R1 + '"' + \
               ' "' + R2 + '"' + \
               ' > "' + tmp_aln + '"' + \
               ' 2> ' + stderr
    subprocess.check_output( cmd_aln, shell=True )

    # Sort and conversion to BAM
    cmd_sort = samtools_exec_path + \
                ' sort' + \
                ' -O BAM' + \
                ' -o "' + aln_file + '"' + \
                ' "' + tmp_aln + '"' + \
                ' 2>> ' + stderr
    subprocess.check_output( cmd_sort, shell=True )

    # Remove mp
    os.remove( tmp_aln )


class BWAmem (Component):

    def define_parameters(self, reference_genome, R1, R2=None): ############################ pb without R2 passage nb variable param
        self.add_input_file_list( "R1", "Which reads files should be used. ***************************", default=R1, required=True )
        self.add_input_file_list( "R2", "Which reads files should be used. ***************************", default=R2, required=True )
        self.add_input_file("reference_genome", "Which reference file should be used", default=reference_genome, required=True)
        self.add_output_file_list("aln_files", "The BWA outputed file", pattern='{basename_woext}.bam', items=self.R1)
        self.add_output_file_list("stderr", "The BWA stderr file", pattern='{basename_woext}.stderr', items=self.R1)

    def process(self):
        bwamem = PythonFunction( bwa, cmd_format='{EXE} ' + self.get_exec_path("bwa") + ' ' + self.get_exec_path("samtools") + ' "' + self.reference_genome + '" {IN} {OUT}')
        MultiMap(bwamem, inputs=[self.R1, self.R2], outputs=[self.aln_files, self.stderr], includes=[self.reference_genome])
