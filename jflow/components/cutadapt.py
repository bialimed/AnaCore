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
__version__ = '1.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class Cutadapt (Component):

    def define_parameters(self, adaptor_type, adaptor_file, R1, R2=None, error_rate=0.1, min_overlap=None, discard_untrimmed=False):
        # Parameters
        self.add_parameter("adaptor_type", "Location of adaptor a: 3' end, g: 5' end or b: 3' or 5' end.", default=adaptor_type, choices=["a", "g", "b"])
        self.add_parameter("discard_untrimmed", "With this option reads that do not contain the adapter are discarded.", default=discard_untrimmed, type=bool)
        self.add_parameter("error_rate", "Maximum allowed error rate (no. of errors divided by the length of the matching region).", default=error_rate, type=float)
        self.add_parameter("min_overlap", "If the overlap between the read and the adapter is shorter than this value, the read is not modified.", default=min_overlap, type=int)

        # Input Files
        self.add_input_file("adaptor_file", "Path to the adapter(s) sequence(s) file (format: fasta).", default=adaptor_file, required=True)
        self.add_input_file_list("in_R1", "Path to the R1 files (format: fastq).", default=R1, required=True)
        self.add_input_file_list("in_R2", "Path to the R2 files (format: fastq).", default=R2, required=True)

        # Output Files
        self.add_output_file_list("out_R1", "Path to the outputted R1 files (format: fastq).", pattern='{basename_woext}_trim.fastq.gz', items=self.in_R1)
        if len(self.in_R2) != 0:
            self.add_output_file_list("out_R2", "Path to the outputted R2 files (format: fastq).", pattern='{basename_woext}_trim.fastq.gz', items=self.in_R2)
        self.add_output_file_list("stdout", "Path to the stdout files (format: txt).", pattern='{basename_woext}.stdout', items=self.in_R1)
        self.add_output_file_list("stderr", "Path to the stderr files (format: txt).", pattern='{basename_woext}.stderr', items=self.in_R1)

    def process(self):
        cmd = self.get_exec_path("cutadapt") + \
            " --error-rate " + str(self.error_rate) + \
            (" --overlap " + str(self.min_overlap) if self.min_overlap == None else "") + \
            (" --discard-untrimmed" if self.discard_untrimmed else "") + \
            " -" + self.adaptor_type + " file:" + self.adaptor_file

        if len(self.in_R2) == 0:  # Process single read cutadapt
            cmd += " --output $2" + \
                " $1" + \
                " > $3" + \
                " 2> $4"
            cutadapt_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
            MultiMap(cutadapt_fct, inputs=[self.in_R1], outputs=[self.out_R1, self.stdout, self.stderr], includes=[self.adaptor_file])
        else:  # Process paired-end cutadapt
            cmd += " --output $3" + \
                " --paired-output $4" + \
                " $1" + \
                " $2" + \
                " > $5" + \
                " 2> $6"
            cutadapt_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
            MultiMap(cutadapt_fct, inputs=[self.in_R1, self.in_R2], outputs=[self.out_R1, self.out_R2, self.stdout, self.stderr], includes=[self.adaptor_file])
