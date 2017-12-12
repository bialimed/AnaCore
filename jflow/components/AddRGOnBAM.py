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

from weaver.function import ShellFunction


class AddRGOnBAM (Component):

    def define_parameters(self, in_aln, platforms, samples, libraries):
        # Parameters
        self.add_parameter_list("libraries", "Libraries.", default=libraries, required=True)
        self.add_parameter_list("platforms", "Platforms/technologies used to produce the reads.", choices=["CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "ONT", "PACBIO"], default=platforms, required=True)
        self.add_parameter_list("samples", "Samples. Use pool name where a pool is being sequenced.", default=samples, required=True)

        # Files
        self.add_input_file_list("in_aln", "Path to the alignments files (format: BAM).", default=aln_files, required=True)
        self.add_output_file_list("out_aln", "Path to the outputted alignments file (format: BAM).", pattern='{basename_woext}.bam', items=self.in_aln)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.in_aln)

    def process(self):
        for idx_aln, curr_aln in enumerate(self.in_aln):
            cmd = self.get_exec_path("addRGPOnBAM.py") + \
                " --pl " + self.platforms[idx_aln] + \
                " --sm " + self.samples[idx_aln] + \
                " --lb " + self.libraries[idx_aln] + \
                " --input-aln $1" + \
                " --output-aln $2" + \
                " 2> $3"
            addRG = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
            addRG(inputs=[self.in_aln[idx_aln]], outputs=[self.out_aln[idx_aln], self.stderr[idx_aln]])
