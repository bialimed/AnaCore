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


class InterOp (Component):

    def define_parameters(self, sequencer_run_folder):
        # Input files
        self.add_parameter("sequencer_run_folder", "Path to the Illumina's run folder. The files RunInfo.xml and runParameters.xml must exist at the same level as Interop folder.", default=sequencer_run_folder, required=True)

        # Output files
        self.add_output_file("interop_text", "The Illumina's run metrics (format: txt).", filename='interop.txt')
        self.add_output_file("interop_stderr", "The stderr files (format: txt).", filename='dumptext.stderr')
        self.add_output_file("summary_file", "The Illumina's run summary (format: JSON)", filename='interop.json')
        self.add_output_file("summary_stderr", "The stderr files (format: txt).", filename='summary.stderr')

    def process(self):
        # Binaries to txt
        dumptext_cmd = self.get_exec_path("dumptext") + \
                       ' "' + self.sequencer_run_folder + '"' + \
                       ' > $1' + \
                       ' 2> $2'
        dumptext_fct = ShellFunction(dumptext_cmd, cmd_format='{EXE} {OUT}')
        dumptext_fct(outputs=[self.interop_text, self.interop_stderr])

        # Summary
        summary_cmd = self.get_exec_path("interopSummary.py") + \
                      ' --interop-file $1' + \
                      ' --output-file $2' + \
                      ' 2> $3'
        summary_fct = ShellFunction(summary_cmd, cmd_format='{EXE} {IN} {OUT}')
        summary_fct(inputs=[self.interop_text], outputs=[self.summary_file, self.summary_stderr])
