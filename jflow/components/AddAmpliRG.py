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
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class AddAmpliRG (Component):

    def define_parameters(self, in_regions, in_aln, summary_format="tsv", min_zoi_cov=10):
        # Parameters
        self.add_parameter("min_zoi_cov", 'The minimum cumulative length of reads pair in zone of interest. If the number of nucleotids coming from R1 on ZOI + the number of nucleotids coming from R2 on ZOI is lower than this value the pair is counted in "only_primers".', type=int, default=min_zoi_cov)
        self.add_parameter("summary_format", "The format of the stdout.", choices=["tsv", "json"], default=summary_format)

        # Input files
        self.add_input_file("in_regions",
                            "Path to the list of amplicons with their primers (format: BED). Each area must have an unique ID in the name field and a strand.",
                            default=in_regions,
                            required=True)
        self.add_input_file_list("in_aln",
                                 "The path to the alignments files (format: BAM). This file must be sorted by coordinates.",
                                 default=in_aln,
                                 required=True)

        # Output files
        self.add_output_file_list("out_aln", "The path to the alignments files (format: BAM).", pattern='{basename_woext}.bam', items=self.in_aln)
        self.add_output_file_list("out_summary",
                                  "The summary files (format: see summary_format). These files contain information about the number of reads out off target, reversed and valid.",
                                  pattern='{basename_woext}.' + self.summary_format,
                                  items=self.in_aln)
        self.add_output_file_list("stderr", "The path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.in_aln)

    def process(self):
        cmd = self.get_exec_path("addAmpliRG.py") + \
            " --min-zoi-cov " + str(self.min_zoi_cov) + \
            " --summary-format " + self.summary_format + \
            " --input-panel " + self.in_regions + \
            " --input-aln $1" + \
            " --output-aln $2" + \
            " --output-summary $3" + \
            " 2> $4"
        add_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(add_fct, inputs=[self.in_aln], outputs=[self.out_aln, self.out_summary, self.stderr], includes=[self.in_regions])
