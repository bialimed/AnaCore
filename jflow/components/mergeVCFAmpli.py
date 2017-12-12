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
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class MergeVCFAmpli (Component):

    def define_parameters(self, libA_regions_wout_primers, libB_regions_wout_primers, libA_variants, libB_variants, libA_aln, libB_aln):
        # Input files
        self.add_input_file("libA_regions_wout_primers", "The path to the amplicons design for the first library. The start and end of the amplicons must be without primers (format: BED).", default=libA_regions_wout_primers, required=True)
        self.add_input_file("libB_regions_wout_primers", "The path to the amplicons design for the second library. The start and end of the amplicons must be without primers (format: BED).", default=libB_regions_wout_primers, required=True)
        self.add_input_file_list("libA_variants", "The path to the variants file produced from the first library (format: VCF).", default=libA_variants, required=True)
        self.add_input_file_list("libB_variants", "The path to the variants file produced from the first second (format: VCF).", default=libB_variants, required=True)
        self.add_input_file_list("libA_aln", "The path to the alignments file produced from the first library (format: BAM).", default=libA_aln, required=True)
        self.add_input_file_list("libB_aln", "The path to the alignments file produced from the second library (format: BAM).", default=libB_aln, required=True)

        # Output files
        self.add_output_file_list("out_variants", "The path to the outputted file (format: VCF).", pattern='{basename_woext}.vcf', items=self.libA_variants)
        self.add_output_file_list("stderr", "The path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.libA_variants)

    def process(self):
        cmd = self.get_exec_path("mergeVCFAmpli.py") + \
            " --input-designs " + self.libA_regions_wout_primers + " " + self.libB_regions_wout_primers + \
            " --input-variants $1 $2" + \
            " --input-aln $3 $4" + \
            " --output-variants $5" + \
            " 2> $6"
        merge_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(
            merge_fct,
            inputs=[self.libA_variants, self.libB_variants, self.libA_aln, self.libB_aln],
            outputs=[self.out_variants, self.stderr],
            includes=[self.libA_regions_wout_primers, self.libB_regions_wout_primers]
        )
