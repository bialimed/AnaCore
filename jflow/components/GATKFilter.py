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


class GATKFilter (Component):

    def define_parameters(self, SNV_files, reference_seq, filters=None, discard=True, memory=10):
        if filters == None:
            filters = [
                ' -window 35 -cluster 3',
                ' -filterName strdBias -filter "FS > 30.0"',
                ' -filterName lowQD -filter "QD < 2.0"',
                ' -filterName LowDP -filter "DP < 100.0"'
            ]

        # Parameters
        self.add_parameter("discard", "True if filtered variants are removed.", type="bool", default=discard)
        self.add_parameter_list("filters", "The list of filters applied on variants.", default=filters)
        self.add_parameter("memory", "The memory size to allocate (in GigaBytes).", type="int", default=memory)

        # Input files
        self.add_input_file_list("SNV_files", "Path to the variants files (format: VCF).", default=SNV_files, required=True)
        self.add_input_file("reference_seq", "Path to the sequences of the reference (format: FASTA).", default=reference_seq, required=True)

        # Output files
        self.add_output_file_list("marked_SNV_files", "The variants files containing all variants and with the filter field filled.", pattern='{basename_woext}.vcf', items=self.SNV_files)
        self.add_output_file_list("marked_stderr", "The VariantFiltration error files.", pattern='{basename_woext}.stderr', items=self.SNV_files)
        if discard:
            self.add_output_file_list("filtered_SNV_files", "The filtered variants files.", pattern='{basename_woext}_filtered.vcf', items=self.SNV_files)
            self.add_output_file_list("filtered_stderr", "The SelectVariants error files.", pattern='{basename_woext}_filtered.stderr', items=self.SNV_files)

    def process(self):
        # Mark filtered
        cmd = "java -Xmx" + str(self.memory) + "g -jar " + self.get_exec_path("GATK") + \
              " --analysis_type VariantFiltration" + \
              " ".join(self.filters) + \
              " --reference_sequence " + self.reference_seq + \
              " --variant $1" + \
              " --out $2" + \
              " 2> $3"
        tag_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(tag_fct, inputs=[self.SNV_files], outputs=[self.marked_SNV_files, self.marked_stderr], includes=self.reference_seq)

        # Discard no passing filters
        if self.discard:
            cmd = "java -Xmx" + str(self.memory) + "g -jar " + self.get_exec_path("GATK") + \
                  " --analysis_type SelectVariants" + \
                  " --excludeFiltered" + \
                  " --reference_sequence " + self.reference_seq + \
                  " --variant $1" + \
                  " --out $2" + \
                  " 2> $3"
            filter_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
            MultiMap(filter_fct, inputs=[self.marked_SNV_files], outputs=[self.filtered_SNV_files, self.filtered_stderr], includes=self.reference_seq)
