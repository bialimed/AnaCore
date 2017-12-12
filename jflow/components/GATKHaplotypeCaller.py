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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class GATKHaplotypeCaller (Component):

    def define_parameters(self, aln_files, genome_seq, genome_SNV=None, intervals=None, std_call_conf=30, std_emit_conf=30, memory=30):
        # Parameters
        self.add_parameter("memory", "The memory size to allocate (in GigaBytes).", type="int", default=memory)
        self.add_parameter("std_call_conf", "The minimum phred-scaled confidence threshold at which variants should be called.", type="float", default=std_call_conf)
        self.add_parameter("std_emit_conf", "The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold).", type="float", default=std_emit_conf)

        # Input files
        self.add_input_file_list("aln_files", "The alignment files after preprocess and recalibration (format: BAM).", default=aln_files, required=True)
        self.add_input_file("genome_seq", "Path to the sequences of the reference genome (format: FASTA).", default=genome_seq, required=True)
        self.add_input_file("genome_SNV", "Path to the variants of the reference genome (format: VCF).", default=genome_SNV)
        self.add_input_file("intervals", "Path to the intervals scanned in caller (format: BED).", default=intervals)

        # Output files
        self.add_output_file_list("SNV_files", "The variants files", pattern='{basename_woext}.vcf', items=self.aln_files)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.aln_files)

    def process(self):
        includes_files = [self.genome_seq] if self.genome_SNV == None else [self.genome_seq, self.genome_SNV]
        cmd = "java -Xmx" + str(self.memory) + "g -jar " + self.get_exec_path("GATK") + \
              " --analysis_type HaplotypeCaller" + \
              " --standard_min_confidence_threshold_for_calling " + str(self.std_call_conf) + \
              " --standard_min_confidence_threshold_for_emitting " + str(self.std_emit_conf) + \
              " --annotateNDA" + \
              " --dontUseSoftClippedBases" + \
              " --reference_sequence " + self.genome_seq + \
              ("" if self.genome_SNV == None else " --dbsnp " + self.genome_SNV) + \
              ("" if self.intervals == None else " --intervals " + self.intervals) + \
              " --input_file $1" + \
              " --out $2" + \
              " 2> $3"
        haplotypecaller_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(haplotypecaller_fct, inputs=[self.aln_files], outputs=[self.SNV_files, self.stderr], includes=includes_files)
