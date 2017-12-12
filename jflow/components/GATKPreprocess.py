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

from weaver.function import PythonFunction


def preproWrapper(jvm_options, picard_path, gatk_path, sequencer, is_RNA, reference_genome, in_aln_file, out_aln_file, stderr):
    """
    """
    import os
    import subprocess

    #~ ("" if sample == "None" else " SM=" + sample) + \
    #~ ("" if library == "None" else" LB=" + library) + \
    #~ ("" if platform_unit == "None" else " PU=" + platform_unit) + \
    sample_name = os.path.basename(in_aln_file).split("_")[0]
    library = sample_name
    platform_unit = "UNKNOWN"

    # Add read group
    addGroup_out = out_aln_file + "_01-addGroup.bam"
    cmd = "java " + jvm_options + " -jar " + picard_path + " AddOrReplaceReadGroups" + \
          " SO=coordinate" + \
          " SM=" + sample_name + \
          " LB=" + library + \
          " PL=" + sequencer + \
          " PU=" + platform_unit + \
          " VALIDATION_STRINGENCY=LENIENT" + \
          " I=" + in_aln_file + \
          " O=" + addGroup_out + \
          " 2> " + stderr
    subprocess.check_output(cmd, shell=True)

    # Mark duplications
    dedup_out = out_aln_file + "_02-dedup.bam" if is_RNA == "True" else out_aln_file
    dedup_stat = out_aln_file + "_02-dedup_stat.txt"
    cmd = "java " + jvm_options + " -jar " + picard_path + " MarkDuplicates" + \
          " CREATE_INDEX=true" + \
          " VALIDATION_STRINGENCY=LENIENT" + \
          " I=" + addGroup_out + \
          " O=" + dedup_out + \
          " M=" + dedup_stat + \
          " 2>> " + stderr
    subprocess.check_output(cmd, shell=True)

    # Remove temporary file
    os.remove(addGroup_out)

    # Split exons and reassign one mapping quality
    if is_RNA == "True":
        splitExons_out = out_aln_file
        cmd = "java " + jvm_options + " -jar " + gatk_path + \
              " --analysis_type SplitNCigarReads" + \
              " -U ALLOW_N_CIGAR_READS" + \
              " -rf ReassignOneMappingQuality" + \
              " --reassign_mapping_quality_from 255" + \
              " --reassign_mapping_quality_to 60" + \
              " --reference_sequence " + reference_genome + \
              " --input_file " + dedup_out + \
              " -o " + splitExons_out + \
              " 2>> " + stderr
        subprocess.check_output(cmd, shell=True)

        # Remove temporary file
        os.remove(dedup_out)


class GATKPreprocess (Component):

    def define_parameters(self, aln_files, genome_seq, is_RNA=False, sequencer="ILLUMINA", memory=20):
        # Parameters
        self.add_parameter("is_RNA", "True if sequences aligned come from RNA-Seq.", type="bool", default=is_RNA)
        self.add_parameter("memory", "The memory size to allocate (in GigaBytes).", type="int", default=memory)
        self.add_parameter("sequencer", "The sequencer used.", choices=["ILLUMINA", "SOLID", "LS454", "HELICOS", "PACBIO"], default=sequencer)

        # Input files
        self.add_input_file_list("aln_files", "The alignment files sorted by coordinates (format: BAM).", default=aln_files, required=True)
        self.add_input_file("genome_seq", "Path to the sequences of the reference genome (format: FASTA).", default=genome_seq, required=True)

        # Output files
        self.add_output_file_list("cleaned_files", "The GATK pre-processed files", pattern='{basename_woext}.bam', items=self.aln_files)
        self.add_output_file_list("stderr", "The GATK pre-processed files", pattern='{basename_woext}.stderr', items=self.aln_files)

    def process(self):
        cmd_fmt = '{EXE}' + \
                  ' "-Xmx' + str(self.memory) + 'g"' + \
                  ' "' + self.get_exec_path("picardtools") + '"' + \
                  ' "' + self.get_exec_path("GATK") + '"' + \
                  ' "' + self.sequencer + '"' + \
                  ' "' + str(self.is_RNA) + '"' + \
                  ' "' + self.genome_seq + '"' + \
                  ' {IN} {OUT}'
        prepro_fct = PythonFunction(preproWrapper, cmd_format=cmd_fmt)
        MultiMap(prepro_fct, inputs=[self.aln_files], outputs=[self.cleaned_files, self.stderr], includes=[self.genome_seq])
