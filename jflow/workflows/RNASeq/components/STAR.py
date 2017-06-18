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


def starWrapper(exec_path, genome_dir, genome_gtf, R1, R2, aln_file, stderr, nb_threads):
    """
    """
    import subprocess

    prefix = aln_file.split(".bam")[0]
    cmd = exec_path + \
          " --runThreadN " + nb_threads + \
          " --readFilesCommand zcat" + \
          " --outSAMtype BAM SortedByCoordinate" + \
          " --outFileNamePrefix " + prefix + \
          " --genomeDir " + genome_dir + \
          " --sjdbGTFfile " + genome_gtf + \
          " --readFilesIn " + R1 + " " + R2 + \
          " 2> " + stderr
    subprocess.check_output( cmd, shell=True )
    os.rename(  prefix + "Aligned.sortedByCoord.out.bam", prefix + ".bam" )


class STAR (Component):

    def define_parameters(self, genome_dir, genome_gtf, R1, R2, nb_threads=1):
        self.add_parameter("genome_dir", "Path to the directory containing the STAR index of the reference genome.", default=genome_dir, required=True)
        self.add_input_file("genome_gtf", "Path to the annotation of the reference genome (format: GTF).", default=genome_gtf, required=True)
        self.add_input_file_list( "R1", "Which R1 files should be used.", default=R1, required=True )
        self.add_input_file_list( "R2", "Which R2 files should be used.", default=R2, required=True )
        self.add_parameter("nb_threads", "Number of threads used", type=int, default=nb_threads)
        self.add_output_file_list("aln_files", "The STAR outputed file", pattern='{basename_woext}.bam', items=self.R1)
        self.add_output_file_list("stderr", "The STAR stderr file", pattern='{basename_woext}.stderr', items=self.R1)

    def process(self):
        star = PythonFunction(
            starWrapper,
            cmd_format="{EXE} " + self.get_exec_path("STAR") + ' "' + self.genome_dir + '" "' + self.genome_gtf + '" ' + " {IN} {OUT} " + str(self.nb_threads) )
        MultiMap( star, inputs=[self.R1, self.R2], outputs=[self.aln_files, self.stderr], includes=[self.genome_gtf] )
