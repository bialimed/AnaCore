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

def tophatWrapper(exec_path, nb_threads, anchor_length, genome_index_base, genome_gtf, aln_file, stderr, R1, R2=None):
    import os
    import subprocess

    output_dir = os.path.dirname( aln_file )
    #~ tmp_dir
    #~ sample
    cmd = exec_path + \
          " --fusion-search" + \
          " --num-threads " + nb_threads + \
          " --bowtie1" + \
          " --fusion-anchor-length " + str(anchor_length) + \
          " --output-dir " + output_dir + \
          (" --GTF " + genome_gtf if genome_gtf != "no_gtf" else "") + \
          " " + genome_index_base + \
          " " + R1 + \
          (" " + R2 if R2 is not None else "") + \
          " 2> " + stderr
    subprocess.check_output( cmd, shell=True )
    # Move files
    # accepted hits.bam, deletions.bed, insertions.bed, junctions.bed, fusions.out
    #~ os.rename( os.path.join(tmp_dir, "fusions.out"), os.path.join(output_dir, sample + "_fusions.out") )
    #~ os.rename( os.path.join(tmp_dir, "fusions.out"), os.path.join(output_dir, sample + "_fusions.out") )
    #~ os.rename( os.path.join(tmp_dir, "fusions.out"), os.path.join(output_dir, sample + "_fusions.out") )
    #~ os.rename( os.path.join(tmp_dir, "fusions.out"), os.path.join(output_dir, sample + "_fusions.out") )

class TophatFusion (Component):

    def define_parameters(self, genome_index_base, genome_gtf=None, R1, R2=None, anchor_length=20, nb_threads=1):
        self.add_parameter("genome_index_base", "***************************", default=genome_index_base, required=True)
        self.add_input_file("genome_gtf", "***************************", default=genome_gtf)
        self.add_input_file_list( "R1", "Which R1 files should be used.", default=R1, required=True )
        self.add_input_file_list( "R2", "Which R2 files should be used.", default=R2 )
        self.add_parameter("anchor_length", "****************", type=int, default=anchor_length)
        self.add_parameter("nb_threads", "Number of threads used", type=int, default=nb_threads)
        self.add_output_file_list("aln_files", "The ********* outputed file", pattern='{basename_woext}.bam', items=self.R1)
        self.add_output_file_list("fusions_files", "The ******* outputed file", pattern='{basename_woext}.tsv', items=self.R1)
        self.add_output_file_list("stderr", "The STAR stderr file", pattern='{basename_woext}.stderr', items=self.R1)

    def process(self):
        gtf_opt = '"no_gtf"'
        if self.genome_gtf != None:
            gtf_opt = ' "' + self.genome_gtf + '"'
        tophat = PythonFunction(
            tophatWrapper,
            cmd_format="{EXE} " + self.get_exec_path("tophat2") + ' ' + str(self.nb_threads) + ' ' + str(self.anchor_length) + ' "' + self.genome_index_base + '"' + gtf_opt + " {OUT} {IN}" 
        )
        MultiMap( 
            tophat,
            inputs=( [self.R1, self.R2] if self.R2 != None else [self.R1] ),
            outputs=[self.aln_files, self.stderr],
            includes=( [self.genome_index_base + ".1.ebwt", self.genome_gtf] if self.genome_gtf != None else [self.genome_index_base + ".1.ebwt"])
        )
