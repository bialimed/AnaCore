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

from weaver.function import ShellFunction


class SplitBAMByRG (Component):

    def define_parameters(self, in_regions, in_aln):
        self.add_input_file( "in_regions", "*****", default=in_regions, required=True )    
        self.add_input_file_list( "in_aln", "*****", default=aln_files, required=True )


        groups_names = set()
        with open(in_amplicon_design) as FH_gp:
           for line in FH_gp:
               if not line.startswith("#"):
                   amplicon_id, group_name = [elt.strip() for elt in line.split("\t")]
                   groups_names.add( group_name )
        out_items = list()
        for curr_aln in in_aln:
			curr_aln_basename = os.path.basename( curr_aln )
			
        self.add_output_file_list( "out_aln", "******************", pattern='{basename_woext}.bam', items=self.in_aln )
        self.add_output_file_list( "stdout", "******************", pattern='{basename_woext}.tsv', items=self.in_aln )
        self.add_output_file_list( "stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.in_aln )

    def process(self):
        cmd = self.get_exec_path("splitBAMByRG.py") + \
            " --input-design " + self.in_regions + \
            " --input-aln $1" + \
            " --output-pattern " + out_pattern + \
            " 2> $2"
        split = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap(addRG, inputs=[self.in_aln], outputs=[self.stderr], includes=[self.in_regions] )
