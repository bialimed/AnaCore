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


class MergeVCFAmpli (Component):

    def define_parameters(self, libA_regions_wout_primers, libB_regions_wout_primers, libA_variants, libB_variants, libA_aln, libB_aln):
        self.add_input_file( "libA_regions_wout_primers", "***** (format: BED)", default=libA_regions_wout_primers, required=True )
        self.add_input_file( "libB_regions_wout_primers", "***** (format: BED)", default=libB_regions_wout_primers, required=True )    
        self.add_input_file_list( "libA_variants", "*****", default=libA_variants, required=True )
        self.add_input_file_list( "libB_variants", "*****", default=libB_variants, required=True )
        self.add_input_file_list( "libA_aln", "*****", default=libA_aln, required=True )
        self.add_input_file_list( "libB_aln", "*****", default=libB_aln, required=True )
        self.add_output_file_list( "out_variants", "******************", pattern='{basename_woext}.vcf', items=self.libA_variants )
        self.add_output_file_list( "stderr", "The stderr file", pattern='{basename_woext}.stderr', items=self.libA_variants )

    def process(self):
        cmd = self.get_exec_path("mergeVCFAmpli.py") + \
            " --input-designs " + self.libA_regions_wout_primers + " " + self.libB_regions_wout_primers + \
            " --input-variants $1 $2" + \
            " --input-aln $3 $4" + \
            " --output-variants $5" + \
            " 2> $6"
        merge_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap(
            merge_fct,
            inputs=[self.libA_variants, self.libB_variants, self.libA_aln, self.libB_aln],
            outputs=[self.out_variants, self.stderr],
            includes=[self.libA_regions_wout_primers, self.libB_regions_wout_primers]
        )
