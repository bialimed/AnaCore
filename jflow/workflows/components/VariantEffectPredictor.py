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

import os

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class VEP (Component):
    
    def define_parameters(self, in_variants, species, assembly=None, out_format="vcf"):
        # Parameters
        self.add_parameter( "assembly", "The assembly version to use.", default=assembly )
        self.add_parameter( "out_format", "The output format.", choices=["json", "vcf", "tab"], default=out_format )
        self.add_parameter( "species", 'The latin name for your species (example: "homo_sapiens")', default=species, required=True  )

        # Input files
        self.add_input_file_list( "in_variants", "Path to the variants file before annotation (format: VCF).", default=in_variants, required=True )

        # Output files
        self.add_output_file_list( "out_variants", "The VEP annotations files (format: see out_format).", pattern='{basename_woext}.' + out_format, items=self.in_variants )
        self.add_output_file_list( "stderr", "The error files.", pattern='{basename_woext}.stderr', items=self.in_variants )
        
    def process(self):
        cache_directory = None
        try: 
            cache_directory = self.get_resource("VEP_cache_directory")
        except: pass
        cmd = self.get_exec_path("vep") + \
              " --cache" + \
              ("" if cache_directory == None else " --dir_cache " + cache_directory) + \
              " --offline" + \
              " --species " + self.species + \
              ("" if self.assembly == None else " --assembly " + self.assembly) + \
              " --sift b" + \
              " --polyphen b" + \
              " --gene_phenotype" + \
              " --numbers" + \
              " --hgvs" + \
              " --protein" + \
              " --symbol" + \
              " --uniprot" + \
              " --tsl" + \
              " --biotype" + \
              " --xref_refseq" + \
              " --af" + \
              " --af_1kg" + \
              " --af_esp" + \
              " --af_exac" + \
              " --pubmed" + \
              " --no_stats" + \
              " --merged" + \
              " --" + self.out_format + \
              " --input_file $1" + \
              " --output_file $2" + \
              " 2> $3"
        vep = ShellFunction( cmd, cmd_format="{EXE} {IN} {OUT}")
        MultiMap(vep, inputs=self.in_variants, outputs=[self.out_variants, self.stderr])
