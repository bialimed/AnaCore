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
__version__ = '2.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import PythonFunction


def VEPWrapper(exec_path, in_vcf, out_file, stderr, out_format, species, assembly, cache_directory):
    """Wrap VEP execution for prevent error when the in_vcf does not contain any records."""
    import subprocess

    # Check if VCF contains records
    is_empty = True
    with open(in_vcf) as FH_in:
        for line in FH_in:
            if not line.startswith("#"):
                is_empty = False
    # Process
    if is_empty:  # Create empty output if VCF does not contain any records
        with open(stderr, "w") as FH_err:
            FH_err.write('[WARN] The input file "' + stderr + '" does not contain any records.')
        if out_format == "json":
            with open(out_file, "w") as FH_out:
                FH_out.write("{}")
        elif out_format == "tab":
            with open(out_file, "w") as FH_out:
                FH_out.write("")
        else:  # vcf
            with open(in_vcf) as FH_in:
                with open(out_file, "w") as FH_out:
                    for line in FH_in:
                        FH_out.write(line)
    else:  # Annotate
        cmd = exec_path + \
              " --cache" + \
              ("" if cache_directory == "none" else " --dir_cache " + cache_directory) + \
              " --offline" + \
              ' --species "' + species + '"' + \
              ("" if assembly == "none" else ' --assembly "' + assembly + '"') + \
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
              " --vcf_info_field ANN " + \
              " --" + out_format + \
              " --input_file " + in_vcf + \
              " --output_file " + out_file + \
              " 2> " + stderr
        subprocess.check_output(cmd, shell=True)


class VEP (Component):

    def define_parameters(self, in_variants, species, assembly=None, out_format="vcf"):
        # Parameters
        self.add_parameter("assembly", "The assembly version to use.", default=assembly)
        self.add_parameter("out_format", "The output format.", choices=["json", "vcf", "tab"], default=out_format)
        self.add_parameter("species", 'The latin name for your species (example: "homo_sapiens")', default=species, required=True)

        # Input files
        self.add_input_file_list("in_variants", "Path to the variants file before annotation (format: VCF).", default=in_variants, required=True)

        # Output files
        self.add_output_file_list("out_variants", "The VEP annotations files (format: see out_format).", pattern='{basename_woext}.' + out_format, items=self.in_variants)
        self.add_output_file_list("stderr", "Path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.in_variants)

    def process(self):
        # Transform parameters for command line
        cache_param = "none"
        try:
            cache_param = self.get_resource("VEP_cache_directory")
        except:
            pass
        assembly_param = "none" if self.assembly == None else self.assembly
        # Command
        vep_fct = PythonFunction(
            VEPWrapper,
            cmd_format="{EXE}" + \
                       " " + self.get_exec_path("vep") + \
                       " {IN}" + \
                       " {OUT}" + \
                       " " + self.out_format + \
                       ' "' + self.species + '"' + \
                       ' "' + assembly_param + '"' + \
                       ' "' + cache_param + '"'
        )
        MultiMap(vep_fct, inputs=self.in_variants, outputs=[self.out_variants, self.stderr])
