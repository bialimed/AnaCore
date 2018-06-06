#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class FilterVCFNoise (Component):

    def define_parameters(self, in_variants, in_noises, mode="tag", tag_name="popConst", tag_desc="The variant correspond to a constitutive detection (sequencing polymerase error rate at this position, workflow artifact, ...)."):
        # Parameters
        self.add_parameter("mode", 'Select the filter mode. In mode "tag" if the variant is noise in all samples present in the VCF the a tag is added in FILTER (cf. "--tag-name"). In mode "remove" if the variant is noise all samples present in the VCF this variant is removed.', choices=["tag", "remove"], default=mode)
        self.add_parameter("tag_desc", "The description of the applied filter. This description is stored in header of the VCF.", default=tag_desc)
        self.add_parameter("tag_name", "The name of the tag added on variant if it correspond to noise.", default=tag_name)

        # Files
        self.add_input_file("in_noises", 'The path to the file containing artifactual variants with their maximum frequency (format: TSV). The header line of the file must be "#Chromosome<tab>Possition<tab>Reference_allele<tab>Alternative_allele<tab>Noise_rate".', default=in_noises, required=True)
        self.add_input_file_list("in_variants", "The path to the variants file (format: VCF).", default=in_variants, required=True)
        self.add_output_file_list("out_variants", "Path to the filtered files (format: VCF).", pattern='{basename_woext}.vcf', items=self.in_variants)
        self.add_output_file_list("stderr", "Path to the stderr files (format: txt).", pattern='{basename_woext}.stderr', items=self.in_variants)

    def process(self):
        cmd = self.get_exec_path("filterVCFNoise.py") + \
            " --mode " + self.mode + \
            " --tag-name '{}'".format(self.tag_name) + \
            " --tag-description '{}'".format(self.tag_desc) + \
            " --input-noises " + self.in_noises + \
            " --input-variants $1" + \
            " --output-variants $2" + \
            " 2> $3"
        filter_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(filter_fct, inputs=[self.in_variants], outputs=[self.out_variants, self.stderr], includes=[self.in_noises])
