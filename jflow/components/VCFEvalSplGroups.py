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

import os
from jflow.component import Component
from weaver.function import ShellFunction


class VCFEvalSplGroups (Component):
    def define_parameters(self, variants_files, groups_file, minimum_AF=0.0, distance_method="euclidean", linkage_method="ward", out_format="human"):
        # Parameters
        self.add_parameter("distance_method", "Used distance (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html#module-scipy.spatial.distance).", default=distance_method, required=True)
        self.add_parameter("linkage_method", "Used linkage (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage).", default=linkage_method, required=True)
        self.add_parameter("minimum_AF", "Only variants with at least one sample with an allele frequency superior or equal to this value are evaluated.", default=minimum_AF, required=True)
        self.add_parameter("out_format", "The output format. With human invalids are in TSV and tree in png ; With computer invalids and computer are in JSON.", choices=["human", "computer"], default=out_format, required=True)

        # Input Files
        self.add_input_file("groups_file", "Path to the file describing links between samples and groups (format: TSV).", default=groups_file, required=True)
        self.add_input_file_list("variants_files", "The path to the variants files (format: VCF).", default=variants_files, required=True)

        # Output Files
        self.add_output_file(
            "out_invalids",
            "Path to the file listing the samples more similar with samples from another group than others samples coming from the same group. (format: see output_format)",
            filename='invalids_groups.{}'.format("tsv" if self.out_format == "human" else "json")
        )
        self.add_output_file(
            "out_tree",
            "Path to the hierarchical clustering tree (format: see output_format).",
            filename='samples_HC.{}'.format("png" if self.out_format == "human" else "json")
        )
        self.add_output_file("stderr", "Path to the stderr files (format: txt).", filename='groups_eval.stderr')

    def process(self):
        # Write list of inputs variants to prevent bug xith size limitation of command line.
        vcf_list = os.path.join(self.output_directory, "variants_files_list.txt")
        with open(vcf_list, "w") as FH_list:
            FH_list.write(
                " ".join([elt.replace(" ", "\ ") for elt in self.variants_files])
            )
        # Submit command
        cmd = self.get_exec_path("VCFEvalSplGroups.py") + \
            " --min-AF " + str(self.minimum_AF) + \
            " --distance-method " + self.distance_method + \
            " --linkage-method " + self.linkage_method + \
            " --output-format " + self.out_format + \
            " --input-groups $1" + \
            " --inputs-variants `cat $2`" + \
            " --output-invalids $3" + \
            " --output-tree $4" + \
            " 2> $5"
        eval_groups_fct = ShellFunction(cmd, cmd_format='{EXE} {IN} {OUT}')
        eval_groups_fct(inputs=[self.groups_file, vcf_list], outputs=[self.out_invalids, self.out_tree, self.stderr], includes=[self.variants_files])
