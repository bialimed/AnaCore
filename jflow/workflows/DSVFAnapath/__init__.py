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
__version__ = '0.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import shutil

from jflow.workflow import Workflow
from workflows.DSVF import DSVF

class DSVFAnapath (DSVF):
    def get_description(self):
        return "Variant analysis on amplicon double strand librairies."

    def define_parameters(self, parameters_section=None):
        DSVF.define_parameters(self)
        self.add_input_directory( "output_dir", "Path to the output folder.", required=True, group="Output data" )

    def pre_process(self):
        DSVF.pre_process(self)
        self.metrics_type = "json"

    def process(self):
        DSVF.process(self)

    def get_cmpt_by_nameid( self, selected_nameid ):
        selected_cpmt = None
        for cmpt in self.components:
            if cmpt.get_nameid() == selected_nameid:
                selected_cpmt = cmpt
        return selected_cpmt

    def post_process(self):
        DSVF.post_process(self)

        # Alignment
        for curr_lib in self.libraries:
            aln_lib_folder = os.path.join( self.output_dir, "alignments_" + curr_lib["name"] )
            if not os.path.exists(aln_lib_folder): os.mkdir(aln_lib_folder)
            for aln in self.get_cmpt_by_nameid("AddAmpliRG." + curr_lib["name"]).out_aln: # Aln
                filename = os.path.basename(aln).split("_")[0] + ".bam"
                shutil.move( aln, os.path.join(aln_lib_folder, filename) )
            for bai in self.get_cmpt_by_nameid("BAMIndex." + curr_lib["name"]).out_idx: # Aln idx
                filename = os.path.basename(bai).split("_")[0] + ".bam.bai"
                shutil.move( bai, os.path.join(aln_lib_folder, filename) )

        # Variants
        variants_folder = os.path.join( self.output_dir, "variants" )
        if not os.path.exists(variants_folder): os.mkdir(variants_folder)
        for vcf in self.get_cmpt_by_nameid("SortVCF.default").out_variants: # VCF
            filename = os.path.basename(vcf).split("_")[0] + ".vcf"
            shutil.move( vcf, os.path.join(variants_folder, filename) )
        for vcf_filtered in self.get_cmpt_by_nameid("FilterVCF.default").out_variants: # VCF filtered
            filename = os.path.basename(vcf_filtered).split("_")[0] + "_filtered.vcf"
            shutil.move( vcf_filtered, os.path.join(variants_folder, filename) )
        for tsv in self.get_cmpt_by_nameid("VCFToTSV.default").out_variants: # TSV
            filename = os.path.basename(tsv).split("_")[0] + "_filtered.tsv"
            shutil.move( tsv, os.path.join(variants_folder, filename) )

        # Data
        data_folder = os.path.join( self.output_dir, "data" )
        if not os.path.exists(data_folder): os.mkdir(data_folder)
        for curr_lib in self.libraries:
            # Reads quality
            #################################################################### TODO
            # Alignment statistics
            for metrics in self.get_cmpt_by_nameid("AddAmpliRG." + curr_lib["name"]).out_summary:
                filename = os.path.basename(metrics).split("_")[0] + "_" + curr_lib["name"] + "_aln.json"
                shutil.move( metrics, os.path.join(data_folder, filename) )
            # Depths
            for depths in self.get_cmpt_by_nameid("DepthsDistribution." + curr_lib["name"]).out_depths:
                filename = os.path.basename(depths).split("_")[0] + "_" + curr_lib["name"] + "_depths.json"
                shutil.move( depths, os.path.join(data_folder, filename) )
        # Variants
        ######################################################################## TODO ?
        # Positives control
        if len(self.pos_ctrl_spl) > 0:
            for ctrl in self.get_cmpt_by_nameid("VariantsCtrlCheck.default").eval_files:
                filename = os.path.basename(ctrl).split("_")[0] + "_ctrl.json"
                shutil.move( ctrl, os.path.join(data_folder, filename) )
        # Sequencer metrics
        if self.sequencer_run_dir != None:
            interop = self.get_cmpt_by_nameid("InterOp.default").summary_file
            filename = os.path.basename(interop)
            shutil.move( interop, os.path.join(data_folder, filename) )
