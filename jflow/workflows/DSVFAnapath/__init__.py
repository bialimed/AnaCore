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
__version__ = '0.3.0'
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
        self.add_input_directory( "output_dir", "Path to the output folder.", group="Output data" )

    def process(self):
        # Store in database
        DSVF.process(self)
        ############################# Convert BAM to CRAM

    def get_cmpt_by_nameid( self, selected_nameid ):
        selected_cpmt = None
        for cmpt in self.components:
            if cmpt.get_nameid() == selected_nameid:
                selected_cpmt = cmpt
        return selected_cpmt

    def post_process(self):
        DSVF.post_process(self)

        # Interop
        if self.sequencer_run_dir != None:
            interop_folder = os.path.join( self.output_dir, "instrument_log" )
            os.mkdir(interop_folder)
            for filename in ["RunInfo.xml", "runParameters.xml"]:
                shutil.copyfile( os.path.join(self.sequencer_run_dir, filename), os.path.join(interop_folder, filename) )
            shutil.copytree( os.path.join(self.sequencer_run_dir, "InterOp"), os.path.join(interop_folder, "InterOp") )

        # Reads stat
        ############################# ReadsStat

        # Alignment
        for curr_lib in self.libraries:
            aln_lib_folder = os.path.join( self.output_dir, "alignments_" + curr_lib["name"] )
            os.mkdir(aln_lib_folder)
            for aln in self.get_cmpt_by_nameid("AddAmpliRG." + curr_lib["name"]).out_aln: # Aln
                filename = os.path.basename(aln).split("_")[0] + ".bam"
                shutil.move( aln, os.path.join(aln_lib_folder, filename) )
            for bai in self.get_cmpt_by_nameid("BAMIndex." + curr_lib["name"]).out_idx: # Aln idx
                filename = os.path.basename(bai).split("_")[0] + ".bam.bai"
                shutil.move( bai, os.path.join(aln_lib_folder, filename) )

        # Variants
        variants_folder = os.path.join( self.output_dir, "variants" )
        os.mkdir(variants_folder)
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
        os.mkdir(data_folder)
        for curr_lib in self.libraries:
            # Alignment statistics
            for metrics in self.get_cmpt_by_nameid("AddAmpliRG." + curr_lib["name"]).out_summary:
                filename = os.path.basename(metrics).split("_")[0] + "_" + curr_lib["name"] + "_aln.tsv"
                shutil.move( metrics, os.path.join(data_folder, filename) )
            # Depths
            for depths in self.get_cmpt_by_nameid("DepthsDistribution." + curr_lib["name"]).out_depths:
                filename = os.path.basename(depths).split("_")[0] + "_" + curr_lib["name"] + "_depths.json"
                shutil.move( depths, os.path.join(data_folder, filename) )
        # Variants
        ########################### nothing
        # Positives control
        for ctrl in self.get_cmpt_by_nameid("VariantsCtrlCheck.default").eval_files:
            filename = os.path.basename(ctrl).split("_")[0] + "_ctrl.tsv"
            shutil.move( ctrl, os.path.join(data_folder, filename) )
        # Sequencer metrics
        if self.sequencer_run_dir != None:
            interop = self.get_cmpt_by_nameid("InterOp.default").summary_file
            filename = os.path.basename(interop)
            shutil.move( interop, os.path.join(data_folder, filename) )
