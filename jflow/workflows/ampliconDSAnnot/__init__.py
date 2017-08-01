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
import shutil
import warnings

from illumina import SampleSheetIO
from jflow.workflow import Workflow


class AmpliconDSAnnot (Workflow):

    def get_description(self):
        return "Process annotations on variants produced by AmpliconDS."

    def define_parameters(self, parameters_section=None):
        # Inputs data
        self.add_input_file( "samplesheet", "Illumina sheet describing run and samples (format: Illumina's samplesheet).", required=True, group="Inputs data" )
        self.add_input_directory( "variants_dir", "Path to the folder containing variants files. [Default: Alignment directory in samplesheet directory]", group="Inputs data" )
        # Control samples
        self.add_parameter_list( "pos_ctrl_names", 'Names of the positives controls. These samples contain a known set of variants with a predetermined AF.', group="Control samples" )
        self.add_input_file( "pos_ctrl_expected", 'Variants expected in positives controls (format: VCF).', group="Control samples" )
        # Analysis
        self.add_input_file( "RNA_selection", "The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.", group="Analysis" )
        self.add_input_file( "filters", "Parameters to filter variants (format: JSON).", required=True, group="Analysis" )
        self.add_parameter( "assembly_version", "Genome assembly version.", default="GRCh37", group="Analysis" )

    def pre_process(self):
        variants_dir = self.variants_dir if self.variants_dir else os.path.join(os.path.dirname(self.samplesheet), "Alignment")
        self.variants_dir = variants_dir.replace(" ", "\ ")
        self.samples_names = list()
        self.in_variants = list()
        already_processed = dict()
        for spl in SampleSheetIO( self.samplesheet ).samples:
            spl_name = spl["Sample_Name"].replace("_", "-").replace(" ", "-")
            if spl_name not in already_processed:
                self.samples_names.append( spl_name )
                self.in_variants.append( os.path.join(variants_dir, spl_name + '.vcf') )
                if not os.path.exists(self.in_variants[-1]):
                    raise Exception( 'The variant file cannot be found for sample "' + spl_name + '" in folder "' + variants_dir + '".' )
            already_processed[spl_name] = 1

    def process(self):
        # Variants annotation
        variants_annot = self.add_component( "VEP", [self.in_variants, "homo_sapiens", self.assembly_version] )

        # Filter variants
        tag_annot_variants = self.add_component( "FilterVCFOnAnnot", [variants_annot.out_variants] ) # tags on consequences and polymorphism
        sorted_variants = self.add_component( "SortVCF", [tag_annot_variants.out_variants] )
        final_vcf = sorted_variants.out_variants
        if self.RNA_selection != None:
            select_annot_variants = self.add_component( "FilterVCFAnnotOnRNA", [self.RNA_selection, sorted_variants.out_variants] )
            final_vcf = select_annot_variants.out_variants
        filtered_variants = self.add_component( "FilterVCF", [self.filters, final_vcf] ) # removes variants on their FILTER tags

        # Convert to TSV
        self.add_component( "VCFToTSV", [filtered_variants.out_variants] )

        # Variant calling control (horizon)
        if len(self.pos_ctrl_names) > 0:
            positive_ctrl_variants = [filtered_variants.out_variants[idx] for idx, spl_name in enumerate(self.samples_names) if spl_name in self.pos_ctrl_names]
            if len(positive_ctrl_variants) > 0:
                self.add_component( "VariantsCtrlCheck", [self.pos_ctrl_expected, positive_ctrl_variants] )
            else:
                warnings.warn( 'No positives controls can be found with the names "' + '" or "'.join(self.pos_ctrl_names) + '"' )

    def get_cmpt_by_nameid( self, selected_nameid ):
        selected_cpmt = None
        for cmpt in self.components:
            if cmpt.get_nameid() == selected_nameid:
                selected_cpmt = cmpt
        return selected_cpmt

    def post_process(self):
        # TSV
        for tsv in self.get_cmpt_by_nameid("VCFToTSV.default").out_variants:
            shutil.move( tsv, self.variants_dir )
        # Positives control
        for ctrl in self.get_cmpt_by_nameid("VariantsCtrlCheck.default").eval_files:
            shutil.move( ctrl, self.variants_dir )
        # Annotation comple
        log_file = os.path.join(self.variants_dir, "annotationComplete.txt")
        with open(log_file, "w") as FH_log:
            FH_log.write("success")
