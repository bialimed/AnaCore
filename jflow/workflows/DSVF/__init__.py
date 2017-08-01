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
__version__ = '0.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import glob
import warnings

from illumina import SampleSheetIO
from jflow.workflow import Workflow


class DSVF (Workflow):
    """
    @todo: - Add quality check
           - Add filter on NM
    """
    def get_description(self):
        return "Variant analysis on amplicon double strand librairies."

    def define_parameters(self, parameters_section=None):
        # Inputs data
        self.add_input_file( "samplesheet", "Illumina sheet describing run and samples (format: Illumina's samplesheet).", required=True, group="Inputs data" )
        self.add_input_directory( "samples_dir", "Path to the folder containing reads files. [Default: samplesheet directory]", group="Inputs data" )
        # Control samples
        self.add_parameter_list( "pos_ctrl_names", 'Names of the positives controls. These samples contain a known set of variants with a predetermined AF.', group="Control samples" )
        self.add_input_file( "pos_ctrl_expected", 'Variants expected in positives controls (format: VCF).', group="Control samples" )
        # Analysis
        self.add_input_file( "RNA_selection", "The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.", group="Analysis" )
        self.add_parameter( "min_AF", 'Under this allele frequency the variant is tagged "lowAF".', type=float, default=0.02, group="Analysis" )
        self.add_parameter( "min_DP", 'Under this depth in each library the variant is tagged "lowDP". Under this value in one library the variant can be tagged "incomplete" or "invalid".', type=int, default=120, group="Analysis" )
        self.add_input_file( "filters", "Parameters to filter variants (format: JSON).", required=True, group="Analysis" )
        self.add_parameter( "assembly_version", "Genome assembly version.", default="GRCh37", group="Analysis" )
        self.add_input_file( "genome_seq", "Path to the sequences of the genome (format: FASTA).", file_format="fasta", required=True, group="Analysis" )
        # Cleaning
        self.add_input_file( "R1_end_adapter", "Path to sequence file containing the start of Illumina P7 adapter (format: fasta). This sequence is trimmed from the end of R1 of the amplicons with a size lower than read length.", file_format="fasta", required=True, group="Cleaning" )
        self.add_input_file( "R2_end_adapter", "Path to sequence file containing the start of reverse complemented Illumina P5 adapter ((format: fasta). This sequence is trimmed from the end of R2 of the amplicons with a size lower than read length.", file_format="fasta", required=True, group="Cleaning" )
        # Libraries design
        self.add_input_directory( "libA_folder", "Path to folder containing files describing libA design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True, group="Librairies design" )
        self.add_input_directory( "libB_folder", "Path to folder containing files describing libB design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True, group="Librairies design")

    def pre_process(self):
        self.metrics_type = "tsv"

        # Samples
        samples_dir = self.samples_dir if self.samples_dir else os.path.dirname(self.samplesheet)
        samples_dir = samples_dir.replace(" ", "\ ")
        self.samples = SampleSheetIO( self.samplesheet ).samples
        for spl_idx, spl in enumerate(self.samples):
            spl_name = spl["Sample_Name"].replace("_", "-").replace(" ", "-")
            lib_name = spl_name + "_S" + str(spl_idx + 1)
            spl["Manifest"] = "lib" + spl["Manifest"]
            spl["R1"] = sorted( glob.glob(os.path.join(samples_dir, lib_name + '_*_R1_*.fastq.gz')) )
            spl["R2"] = sorted( glob.glob(os.path.join(samples_dir, lib_name + '_*_R2_*.fastq.gz')) )
            if len(spl["R1"]) == 0 or len(spl["R2"]) == 0:
                raise Exception( 'All reads files cannot be found for sample "' + lib_name + '" in folder "' + samples_dir + '".' )

        # InterOp directory
        self.sequencer_run_dir = None
        run_dir = os.path.dirname(os.path.dirname(os.path.dirname(samples_dir)))
        interop_dir = os.path.join(run_dir, "InterOp")
        runParam_file = os.path.join(run_dir, "runParameters.xml")
        runInfo_file = os.path.join(run_dir, "RunInfo.xml")
        if os.path.exists(run_dir) and os.path.exists(interop_dir) and os.path.exists(runParam_file) and os.path.exists(runInfo_file):
            self.sequencer_run_dir = run_dir

        # Libraries designs
        self.libraries = [{
            "name": "libA",
            "design_with_primers": os.path.join( self.libA_folder, "withPrimers.bed" ),
            "design_wout_primers": os.path.join( self.libA_folder, "woutPrimers.bed" ),
            "non_overlap_groups": os.path.join( self.libA_folder, "nonOverlappingGroups.tsv" )
        }, {
            "name": "libB",
            "design_with_primers": os.path.join( self.libB_folder, "withPrimers.bed" ),
            "design_wout_primers": os.path.join( self.libB_folder, "woutPrimers.bed" ),
            "non_overlap_groups": os.path.join( self.libB_folder, "nonOverlappingGroups.tsv" )
        }]

    def process(self):
        R1 = [spl["R1"][0] for spl in self.samples]
        R2 = [spl["R2"][0] for spl in self.samples]

        # Sequencer quality metrics
        if self.sequencer_run_dir != None:
            self.add_component( "InterOp", [self.sequencer_run_dir] )

        # Samples quality control
        self.add_component( "ReadsStat", [(R1 + R2), True] )

        # Prepare reads by amplicon
        cleaned_R1 = self.add_component( "Cutadapt", ["a", self.R1_end_adapter, R1, None, 0.01, 10], component_prefix="R1" )
        cleaned_R2 = self.add_component( "Cutadapt", ["a", self.R2_end_adapter, R2, None, 0.01, 10], component_prefix="R2" )
        
        # Align reads by amplicon
        ################################################################
        for idx in range(len(R1)):
            curr_spl_R1 = os.path.basename( cleaned_R1.out_R1[idx] ).split("_")[0]
            curr_spl_R2 = os.path.basename( cleaned_R2.out_R1[idx] ).split("_")[0]
            if curr_spl_R1 != curr_spl_R2:
                print( curr_spl_R1, curr_spl_R2 )
                raise Exception("Aln list are not consistent")    
        ################################################################
        bwa = self.add_component( "BWAmem", [self.genome_seq, cleaned_R1.out_R1, cleaned_R2.out_R1] )
        idx_aln = self.add_component( "BAMIndex", [bwa.aln_files] )

        for curr_lib in self.libraries:
            curr_lib_aln = [idx_aln.out_aln[idx_spl] for idx_spl, spl in enumerate(self.samples) if spl["Manifest"] == curr_lib["name"]]

            # Add read group by amplicon
            aln_RG = self.add_component( "AddAmpliRG", [curr_lib["design_with_primers"], curr_lib_aln, self.metrics_type], component_prefix=curr_lib["name"] )
            idx_aln_RG = self.add_component( "BAMIndex", [aln_RG.out_aln], component_prefix=curr_lib["name"] )
            curr_lib["aln"] = idx_aln_RG.out_aln
            
            # Coverage
            coverage = self.add_component( "Coverage", [idx_aln_RG.out_aln, curr_lib["design_with_primers"], 1000000], component_prefix=curr_lib["name"] )
            self.add_component( "DepthsDistribution", [curr_lib["design_wout_primers"], coverage.depth_files, self.metrics_type], component_prefix=curr_lib["name"] )
            
            # Variant Calling
            variant_calling = self.add_component( "VarDictAmpli", [self.genome_seq, curr_lib["design_with_primers"], idx_aln_RG.out_aln, self.min_AF], component_prefix=curr_lib["name"] )
            curr_lib["vcf"] = variant_calling.out_variants

        # Merge DS variants
        libA, libB = self.libraries
        libA_spl = [spl["Sample_Name"] for spl in self.samples if spl["Manifest"] == libA["name"]]
        libA["aln"] = [x for (y, x) in sorted(zip(libA_spl, libA["aln"]))]
        libA["vcf"] = [x for (y, x) in sorted(zip(libA_spl, libA["vcf"]))]
        libB_spl = [spl["Sample_Name"] for spl in self.samples if spl["Manifest"] == libB["name"]]
        libB["aln"] = [x for (y, x) in sorted(zip(libB_spl, libB["aln"]))]
        libB["vcf"] = [x for (y, x) in sorted(zip(libB_spl, libB["vcf"]))]
        ################################################################
        for idx in range(len(libA["aln"])):
            curr_spl_libA_aln = os.path.basename( libA["aln"][idx] ).split("_")[0]
            curr_spl_libB_aln = os.path.basename( libB["aln"][idx] ).split("_")[0]
            curr_spl_libA_vcf = os.path.basename( libA["vcf"][idx] ).split("_")[0]
            curr_spl_libB_vcf = os.path.basename( libB["vcf"][idx] ).split("_")[0]
            if curr_spl_libA_aln == curr_spl_libB_aln and curr_spl_libB_aln == curr_spl_libA_vcf and curr_spl_libB_vcf:
                pass
            else:
                print( curr_spl_libA_aln, curr_spl_libB_aln, curr_spl_libA_vcf, curr_spl_libB_vcf)
                raise Exception("Merge list are not consistent")       
        ################################################################
        spl_merging = self.add_component( "MergeVCFAmpli", [libA["design_wout_primers"], libB["design_wout_primers"], libA["vcf"], libB["vcf"], libA["aln"], libB["aln"]] )
        
        # Variants annotation
        variants_annot = self.add_component( "VEP", [spl_merging.out_variants, "homo_sapiens", self.assembly_version] )

        # Filter variants
        tag_annot_variants = self.add_component( "FilterVCFOnAnnot", [variants_annot.out_variants] ) # tags on consequences and polymorphism
        tag_count_variants = self.add_component( "FilterVCFOnCount", [tag_annot_variants.out_variants, self.min_AF, self.min_DP] ) # tags on AF, DP and on presence in two lib
        sorted_variants = self.add_component( "SortVCF", [tag_count_variants.out_variants] )
        final_vcf = sorted_variants.out_variants
        if self.RNA_selection != None:
            select_annot_variants = self.add_component( "FilterVCFAnnotOnRNA", [self.RNA_selection, sorted_variants.out_variants] )
            final_vcf = select_annot_variants.out_variants
        filtered_variants = self.add_component( "FilterVCF", [self.filters, final_vcf] ) # removes variants on their FILTER tags

        # Convert to TSV
        self.add_component( "VCFToTSV", [filtered_variants.out_variants] )

        # Variant calling control (horizon)
        if len(self.pos_ctrl_names) > 0:
            spl_names = sorted([spl["Sample_Name"] for spl in self.samples if spl["Manifest"] == libA["name"]])
            positive_ctrl_variants = [filtered_variants.out_variants[idx] for idx, spl_name in enumerate(spl_names) if spl_name in self.pos_ctrl_names]
            if len(positive_ctrl_variants) > 0:
                self.add_component( "VariantsCtrlCheck", [self.pos_ctrl_expected, positive_ctrl_variants, self.metrics_type] )
            else:
                warnings.warn( 'No positives controls can be found with the names "' + '" or "'.join(self.pos_ctrl_names) + '"' )
