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
import glob
import shutil
import warnings


from illumina import SampleSheetIO
from jflow.workflow import Workflow


class ADIVaR (Workflow):
    def get_description(self):
        return "Variant analysis on amplicon double strand librairies."

    def get_cmpt_by_nameid( self, selected_nameid ):
        selected_cpmt = None
        for cmpt in self.components:
            if cmpt.get_nameid() == selected_nameid:
                selected_cpmt = cmpt
        return selected_cpmt

    def define_parameters(self, parameters_section=None):
        # Inputs data
        self.add_input_file( "samplesheet", "Illumina sheet describing run and samples (format: Illumina's samplesheet).", required=True, group="Inputs data" )
        self.add_input_directory( "samples_dir", "Path to the folder containing reads files. [Default: samplesheet directory]", group="Inputs data" )
        # Output data
        self.add_parameter( "metrics_type", 'Type of metrics files: "tsv" for human readable or "json" for machine readable.', choices=["json", "tsv"], default="json", group="Output data" )
        self.add_input_directory( "output_dir", "Path to the output folder.", required=True, group="Output data" )
        # Control samples
        self.add_parameter_list( "pos_ctrl_names", 'Names of the positives controls. These samples contain a known set of variants with a predetermined AF.', group="Control samples" )
        self.add_input_file( "pos_ctrl_expected", 'Variants expected in positives controls (format: VCF).', group="Control samples" )
        # Variant calling
        self.add_parameter( "min_call_AF", 'Variants with an allele frequency under this value are not emitted by variant caller.', type=float, default=0.02, group="Variant calling" )
        self.add_parameter( "min_base_qual", 'The phred score for a base to be considered a good call.', type=int, default=25, group="Variant calling" )
        # Variant filters
        self.add_parameter_list( "kept_consequences", "The variants without one of these consequences are removed (see http://www.ensembl.org/info/genome/variation/predicted_data.html).",
            default=["TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", "splice_region_variant"],
            group="Variants filters" )
        self.add_parameter( "min_AF", 'Under this allele frequency the variant is tagged "lowAF".', type=float, default=0.015, group="Variants filters" )
        self.add_parameter( "min_DP", 'Under this depth in each library the variant is tagged "lowDP". Under this value in one library the variant can be tagged "incomplete" or "invalid".', type=int, default=100, group="Variants filters" )
        self.add_input_file( "filters", "Parameters to filter variants (format: JSON).", required=True, group="Variants filters" )
        # Analysis
        self.add_input_file( "RNA_selection", "The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.", group="Analysis" )
        self.add_parameter( "assembly_version", "Genome assembly version.", default="GRCh37", group="Analysis" )
        self.add_input_file( "genome_seq", "Path to the sequences of the genome (format: FASTA).", file_format="fasta", required=True, group="Analysis" )
        # Cleaning
        self.add_input_file( "R1_end_adapter", "Path to sequence file containing the start of Illumina P7 adapter (format: fasta). This sequence is trimmed from the end of R1 of the amplicons with a size lower than read length.", file_format="fasta", required=True, group="Cleaning" )
        self.add_input_file( "R2_end_adapter", "Path to sequence file containing the start of reverse complemented Illumina P5 adapter ((format: fasta). This sequence is trimmed from the end of R2 of the amplicons with a size lower than read length.", file_format="fasta", required=True, group="Cleaning" )
        # Libraries design
        self.add_input_directory( "libA_folder", "Path to folder containing files describing libA design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True, group="Librairies design" )
        self.add_input_directory( "libB_folder", "Path to folder containing files describing libB design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True, group="Librairies design")

    def pre_process(self):
        # Samples
        samples_dir = self.samples_dir if self.samples_dir else os.path.dirname(self.samplesheet)
        samples_dir = samples_dir.replace(" ", "\ ")
        self.pos_ctrl_spl = set()
        self.samples = SampleSheetIO( self.samplesheet ).samples
        for spl_idx, spl in enumerate(self.samples):
            if spl["Sample_Name"] in self.pos_ctrl_names:
                self.pos_ctrl_spl.add( spl["Sample_Name"] )
            spl_name = spl["Sample_Name"].replace("_", "-").replace(" ", "-")
            lib_name = spl_name + "_S" + str(spl_idx + 1)
            spl["Manifest"] = "lib" + spl["Manifest"]
            spl["R1"] = sorted( glob.glob(os.path.join(samples_dir, lib_name + '_*_R1_*.fastq.gz')) )
            spl["R2"] = sorted( glob.glob(os.path.join(samples_dir, lib_name + '_*_R2_*.fastq.gz')) )
            if len(spl["R1"]) == 0 or len(spl["R2"]) == 0:
                raise Exception( 'All reads files cannot be found for sample "' + lib_name + '" in folder "' + samples_dir + '".' )
        if len(self.pos_ctrl_names) > 0 and len(self.pos_ctrl_spl) == 0:
            warnings.warn( 'No positives controls can be found with the names "' + '" or "'.join(self.pos_ctrl_names) + '"' )

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
            variant_calling = self.add_component( "AmpliVariantCalling", [self.genome_seq, curr_lib["design_with_primers"], curr_lib["design_wout_primers"], curr_lib["non_overlap_groups"], idx_aln_RG.out_aln, self.min_call_AF, self.min_base_qual], component_prefix=curr_lib["name"] )
            curr_lib["vcf"] = variant_calling.out_variants

        # Merge DS variants
        libA, libB = self.libraries
        libA_spl = [spl["Sample_Name"] for spl in self.samples if spl["Manifest"] == libA["name"]]
        libA["aln"] = [x for (y, x) in sorted(zip(libA_spl, libA["aln"]))]
        libA["vcf"] = [x for (y, x) in sorted(zip(libA_spl, libA["vcf"]))]
        libB_spl = [spl["Sample_Name"] for spl in self.samples if spl["Manifest"] == libB["name"]]
        libB["aln"] = [x for (y, x) in sorted(zip(libB_spl, libB["aln"]))]
        libB["vcf"] = [x for (y, x) in sorted(zip(libB_spl, libB["vcf"]))]
        spl_merging = self.add_component( "MergeVCFAmpli", [libA["design_wout_primers"], libB["design_wout_primers"], libA["vcf"], libB["vcf"], libA["aln"], libB["aln"]] )

        # Variants annotation
        variants_annot = self.add_component( "VEP", [spl_merging.out_variants, "homo_sapiens", self.assembly_version] )

        # Filter variants
        tag_annot_variants = self.add_component( "FilterVCFOnAnnot", [variants_annot.out_variants, self.kept_consequences] ) # tags on consequences and polymorphism
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

    def post_process(self):
        # Move alignment files
        for curr_lib in self.libraries:
            aln_lib_folder = os.path.join( self.output_dir, "alignments_" + curr_lib["name"] )
            if not os.path.exists(aln_lib_folder): os.mkdir(aln_lib_folder)
            for aln in self.get_cmpt_by_nameid("AddAmpliRG." + curr_lib["name"]).out_aln: # Aln
                filename = os.path.basename(aln).split("_")[0] + ".bam"
                shutil.move( aln, os.path.join(aln_lib_folder, filename) )
            for bai in self.get_cmpt_by_nameid("BAMIndex." + curr_lib["name"]).out_idx: # Aln idx
                filename = os.path.basename(bai).split("_")[0] + ".bam.bai"
                shutil.move( bai, os.path.join(aln_lib_folder, filename) )

        # Move variants files
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

        # Move metrics files
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
