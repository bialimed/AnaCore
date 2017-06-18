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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import glob

from jflow.featureio import SampleSheetIO
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
        self.add_input_file( "samplesheet", "Illumina sheet describing run and samples (format: Illumina's samplesheet).", required=True )
        self.add_input_directory( "samples_dir", "Path to the folder containing reads files. [Default: samplesheet directory]" )
        # Analysis
        self.add_parameter( "min_AF", 'Under this allele frequency the variant is tagged "lowAF".', type=float, default=0.02 )
        self.add_parameter( "min_DP", 'Under this depth in each library the variant is tagged "lowDP". Under this value in one library the variant can be tagged "incomplete" or "invalid".', type=int, default=120 )
        self.add_input_file( "filters", "Parameters to filter variants (format: JSON).", required=True )
        self.add_parameter( "assembly_version", "Genome assembly version.", default="GRCh37" )
        self.add_input_file( "genome_seq", "Path to the sequences of the genome (format: FASTA).", file_format="fasta", required=True )
        # Cleaning
        self.add_input_file( "R1_end_adapter", "Path to sequence file containing the start of Illumina P7 adapter (format: fasta). This sequence is trimmed from the end of R1 of the amplicons with a size lower than read length.", file_format="fasta", required=True )
        self.add_input_file( "R2_end_adapter", "Path to sequence file containing the start of reverse complemented Illumina P5 adapter ((format: fasta). This sequence is trimmed from the end of R2 of the amplicons with a size lower than read length.", file_format="fasta", required=True )
        # Libraries
        self.add_input_directory( "libA_folder", "Path to folder containing files describing libA design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True )
        self.add_input_directory( "libB_folder", "Path to folder containing files describing libB design (regions with primers, regions without primers and groups of non-overlapping amplicons).", required=True )

    def pre_process(self):
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
                raise Esception("Aln list are not consistent")                
        ################################################################
        bwa = self.add_component( "BWAmem", [self.genome_seq, cleaned_R1.out_R1, cleaned_R2.out_R1] )
        idx_aln = self.add_component( "BAMIndex", [bwa.aln_files] )
        #~ self.add_component( "SamtoolsFlagStat", [bwa.aln_files] )

        for curr_lib in self.libraries:
            curr_lib_aln = [idx_aln.out_aln[idx_spl] for idx_spl, spl in enumerate(self.samples) if spl["Manifest"] == curr_lib["name"]]

            # Add read group by amplicon
            aln_RG = self.add_component( "AddAmpliRG", [curr_lib["design_with_primers"], curr_lib_aln], component_prefix=curr_lib["name"] )
            idx_aln_RG = self.add_component( "BAMIndex", [aln_RG.out_aln], component_prefix=curr_lib["name"] )
            curr_lib["aln"] = idx_aln_RG.out_aln
            
            # Coverage
            self.add_component( "Coverage", [idx_aln_RG.out_aln, curr_lib["design_with_primers"]], component_prefix=curr_lib["name"] )
            
            # Variant Calling
            variant_calling = self.add_component( "VarDictAmpli", [self.genome_seq, curr_lib["design_with_primers"], idx_aln_RG.out_aln, self.min_AF], component_prefix=curr_lib["name"] )
            curr_lib["vcf"] = variant_calling.out_variants
        
            #~ # Split BAM in non overlapping regions
            #~ groups_names = set()
            #~ with open(self.non_overlapping_design) as FH_gp:
               #~ for line in FH_gp:
                   #~ if not line.startswith("#"):
                       #~ amplicon_id, group_name = [elt.strip() for elt in line.split("\t")]
                       #~ groups_names.add( group_name )
            #~ groups_names = list(groups_names)
            #~ nb_split_gp = len(groups_names)
            
            #~ splitted_aln = self.add_component( SplitBAMByRG, [non_overlapping_design, idx_aln_RG.out_aln], component_prefix=curr_lib ) ############################

            #~ # Add RG on BAM
            #~ splitted_lib_names = [name for name in curr_lib_names for elt in range(0, nb_split_gp)] 
            #~ splitted_sequencer = ["ILLUMINA" for aln in splitted_aln.out_aln]
            #~ aln_RG = self.add_component( "AddRGOnBAM", [splitted_aln.out_aln, splitted_sequencer, splitted_lib_names, splitted_lib_names], component_prefix=curr_lib  )
            #~ idx_aln_RG = self.add_component( "BAMIndex", [aln_RG.out_aln], component_prefix=curr_lib + "_2" )
            
            #~ # Variant calling
            #~ #echo "freebayes --fasta-reference ${homo_sapiens_genome} ${library}_sorted_RG_grp0_RG.bam > ${library}_grp0.vcf
            #~ variant_calling = self.add_component( "GATKHaplotypeCaller", [idx_aln_RG.out_aln, self.genome_seq, self.genome_SNV, regions_with_primers, 30, 30] )
            #~ ############################################################### Manage empty BAM

            #~ # Remove variants on primers
            #~ zoi_selection = self.add_component( "FilterVCFRegions", [regions_wout_primers, variant_calling.SNV_files], component_prefix=curr_lib )

            #~ # Merge overlapping amplicons
            #~ grouped_zoi = [zoi_selection.out_variants[x:x+nb_split_gp] for elt in range(0, len(zoi_selection.out_variants), nb_split_gp)]
            #~ grouped_aln = [idx_aln_RG.out_aln[x:x+nb_split_gp] for elt in range(0, len(idx_aln_RG.out_aln), nb_split_gp)]
            #~ variants = self.add_component( "MergeOverlappingFragments", [regions_wout_primers, grouped_zoi, grouped_aln], component_prefix=curr_lib ) ############################

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
                raise Esception("Merge list are not consistent")                
        ################################################################
        spl_merging = self.add_component( "MergeVCFAmpli", [libA["design_wout_primers"], libB["design_wout_primers"], libA["vcf"], libB["vcf"], libA["aln"], libB["aln"]] )
        
        # Variants annotation
        variants_annot = self.add_component( "VEP", [spl_merging.out_variants, "homo_sapiens", self.assembly_version] )

        # Filter variants
        tag_annot_variants = self.add_component( "FilterVCFOnAnnot", [variants_annot.out_variants] ) # tags on consequences and polymorphism
        tag_count_variants = self.add_component( "FilterVCFOnCount", [tag_annot_variants.out_variants, self.min_AF, self.min_DP] ) # tags on AF, DP and on presence in two lib
        sorted_variants = self.add_component( "SortVCF", [tag_count_variants.out_variants] )
        filtered_variants = self.add_component( "FilterVCF", [self.filters, sorted_variants.out_variants] ) # removes variants on their FILTER tags

        # Convert to TSV
        self.add_component( "VCFToTSV", [filtered_variants.out_variants] )
