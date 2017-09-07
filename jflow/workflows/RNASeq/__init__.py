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
__copyright__ = 'Copyright (C) 2017 IUCT'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import os
import sys

from jflow.workflow import Workflow
from jflow.parameter import InputFileList, InputFile


class RNASeq (Workflow):
    def get_description(self):
        return "RNASeq exploration"

    def define_parameters(self, parameters_section=None):
        self.add_input_file_list( "R1", "The list of R1 files", file_format="fastq", required=True )
        self.add_input_file_list( "R2", "The list of R2 files", file_format="fastq", required=True )
        # Reference
        self.add_parameter( "ref_species", "Reference species used.", default="homo_sapiens", group="Reference" )
        self.add_parameter( "ref_assemby", "Assembly version.", default="GRCh37", group="Reference" )
        # Reference genome
        self.add_parameter( "genome_dir", "Path to the directory containing the STAR index of the human genome.", required=True, group="Genome" )
        self.add_input_file( "genome_annot", "Path to the annotation of the genome (format: GTF).", required=True, group="Genome" )
        self.add_input_file( "genome_seq", "Path to the sequences of the genome (format: FASTA).", required=True, group="Genome" )
        self.add_input_file( "genome_SNV", "Path to the variants of the genome (format: VCF).", required=True, group="Genome" )
        # Reference transcriptome
        self.add_parameter( "transcriptome_idx", "Path to the salmon index for the human transcriptome (format: salmon index).", group="Transcrptome" )

    def process(self):
        # Samples quality control
        self.add_component( "ReadsStat", [(self.R1 + self.R2), True] )

        # Alignment
        star = self.add_component( "STAR", [self.genome_dir, self.genome_annot, self.R1, self.R2] )
        # self.add_component( "SamtoolsFlagStat", [star.aln_files] ) # Problem not aligned are removed from BAM ###############################
        # self.add_component( "CigarlineGraph", [star.aln_files] ) # Require samtools calmd ##############################

        # SNV detection
        gatk_prepro = self.add_component( "GATKPreprocess", [star.aln_files, self.genome_seq, True] )
        gatk_recalib = self.add_component( "GATKRecalibration", [gatk_prepro.cleaned_files, self.genome_seq, self.genome_SNV] )
        variant_calling = self.add_component( "GATKHaplotypeCaller", [gatk_recalib.output_files, self.genome_seq, self.genome_SNV, None, 20] )
        variant_filtering = self.add_component( "GATKFilter", [variant_calling.SNV_files, self.genome_seq, None, False] )
        variant_annotation = self.add_component( "VEP", [variant_filtering.marked_SNV_files, self.ref_species, self.ref_assembly] )

        # Quantification
        if self.transcriptome_idx != None:
            salmon = self.add_component( "Salmon", [self.transcriptome_idx, self.R1, self.R2, "IU"] )
