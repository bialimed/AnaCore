#
# Copyright (C) 2012 INRA
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
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class GATKRecalibration (Component):
    """
    @summary : Reprocess quality score for each nucleotids in reads.
    """
    def define_parameters(self, input_files, reference, vcf=None, memory=10):
        # Parameters
        self.add_parameter( "memory", "The memory size to allocate (in GigaBytes).", type="int", default=memory )

        # Files
        self.add_input_file_list( "input_files", "The BAM files.", default=input_files, file_format="bam", required=True )
        self.add_input_file( "reference", "The reference file.", default=reference, file_format="fasta", required=True )
        self.add_input_file( "known_sites", "The VCF of known sites.", default=vcf, file_format="vcf" )

        self.add_output_file_list( "recab_files", "Recalibration table files. Content : several covariate values, num observations, num mismatches and empirical quality score.", pattern="{basename_woext}_recalib.grp", items=input_files )
        self.add_output_file_list( "recab_stderrs", "Error file for recalibration calculation step.", pattern="{basename_woext}_recalib.stderr", items=input_files )

        self.add_output_file_list( "output_files", "Recalibrated bams.", pattern="{basename_woext}.bam", items=input_files, file_format="bam" )
        self.add_output_file_list( "print_reads_stderrs", "Error file for recalibration step.", pattern="{basename_woext}_print_reads.stderr", items=input_files )

    def process(self):
        xmx_option = "-Xmx" + str(self.memory) + "g"

        # Recalibration step 1: base recalibration
        dbsnp_options =  " --run_without_dbsnp_potentially_ruining_quality " if self.known_sites == None else " -knownSites " + self.known_sites
        baserecalib = ShellFunction( self.get_exec_path("java") + " " + xmx_option + " -jar " + self.get_exec_path("gatk") + " -T BaseRecalibrator -S SILENT --allow_potentially_misencoded_quality_scores"\
                                     " -I $1 -R " + self.reference + dbsnp_options + " -o $2 2>> $3",  cmd_format='{EXE} {IN} {OUT}' )
        if self.known_sites != None:
            baserecalib_map = MultiMap( baserecalib, inputs=[self.input_files], outputs=[self.recab_files , self.recab_stderrs], includes=[self.reference,self.known_sites] )
        else:
            baserecalib_map = MultiMap( baserecalib, inputs=[self.input_files], outputs=[self.recab_files , self.recab_stderrs], includes=self.reference )

        # Recalibration step 2: print reads
        print_reads = ShellFunction( self.get_exec_path("java") + " " + xmx_option + " -jar " + self.get_exec_path("gatk") + " -T PrintReads --allow_potentially_misencoded_quality_scores "\
                                     " -I $1 -R "+ self.reference + " -BQSR $2 -o $3 2>> $4", cmd_format='{EXE} {IN} {OUT}' )
        print_reads_map = MultiMap( print_reads, inputs=[self.input_files, self.recab_files], outputs=[self.output_files , self.print_reads_stderrs], includes=[self.reference] )
