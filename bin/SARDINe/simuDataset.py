#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT
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
__status__ = 'prod'

import os
import sys
import logging
import argparse
import subprocess
from subprocess import Popen, PIPE


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser( description='**************************************************************.' )
	parser.add_argument( '-n', '--samples-names', nargs="+", default=["withVar-splA", "withVar-splB", "withVar-splC", "withVar-splD", "withVar-splE", "withVar-splF", "withVar-splG", "withVar-splH", "withVar-splI", "withVar-splJ"], help="The *****************************************. [Default: %(default)s]" )
	parser.add_argument( '-d', '--min-distance', type=int, default=2, help="The *****************************************. [Default: %(default)s]" )
	parser.add_argument( '-v', '--version', action='version', version=__version__ )
	group_input = parser.add_argument_group( 'Inputs' ) # Inputs
	group_input.add_argument( '-m', '--input-manifest', default="INCa_V1_A_TC0070336.txt", help='*********************************** (format: *****************). [Default: %(default)s]' )
	group_input.add_argument( '-p', '--input-profile', default="variants_profiles.tsv", help='Path to the variants profile file (format: TSV). [Default: %(default)s]' )
	group_input.add_argument( '-r', '--input-runs', nargs="+", default=["/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_160928_M70265_0018_000000000-ANC7J", "/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_161026_M70265_0023_000000000-ANBWJ", "/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_161116_M70265_0027_000000000-ANC6L", "/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_161207_M70265_0031_000000000-AU8FB", "/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_170118_M70265_0038_000000000-AYAKV", "/Anapath/Illumina\ Run\ datas/Routine/INCa-V1_170113_M02726_0179_000000000-AY3E6"], help='Path to the variants profile file (format: TSV). [Default: %(default)s]' )
	group_output = parser.add_argument_group( 'Outputs' ) # Outputs
	group_output.add_argument( '-o', '--output-folder', default=os.getcwd(), help='***************************************. [Default: %(default)s]' )
	args = parser.parse_args()

	# Logger
	logging.basicConfig(
		level=logging.DEBUG,
		format='%(asctime)s -- [%(name)s][pid:%(process)d][%(levelname)s] -- %(message)s' )
	logger = logging.getLogger( os.path.basename(__file__) )

	for current_spl in args.samples_names:
		logging.info(current_spl)
		cmd = "python3 /save/fescudie/scripts/SARDINe/simuAmpliDS.py --sample-name " + current_spl + " --min-distance " + str(args.min_distance) + " --input-manifests " + args.input_manifest + " --input-genome /work/fescudie/bank/Homo_sapiens/DNA/GRCh19_Broad/ucsc.hg19.fasta --input-profile " + args.input_profile + " --output-variants " + os.path.join(args.output_folder, current_spl + ".vcf") + " --output-folder " + args.output_folder
		subprocess.check_output( cmd, shell=True )
		for sample_fasta in os.listdir( args.output_folder ):
			if sample_fasta.startswith(current_spl) and sample_fasta.endswith(".fasta"):
				logging.info("\t" + sample_fasta)
				sample_name = os.path.basename(sample_fasta).split(".")[0]
				read = "R1"
				if sample_name.endswith("R2"):
					read = "R2"
				for idx_model, model in enumerate(args.input_runs):
					model_files_regexp = model + "/Data/Intensities/BaseCalls/*_" + read + "_00*.fastq.gz"
					out_fastq_path = os.path.join( args.output_folder, "run-" + str(idx_model) + "-" + sample_name + ".fastq" )
					out_trace_path = os.path.join( args.output_folder, "run-" + str(idx_model) + "-" + sample_name + "_error.tsv" )
					cmd = "python3 /save/fescudie/scripts/SARDINe/applyQualityProfile.py --input-models " + model_files_regexp + " --input-sequences " + os.path.join(args.output_folder, sample_fasta) + " --output-sequences " + out_fastq_path + " --output-trace " + out_trace_path
					subprocess.check_output( cmd, shell=True )
					cmd = "gzip " + out_fastq_path
					subprocess.check_output( cmd, shell=True )
				os.remove( os.path.join(args.output_folder, sample_fasta) )
