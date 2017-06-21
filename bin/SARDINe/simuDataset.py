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
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse
import subprocess


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Creates datasets (several samples and runs) of simulated reads corresponding to Illumina sequencing of amplicon double strand experiment.***************' )
    parser.add_argument( '-n', '--samples-names', nargs="+", default=["splA", "splB", "splC", "splD", "splE", "splF", "splG", "splH", "splI", "splJ"], help="The names of simulated samples. [Default: %(default)s]" )
    parser.add_argument( '-s', '--id-start', type=int, default=1, help="********************************************. [Default: %(default)s]" )
    parser.add_argument( '-d', '--min-distance', type=int, default=2, help="The minimum distance between two variants. [Default: %(default)s]" )
    parser.add_argument( '-q', '--qual-penalty', type=int, default=0, help="The penalty applied to reduce the quality of the sequences produced. With 2, the quality of each base is decrease of 2 compared to the model. [Default: %(default)s]" )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--input-manifest', required=True, help='Path to the definition of the amplicons (format: Illumina manifest).' )
    group_input.add_argument( '-p', '--input-profile', required=True, help='Path to the variants profile file (format: TSV).' )
    group_input.add_argument( '-g', '--input-genome', required=True, help='Path to the genome sequences file (format: fasta).' )
    group_input.add_argument( '-r', '--input-runs', required=True, nargs="+", help='Path to the variants profile file (format: TSV).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-folder', default=os.getcwd(), help='The output folder for the dataset. [Default: %(default)s]' )
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s -- [%(name)s][pid:%(process)d][%(levelname)s] -- %(message)s' )
    logger = logging.getLogger( os.path.basename(__file__) )
    logging.info( " ".join(sys.argv) )

    # Dataset creation
    max_lib_id = args.id_start
    lib_by_spl = dict()
    for current_spl in args.samples_names:
        # Simulate variants and perfect reads
        logging.info(current_spl)
        cmd = "simuAmpliDS.py" + \
              " --is-double-strand" + \
              " --sample-name " + current_spl + \
              " --min-distance " + str(args.min_distance) + \
              " --input-manifests " + args.input_manifest + \
              " --input-genome " + args.input_genome + \
              " --input-profile " + args.input_profile + \
              " --output-variants " + os.path.join(args.output_folder, current_spl + ".vcf") + \
              " --output-folder " + args.output_folder
        logging.info( "\t[SPL " + current_spl + "] " + cmd )
        subprocess.check_output( cmd, shell=True )
        # Apply sequencer error model on reads
        for sample_fasta in os.listdir( args.output_folder ):
            if sample_fasta.startswith(current_spl) and sample_fasta.endswith(".fasta"):
                logging.info("\t" + sample_fasta)
                filename_wot_ext = os.path.basename(sample_fasta).split(".")[0]
                sample_name = filename_wot_ext[:-5]
                lib = filename_wot_ext[-4:-3]
                read = "R2" if filename_wot_ext.endswith("R2") else "R1"
                curr_lib_id = None
                if filename_wot_ext[:-3] in lib_by_spl:
                    curr_lib_id = lib_by_spl[filename_wot_ext[:-3]]
                else:
                    curr_lib_id = str(max_lib_id)
                    lib_by_spl[filename_wot_ext[:-3]] = curr_lib_id
                    max_lib_id += 1
                for idx_model, model in enumerate(args.input_runs):
                    run_folder = os.path.join( args.output_folder, "run_" + str(idx_model) )
                    if not os.path.exists(run_folder):
                        os.makedirs(run_folder)
                    model_files_regexp = model + "/Data/Intensities/BaseCalls/*_" + read + "_00*.fastq.gz"
                    out_fastq_path = os.path.join( run_folder, sample_name + "_S" + curr_lib_id + "_L001_" + read + "_001.fastq" )
                    out_trace_path = os.path.join( run_folder, sample_name + "_S" + curr_lib_id + "_L001_" + read + "_001_error.tsv" )
                    cmd = "applyQualityProfile.py" + \
                          " --qual-penalty " + str(args.qual_penalty) + \
                          " --input-models " + model_files_regexp + \
                          " --input-sequences " + os.path.join(args.output_folder, sample_fasta) + \
                          " --output-sequences " + out_fastq_path + \
                          " --output-trace " + out_trace_path
                    logging.info( "\t[SPL " + current_spl + "] " + cmd )
                    subprocess.check_output( cmd, shell=True )
                    cmd = "gzip " + out_fastq_path
                    subprocess.check_output( cmd, shell=True )
                os.remove( os.path.join(args.output_folder, sample_fasta) )
