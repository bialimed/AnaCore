#!/usr/bin/env python3
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
import time
import datetime
import argparse
import subprocess



########################################################################
#
# FUNCTIONS
#
########################################################################
def getRunAmpliDSCmd( in_spl_folder, out_run_folder, design ):
    """
    @summary: Returns the command to launch the amplicon double strand workflow.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param out_run_folder: [str] Path to the workflow output folder.
    @param design: [str] Path to the folder containing the files describing the amplicons.
    @return: [list] The command.
    """
    ressources_folder = "/save/fescudie/softwares/jflow/workflows/DSVF/data/"
    genome = {
        assembly: "GRCh37",
        sequences: "/work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa"
    }
    cmd = [
        "python3", "/save/fescudie/softwares/jflow/jflow_cli", "dsvfanapath",
        "--R1-end-adapter", os.path.join(ressources_folder, "Illumina_3prim_adapter.fasta"),
        "--R2-end-adapter", os.path.join(ressources_folder, "Illumina_5prim_adapter_rvc.fasta"),
        "--libA-folder", os.path.join(ressources_folder, design, "GRCh37", "libA"),
        "--libB-folder", os.path.join(ressources_folder, design, "GRCh37", "libB"),
        "--pos-ctrl-names", "HORIZON", "--pos-ctrl-names", "horizon", "--pos-ctrl-names", "Horizon",
        "--pos-ctrl-expected", os.path.join(ressources_folder, design, "GRCh37", "pos_ctrl_expected.vcf"),
        "--assembly-version", genome["assembly"],
        "--genome-seq", genome["sequences"],
        "--filters", os.path.join(ressources_folder, "ampliDS_filters.json"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    return cmd


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="This script uses an infinite loop to listen the Illumina's sequencer output folder and launch the appropriate workflow when a run is ended." )
    parser.add_argument( '-r', '--roll-time', type=int, default=(60*30), help="The time between each sequencer output folder evaluation (in seconds). [Default: %(default)s]"  )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-l', '--listened-folder', required=True, help="The sequencer output folder." )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-folder', required=True, help='Path to the workflows output folder.')
    args = parser.parse_args()

    # Process
    while True:
        time.sleep( 60 * 30 ) # 30 min
        for filename in os.listdir( args.listened_folder ):
            in_run_folder = os.path.join( args.listened_folder, filename )
            if os.path.isdir(in_run_folder):
                out_run_folder = os.path.join( args.output_folder, curr_filename )
                if not os.path.exists(out_run_folder) and os.path.exists(os.path.join(in_run_folder, "CompletedJobInfo.xml")): # The run must be process
                    design = "INCA_V1" ##################################################### Should be retrieve from run information
                    in_spl_folder = os.path.join( in_run_folder, "Data", "Intensities", "BaseCalls" )
                    cmd = getRunAmpliDSCmd( in_spl_folder, out_run_folder, design )
                    print(" ".join(cmd))
                    #~ os.mkdir( out_run_folder )
                    #~ with open(os.path.join(out_run_folder, "log.txt"), "w") as FH_log:
                        #~ FH_log.write( "[START] " + datetime.datetime.now().isoformat() + "\n" )
                        #~ FH_log.write( "[CMD]" + " ".join(cmd) + "\n" )
                    #~ subprocess.check_call( cmd )
                    #~ with open(os.path.join(out_run_folder, "launch.sh"), "a") as FH_log:
                        #~ FH_log.write( "[END] " + datetime.datetime.now().isoformat() + "\n" )
