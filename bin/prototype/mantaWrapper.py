#!/usr/bin/env python2.7
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
import argparse
import subprocess


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='******************************************.' )
    parser.add_argument( '-t', '--nb-threads', type=int, default=4, help='******************************************.' )
    parser.add_argument( '-m', '--max-mem', type=int, default=10, help='******************************************.' )
    parser.add_argument( '-r', '--is-RNA', action='store_true', help='**************************************************.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-n', '--input-normal-aln', help='*******************************.' )
    group_input.add_argument( '-i', '--input-tumor-aln', help='*******************************.' )
    group_input.add_argument( '-s', '--input-reference-seq', required=True, help='*******************************.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-folder', default="manta_out", help='*****************************************. [Default: %(default)s]' )
    args = parser.parse_args()

    # Process configuration
    cmd_cfg = None
    if args.is_RNA:
        cmd_cfg = [ "configManta.py",
            "--rna",
            "--bam", (args.input_tumor_aln if args.input_tumor_aln != None else args.input_normal_aln),
            "--referenceFasta", args.input_reference_seq,
            "--runDir", args.output_folder ]
    else:
        cmd_cfg = [ "configManta.py",
            ("" if args.input_normal_aln == None else "--normalBam"), ("" if args.input_normal_aln == None else args.input_normal_aln),
            "--tumorBam", args.input_tumor_aln,
            "--referenceFasta", args.input_reference_seq,
            "--runDir", args.output_folder ]
    subprocess.check_call( cmd_cfg )

    # Run manta
    cmd_run = [ os.path.join(args.output_folder, "runWorkflow.py"),
        "--mode", "local",
        "--jobs", str(args.nb_threads),
        "--memGb", str(args.max_mem) ]
    subprocess.check_call( cmd_run )
