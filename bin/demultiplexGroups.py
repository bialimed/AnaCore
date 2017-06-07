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
import re
import sys
import uuid
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from sequenceIO import FastqIO, FastaIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def get_seq_ids( fastq_path ):
    ids = list()
    #~ with FastqIO( fastq_path ) as FH:
        #~ for record in FH:
            #~ ids.append( record.id )
    FH = FastqIO( fastq_path )
    for record in FH:
        ids.append( record.id )
    FH.close()
    return ids

def addStartTag( in_path, out_path ):
    FH_in = FastaIO( in_path )
    FH_out = FastaIO( out_path, "w" )
    for record in FH_in:
        record.string = "^" + record.string
        FH_out.write( record )
    FH_in.close()
    FH_out.close()

def cutadapt( in_fastq, out_pattern, barcode_file, error_rate=0.1 ):
    cmd = [
        "cutadapt",
        "--error-rate", str(error_rate),
        "-g", "file:" + barcode_file,
        "--discard-untrimmed",
        "-o", out_pattern,
        in_fastq
    ]
    subprocess.check_call( cmd, stdout=subprocess.DEVNULL )

def getLibNameFromReadPath( fastq_path ):
    library_name = os.path.basename(fastq_path).split(".")[0]
    if re.search('_[rR][1-2]$', library_name):
        library_name = library_name[:-3]
    elif re.search('_[rR][1-2]_\d\d\d$', library_name):
        library_name = library_name[:-7]
    return library_name


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='**************************************.' )
    parser.add_argument( '-e', '--error-rate', default=0.01, type=float, help='************************************. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-f', '--fwd-barcodes', required=True, help='************ (format: fasta).' )
    group_input.add_argument( '-r', '--rvs-barcodes', required=True, help='************ (format: fasta).' )
    group_input.add_argument( '-R1', '--R1-path', required=True, help='************ (format: fastq).' )
    group_input.add_argument( '-R2', '--R2-path', required=True, help='************ (format: fastq).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-groups', required=True, help='**********************. (format: TSV)')
    args = parser.parse_args()

    tmp_prefix = getLibNameFromReadPath(args.R1_path) + "_" + str(uuid.uuid1())
    output_dir = os.path.dirname(args.output_groups)
    to_clean = list()

    # Demultiplex
    fwd_barcodes = os.path.join( output_dir, tmp_prefix + "_fwdBarcode.fasta" )
    addStartTag( args.fwd_barcodes, fwd_barcodes )
    rvs_barcodes = os.path.join( output_dir, tmp_prefix + "_rvsBarcode.fasta" )
    addStartTag( args.rvs_barcodes, rvs_barcodes )

    # Find barcodes
    cutadapt_out_R1 = os.path.join( output_dir, tmp_prefix + "_{name}_R1.fastq.gz" )
    cutadapt_out_R2 = os.path.join( output_dir, tmp_prefix + "_{name}_R2.fastq.gz" )
    cutadapt( args.R1_path, cutadapt_out_R1, fwd_barcodes, args.error_rate )
    cutadapt( args.R2_path, cutadapt_out_R2, rvs_barcodes, args.error_rate )
    
    # Retrieve reads by barcode
    with open(args.output_groups, "w") as FH_out:
        for filename in os.listdir(output_dir):
            if filename.startswith(tmp_prefix) and filename.endswith("_R1.fastq.gz"):
                # Select files corresponding to the current barcode
                filepath_R1 = os.path.join( output_dir, filename )
                filepath_R2 = filepath_R1[:-10] + "2.fastq.gz"
                barcode_name = filename[len(tmp_prefix)+1:-len("_R1.fastq.gz")]
                if os.path.exists( os.path.join(args.output_groups, filepath_R2) ):
                    # Select reads for the amplicon
                    R1 = set( get_seq_ids(filepath_R1) )
                    R2 = set( get_seq_ids(filepath_R2) )
                    retained_ids = R1.intersection(R2)
                    
                    # Write reads ID
                    for read_id in retained_ids:
                        FH_out.write( read_id + "\t" + barcode_name + "\n" )
                    
                    to_clean.append( filepath_R1 )
                    to_clean.append( filepath_R2 )
    for current_tmp_path in to_clean:
        os.remove( current_tmp_path )
