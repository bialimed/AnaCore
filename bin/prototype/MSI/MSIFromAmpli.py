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
import json
import argparse
import subprocess
from sequenceIO import FastqIO



########################################################################
#
# FUNCTIONS
#
########################################################################
def getBarcodes( in_barcodes ):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons information.
    """
    barcodes = list()
    with open(in_barcodes) as FH_in:
        for line in FH_in:
            if not line.startswith("#"):
                barcode_id, barcode_fwd, barcode_rvs = [elt.strip() for elt in line.split("\t")]
                barcodes.append({
                    "id": barcode_id,
                    "fwd": barcode_fwd,
                    "rvs": barcode_rvs
                })
    return barcodes

def getSampleFromFilename( filename ):
    """
    @summary: Returns the sample name from the R1 or R2 fastq filename.
    @param filename: [str] Path to the sample R1 or R2 file.
    @return: [str] The sample name.
    """
    miseq_spl_name = None
    match = re.search("(.+)_[rR][12](_\d{3})?.fastq.gz$", filename)
    if match is None:
        raise Exception("The sample name cannot be retrieved from the file '" + filename + "'.")
    else:
        miseq_spl_name = match.group(1)
        if re.search("_L\d{3}$", miseq_spl_name): # The filename contains the lane number
            miseq_spl_name = miseq_spl_name[:-5]
            if re.search("_[ATGCN]{6}$", miseq_spl_name): # The filename contains the index
                miseq_spl_name = miseq_spl_name[:-7]
    return miseq_spl_name


def revcom( seq ):
    """
    @summary: Returns the reverse complement the sequence.
    @param seq: [str] The sequence.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}
    return( "".join([complement_rules[base] for base in seq[::-1]]) )

def combinePairs( R1_file, R2_file, out_file, min_overlap=20, max_contradict_ratio=0.25 ):
    cmd = [ "./mergePairs.py",
        "--min-overlap", str(min_overlap),
        "--max-contradict-ratio", str(max_contradict_ratio),
        "--input-R1", R1_file,
        "--input-R2", R2_file,
        "--output-combined", out_file ]
    subprocess.check_call( cmd )

def demultiplex( R1_file, R2_file, barcodes, out_pattern, out_dir, error_rate=0.1 ):
    cmd = [ "demultiplex.py",
        "--error-rate", str(error_rate),
        "--barcodes", barcodes,
        "--R1-path", R1_file,
        "--R2-path", R2_file,
        "--output-filename-pattern", out_pattern,
        "--output-dir", out_dir ]
    subprocess.check_call( cmd )

def writeGraph( barcodes, barcodes_lengths_tsv, in_template, out_html, out_data ):
    data = list()
    for idx, curr_barcode in enumerate(barcodes):
        data.append( curr_barcode.copy() )
        with open(barcodes_lengths_tsv[idx]) as FH_in:
            lengths = list()
            for line in FH_in:
                lengths.append( [int(elt.strip()) for elt in line.split("\t")] )
            data[-1]["lengths"] = lengths
    with open(out_data, "w") as FH_out:
        FH_out.write( json.dumps(data, default=lambda o: o.__dict__, sort_keys=True ) )
    with open(in_template) as FH_in:
        with open(out_html, "w") as FH_out:
            for line in FH_in:
                if '"###PATH###"' in line:
                    line = line.replace( '###PATH###', os.path.basename(out_data) )
                FH_out.write( line )

def getNbSeq( in_fastq ):
    nb_seq = 0
    with FastqIO(in_fastq) as FH_in:
        for record in FH_in:
            nb_seq += 1
    return nb_seq


def getNbTrimmedSeq( in_fastq, reads_size ):
    nb_trimmed_seq = 0
    with FastqIO(in_fastq) as FH_in:
        for record in FH_in:
            if len(record.string) < reads_size:
                nb_trimmed_seq += 1
    return nb_trimmed_seq


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Produces microsatellites instability profile from NGS amplicon data.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-b', '--input-barcodes', required=True, help='********************* (format: TSV).' )
    group_input.add_argument( '-1', '--input-R1', required=True, help='The path to the R1 file (format: fastq).' )
    group_input.add_argument( '-2', '--input-R2', required=True, help='The path to the R2 file (format: fastq).' )
    group_input.add_argument( '-t', '--input-template', default="MSI_template.html", help='***************** (format: HTML).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-p', '--output-filename-pattern', help='****************************************************.')
    group_output.add_argument( '-d', '--output-directory', default=os.getcwd(), help='****************************************************.')
    args = parser.parse_args()

    # Process
    out_prefix_pattern = args.output_filename_pattern
    if out_prefix_pattern is None:
        out_prefix_pattern = getSampleFromFilename(os.path.basename(args.input_R1)) + "_<AMPLI>_<R>"
    demultiplex( args.input_R1, args.input_R2, args.input_barcodes, out_prefix_pattern + ".fastq.gz", args.output_directory, error_rate=0.15 )
    lenghts_tsv = list()
    barcodes = getBarcodes(args.input_barcodes)
    for curr_barcode in barcodes:
        curr_out_path_pattern = os.path.join( args.output_directory, out_prefix_pattern.replace("<AMPLI>", curr_barcode["id"]) )
        barcode_R1_file = curr_out_path_pattern.replace("<R>", "R1") + ".fastq.gz"
        barcode_R2_file = curr_out_path_pattern.replace("<R>", "R2") + ".fastq.gz"
        # Merge pairs
        combined_file = os.path.join(args.output_directory, out_prefix_pattern.replace("<AMPLI>", curr_barcode["id"]).replace("_<R>", "") + ".extendedFrags.fastq.gz")
        combinePairs( barcode_R1_file, barcode_R2_file, combined_file, 20, 0.25 )
        # Lengths
        hist_file = os.path.join(args.output_directory, out_prefix_pattern.replace("<AMPLI>", curr_barcode["id"]).replace("_<R>", "") + ".hist")
        nb_by_len = dict()
        with FastqIO(combined_file) as FH_in:
            for record in FH_in:
                seq_len = len(record.string)
                if seq_len not in nb_by_len:
                    nb_by_len[seq_len] = 1
                else:
                    nb_by_len[seq_len] += 1
        with open(hist_file, "w") as FH_out:
            for length in sorted(nb_by_len):
                FH_out.write( str(length) + "\t" + str(nb_by_len[length]) + "\n")
        lenghts_tsv.append( hist_file )
        # Stat
        curr_barcode["nb_demultiplexed_seq"] = getNbSeq( barcode_R1_file ) * 2
        curr_barcode["nb_trimmed_R1"] = None
        curr_barcode["nb_trimmed_R2"] = None
        curr_barcode["nb_combined_seq"] = getNbSeq( combined_file ) * 2
    writeGraph(
        barcodes,
        lenghts_tsv,
        args.input_template,
        os.path.join( args.output_directory, out_prefix_pattern.replace( "_<R>", "" ).replace( "_<AMPLI>", "" ) + "_MSI.html" ),
        os.path.join( args.output_directory, out_prefix_pattern.replace( "_<R>", "" ).replace( "_<AMPLI>", "" ) + "_lengths.json" )
    )

    # Clean temporary files
    for curr_file in [barcode_R1_file, barcode_R2_file, combined_file]:
        if os.path.exists(curr_file):
            os.remove(curr_file)
