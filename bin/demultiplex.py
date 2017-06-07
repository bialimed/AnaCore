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
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from sequenceIO import FastqIO



########################################################################
#
# FUNCTIONS
#
########################################################################
class Amplicon(object):
    def __init__(self, reference, strand, up_primer, down_primer, start, end, name):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.reference = reference
        self.up_primer = up_primer.upper()
        self.down_primer = down_primer.upper()

    def getInterestStart(self):
        interest_start = self.start + len(self.up_primer)
        if self.strand == "-":
            interest_start = self.start + len(self.down_primer)
        return( interest_start )

    def getInterestEnd(self):
        interest_end = self.end - len(self.down_primer)
        if self.strand == "-":
            interest_end = self.end - len(self.up_primer)
        return( interest_end )

def revcom( seq ):
    """
    @summary: Returns the reverse complement the sequence.
    @param seq: [str] The sequence.
    @return: [str] The reverse complement of the sequence.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N','W':'W','S':'S','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','D':'H','H':'D',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n','w':'w','s':'s','m':'k','k':'m','r':'y','y':'r','b':'v','v':'b','d':'h','h':'d'}

    return( "".join([complement_rules[base] for base in seq[::-1]]) )

def getAmpliconsFromManifest( manifest_path ):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons information.
    """
    amplicons = list()
    with open(manifest_path) as FH_manifest:
        section_probe = False
        probes_header = list()
        for line in FH_manifest:
            if line.strip() != "":
                if re.search('^\[\w+\]$', line.strip()) is not None:
                    if line.strip() != "[Probes]":
                        section_probe = False
                    else:
                        section_probe = True
                        probes_header = [field.strip().lower().replace(" ", "_") for field in FH_manifest.readline().split("\t")]
                elif section_probe:
                    fields = { probes_header[idx]:field.strip() for idx, field in enumerate(line.split("\t"))}
                    amplicons.append( 
                        Amplicon(
                            fields["chromosome"],
                            fields["probe_strand"],
                            fields["ulso_sequence"],
                            fields["dlso_sequence"],
                            None,
                            None,
                            fields["target_region_name"]
                        )
                    )
    return amplicons

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

def pickSequences( in_path, out_path, retained_ids ):
    dict_retained_ids = {curr_id: 1 for curr_id in retained_ids}
    #~ with FastqIO( out_path, "w" ) as FH_out:
        #~ with FastqIO( in_path ) as FH_in:
            #~ for record in FH_in:
                #~ if record.id in dict_retained_ids:
                    #~ FH_out.write( record )
    FH_out = FastqIO( out_path, "w" )
    FH_in = FastqIO( in_path )
    for record in FH_in:
        if record.id in dict_retained_ids:
            FH_out.write( record )
    FH_out.close()
    FH_in.close()

def cutadapt( in_fastq, out_fastq, adapter_seq, error_rate=0.1 ):
    cmd = [
        #"cutadapt",
        "/home/fescudie/.local/bin/cutadapt", ##########################"
        "--error-rate", str(error_rate),
        "-g", adapter_seq,
        "--discard-untrimmed",
        "-o", out_fastq,
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
    # design avec un des 2 primers proche d'un autre a 10% près c'est le premier dans la liste qui va prendre le pas ? Dans ce cas faire de la recup par couple de prime ?
    # au moins un primer si non ambigu ? ############################

    # Manage parameters
    parser = argparse.ArgumentParser( description='**************************************.' )
    parser.add_argument( '-e', '--error-rate', default=0.01, type=float, help='************************************. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--manifest-path', required=True, help='************ (format: Illumina\'s manifest).' )
    group_input.add_argument( '-R1', '--R1-path', required=True, help='************ (format: fastq).' )
    group_input.add_argument( '-R2', '--R2-path', required=True, help='************ (format: fastq).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-filename-pattern', help='**********************.')
    group_output.add_argument( '-d', '--output-dir', default=os.getcwd(), help='**********************. [Default: %(default)s]')
    args = parser.parse_args()

    # Set output variables
    if args.output_filename_pattern is None:
        args.output_filename_pattern = getLibNameFromReadPath(args.R1_path) + "_<AMPLI>_<R>.fastq.gz"

    # Load amplicons param
    amplicons = getAmpliconsFromManifest( args.manifest_path )

    # Demultiplex
    for amplicon_id, amplicon in enumerate(amplicons):
        ampl_out_R1 = os.path.join( args.output_dir, args.output_filename_pattern.replace("<R>", "R1").replace("<AMPLI>", str(amplicon_id)) )
        ampl_out_R2 = os.path.join( args.output_dir, args.output_filename_pattern.replace("<R>", "R2").replace("<AMPLI>", str(amplicon_id)) )

        # Find adapter
        cutadapt_out_R1 = os.path.join( args.output_dir, "tmp_" + args.output_filename_pattern.replace("<R>", "R1").replace("<AMPLI>", str(amplicon_id)) )
        cutadapt_out_R2 = os.path.join( args.output_dir, "tmp_" + args.output_filename_pattern.replace("<R>", "R2").replace("<AMPLI>", str(amplicon_id)) )
        cutadapt( args.R1_path, cutadapt_out_R1, "^" + amplicon.up_primer, args.error_rate )
        cutadapt( args.R2_path, cutadapt_out_R2, "^" + revcom(amplicon.down_primer), args.error_rate )

        # Select reads for the amplicon
        R1 = set( get_seq_ids(cutadapt_out_R1) )
        R2 = set( get_seq_ids(cutadapt_out_R2) )
        retained_ids = R1.intersection(R2)

        # Filter
        pickSequences( cutadapt_out_R1, ampl_out_R1, retained_ids )
        pickSequences( cutadapt_out_R2, ampl_out_R2, retained_ids )

        # Clean
        os.remove( cutadapt_out_R1 )
        os.remove( cutadapt_out_R2 )

        print( amplicon_id, len(retained_ids), amplicon.up_primer, amplicon.down_primer, amplicon.name, sep="\t" )
