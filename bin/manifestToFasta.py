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

import re
import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from sequenceIO import *



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
                    amplicons.append( fields )
    return amplicons

def getAmplicons(  manifest_path ):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons objects.
    """
    amplicons_obj = list()

    # Get amplicons by chr
    amplicons = getAmpliconsFromManifest( manifest_path )
    for ampli_idx, ampli in enumerate(amplicons):
        amplicons_obj.append(
            Amplicon(
                ampli["chromosome"],
                ampli["probe_strand"],
                ampli["ulso_sequence"],
                ampli["dlso_sequence"],
                None,
                None,
                "ampl" + str(ampli_idx) + "_" + ampli["target_region_name"].replace(" ", "_")
            )
        )

    return( amplicons_obj )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Extract primers sequence from Illumina's amplicons manifest." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--input-manifest', required=True, help='Path to the definition of the amplicons (format: Illumina manifest).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-f', '--fwd-barcodes', default="fwd_barcodes.fasta", help='******************************** (format: fasta). [Default: %(default)s]' )
    group_output.add_argument( '-r', '--rvs-barcodes', default="rvs_barcodes.fasta", help='******************************** (format: fasta). [Default: %(default)s]' )
    args = parser.parse_args()

    # Process
    amplicons = getAmplicons( args.input_manifest )
    FH_fwd = FastaIO(args.fwd_barcodes, "w")
    FH_rvs = FastaIO(args.rvs_barcodes, "w")
    for ampl in amplicons:
        record_fwd = Sequence( ampl.name, ampl.up_primer )
        FH_fwd.write( record_fwd )
        record_rvs = Sequence( ampl.name, revcom(ampl.down_primer) )
        FH_rvs.write( record_rvs )
    FH_fwd.close()
    FH_rvs.close()
