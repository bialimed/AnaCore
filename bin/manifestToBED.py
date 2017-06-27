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

import re
import os
import sys
import warnings
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

def getAmplicons( reference_path, manifest_path ):
    """
    @summary: Returns the list of amplicons from a Illumina's manifest.
    @param reference_path: [str] Path to the genome assembly where the amplicon have been defined.
    @param manifest_path: [str] Path to the manifest.
    @return: [list] The amplicons objects.
    """
    complete_amplicons = list()

    # Get amplicons by chr
    amplicons = getAmpliconsFromManifest( manifest_path )
    amplicons_by_chr = dict()
    for ampli_idx, ampli in enumerate(amplicons):
        chr = ampli["chromosome"]
        if chr not in amplicons_by_chr:
            amplicons_by_chr[chr] = list()
        amplicons_by_chr[chr].append(
            Amplicon(
                chr,
                ampli["probe_strand"],
                ampli["ulso_sequence"],
                ampli["dlso_sequence"],
                None,
                None,
                "ampl" + str(ampli_idx) + "_" + ampli["target_region_name"].replace(" ", "_")
            )
        )

    # Find amplicons coord
    FH_ref = SequenceFileReader.factory( reference_path )
    try:
        for record in FH_ref:
            record.id = "chr" + record.id ###############################################################################################
            if record.id in amplicons_by_chr:
                chr_str = record.string.upper()
                for ampli in amplicons_by_chr[record.id]:
                    up_primer = ampli.up_primer
                    down_primer = ampli.down_primer
                    if ampli.strand == "-":
                        up_primer = revcom(ampli.down_primer)
                        down_primer = revcom(ampli.up_primer)
                    # Find positions on chr
                    upstream_matches = list()
                    up_pattern = re.compile(up_primer)
                    for curr_match in up_pattern.finditer(chr_str):
                        upstream_matches.append({ "start":curr_match.start()+1, "end":curr_match.end() })
                    downstream_matches = list()
                    down_pattern = re.compile(down_primer)
                    for curr_match in down_pattern.finditer(chr_str):
                        downstream_matches.append({ "start":curr_match.start()+1, "end":curr_match.end() })
                    if len(upstream_matches) == 0 or len(downstream_matches) == 0:
                        raise Exception( "The primers '" + up_primer + "' and '" + down_primer + "' cannot be found in " + record.id )
                    # Check multiple target in chr
                    if len(upstream_matches) > 1:
                        warnings.warn( "The primer '" + up_primer + "' is found multiple twice in " + record.id )
                    if len(downstream_matches) > 1:
                        warnings.warn( "The primer '" + down_primer + "' is found multiple twice in " + record.id )
                    # Select smaller amplified fragment
                    prev_length = None
                    for curr_up in upstream_matches:
                        for curr_down in downstream_matches:
                            curr_length = curr_down["start"] - curr_up["end"]
                            if curr_length >= 0:
                                if prev_length is None or prev_length > curr_length:
                                    ampli.start = curr_up["start"]
                                    ampli.end = curr_down["end"]
                amplicons_on_chr = sorted( amplicons_by_chr[record.id], key=lambda ampl: (ampl.start, ampl.end) )
                complete_amplicons.extend( amplicons_on_chr )
                del(amplicons_by_chr[record.id])
    finally:
        FH_ref.close()

    return( complete_amplicons )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Converts an Illumina's amplicons manifest in BED format." )
    parser.add_argument( '-p', '--without-primers', action='store_true', help='Start and end position include only interest area (primers are excluded).' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-m', '--input-manifest', required=True, help='Path to the definition of the amplicons (format: Illumina manifest).' )
    group_input.add_argument( '-g', '--input-genome', required=True, help='Path to the genome assembly where the amplicon have been defined (format: fasta).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-BED', default="amplicons.bed", help='The amplicons description (format: BED). [Default: %(default)s]' )
    args = parser.parse_args()

    amplicons = getAmplicons( args.input_genome, args.input_manifest )
    with open(args.output_BED, "w") as FH_out:
        if args.without_primers:
            for ampl in amplicons:
                print(
                    ampl.reference,
                    (ampl.getInterestStart() - 1), # In BED start is zero-based starting position
                    ampl.getInterestEnd(), # In BED end is one-based starting position
                    ampl.name,
                    0,
                    ampl.strand,
                    sep="\t", file=FH_out )
        else:
            for ampl in amplicons:
                print(
                    ampl.reference,
                    (ampl.start - 1), # In BED start is zero-based starting position
                    ampl.end, # In BED end is one-based starting position
                    ampl.name,
                    0,
                    ampl.strand,
                    (ampl.getInterestStart() - 1), # In BED start is zero-based starting position
                    ampl.getInterestEnd(), # In BED end is one-based starting position
                    sep="\t", file=FH_out )
