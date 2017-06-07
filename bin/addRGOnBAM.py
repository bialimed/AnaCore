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
import sys
import pysam
import argparse



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Adds RG on all reads.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_RG = parser.add_argument_group( 'Reads groups' ) # Reads groups
    group_RG.add_argument( '--id', default="1", help='Read group identifier. The value of ID is used in the RG tags of alignment records. Must be unique among all read groups in header section. Read group IDs may be modified when merging SAM files in order to handle collisions. [Default: %(default)s]' )
    group_RG.add_argument( '--cn', help='Name of sequencing center producing the read.' )
    group_RG.add_argument( '--ds', help='Description.' )
    group_RG.add_argument( '--dt', help='Date the run was produced (ISO8601 date or date/time).' )
    group_RG.add_argument( '--fo', help='Flow order. The array of nucleotide bases that correspond to the nucleotides used for each flow of each read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other characters. Format: /\*|[ACMGRSVTWYHKDBN]+/' )
    group_RG.add_argument( '--ks', help='The array of nucleotide bases that correspond to the key sequence of each read.' )
    group_RG.add_argument( '--lb', help='Library.' )
    group_RG.add_argument( '--pg', help='Programs used for processing the read group.' )
    group_RG.add_argument( '--pi', help='Predicted median insert size.' )
    group_RG.add_argument( '--pl', choices=["CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "ONT", "PACBIO"], help='Platform/technology used to produce the reads.' )
    group_RG.add_argument( '--pm', help='Platform model. Free-form text providing further details of the platform/technology used.' )
    group_RG.add_argument( '--pu', help='Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.' )
    group_RG.add_argument( '--sm', help='Sample. Use pool name where a pool is being sequenced.' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-a', '--input-aln', required=True, help='The path to the alignments files (format: BAM).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-aln', required=True, help='The path to the alignments file (format: BAM).' )
    args = parser.parse_args()

    with pysam.AlignmentFile( args.input_aln, "rb" ) as FH_in:
        # Replace RG in header
        new_header = FH_in.header.copy()
        RG = { "ID": args.id }
        if args.cn is not None: RG["CN"] = args.cn
        if args.ds is not None: RG["DS"] = args.ds
        if args.dt is not None: RG["DT"] = args.dt
        if args.fo is not None: RG["FO"] = args.fo
        if args.ks is not None: RG["KS"] = args.ks
        if args.lb is not None: RG["LB"] = args.lb
        if args.pg is not None: RG["PG"] = args.pg
        if args.pi is not None: RG["PI"] = args.pi
        if args.pl is not None: RG["PL"] = args.pl
        if args.pm is not None: RG["PM"] = args.pm
        if args.pu is not None: RG["PU"] = args.pu
        if args.sm is not None: RG["SM"] = args.sm
        new_header["RG"] = [RG]
        # Replace RG in reads
        with pysam.AlignmentFile( args.output_aln, "wb", header=new_header ) as FH_out:
            for curr_read in FH_in.fetch(until_eof=True):
                curr_read.set_tag( "RG", args.id )
                FH_out.write( curr_read )
