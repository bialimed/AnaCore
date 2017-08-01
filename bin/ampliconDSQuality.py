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
import time
import logging
import argparse
import subprocess
from subprocess import Popen, PIPE

from sequenceIO import *



########################################################################
#
# FUNCTIONS
#
########################################################################

def getSamplesFromSampleSheet( samplesheet_path ):
	samples = list()
	with open( samplesheet_path ) as FH_sheet:
		data_titles = None
		is_section_data = False
		for line in FH_sheet:
			if line.strip() == "[Data]":
				is_section_data = True
			elif is_section_data:
				if data_titles is None:
					data_titles = [field.strip() for field in line.split(",")]
				else:
					samples.append( { data_titles[idx]:field.strip() for idx, field in enumerate(line.split(",")) } )
	return( samples )

def getQualityDropMetrics( sequence_file, threshold=[0, 10, 20, 30], qual_offset=33 ):
	FH = FastqIO( sequence_file )
	desc_threshold = sorted(threshold, reverse=True)
	nb_always_at_threshold = {threshold: 0 for threshold in desc_threshold}
	for record in FH:
		integer_qualities = list()
		for qual in record.quality:
			integer_qualities.append( ord(qual) + qual_offset )
		min_quality = min(integer_qualities)
		for threshold in desc_threshold:
			if min_quality >= threshold:
				nb_always_at_threshold[threshold] += 1
	FH.close()
	return( nb_always_at_threshold )

########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Load Illumina run in database.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-s', '--samplesheet', required=True, help='The alignment file processed (format: BAM).' )
    group_input.add_argument( '-f', '--fastq-folder', help='Compute depth at list of positions or regions in specified file (format: BED).' )
    group_input.add_argument( '-r', '--run-folder', help='Compute depth at list of positions or regions in specified file (format: BED).' )
    args = parser.parse_args()

	# Get run informations
	run = {
		"mask": None, #[{"type": "R1", "length":151}, {"type": "I1", "length":7}, {"type": "I2", "length":7}, {"type": "R2", "length":151}]
		"date": None,
		"instrument": None
	}
	librairies = getSamplesFromSampleSheet( args.sample_sheet )
	
	
