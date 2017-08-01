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

import re
import json
import argparse



########################################################################
#
# FUNCTIONS
#
########################################################################
def parseFlagStat(flagstat_path):
    """
    @summary: Returns the alignment metrics from flagstat output.
    @param flagstat_path: [str] Path to the samtools flagstat output.
    @return: [dict] The alignment metrics.
    """
    # Compile regexp
    total_regex = re.compile( "^(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)\s*$" )
    secondary_regex = re.compile( "^(\d+) \+ (\d+) secondary\s*$" )
    mapped_regex = re.compile( "^(\d+) \+ (\d+) mapped \(" )
    paired_regex = re.compile( "^(\d+) \+ (\d+) paired in sequencing\s*$" )
    r1_regex = re.compile( "^(\d+) \+ (\d+) read1\s*$" )
    r2_regex = re.compile( "^(\d+) \+ (\d+) read2\s*$" )
    properlyPaired_regex = re.compile( "^(\d+) \+ (\d+) properly paired \(" )
    mateOtherChr_regex = re.compile( "^(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)\s*$" )
    
    # Parse
    aln_metrics = dict()
    nb_lines = None
    nb_secondary = None
    nb_reads = None
    with open(flagstat_path) as FH_in:
        for line in FH_in:
            line = line.strip()
            parsed = False
            # Total line
            if not parsed:
                match = total_regex.match(line)
                if match is not None:
                    parsed = True
                    nb_lines = int(match.group(1)) + int(match.group(2))
            # Secondary line
            if not parsed:
                match = secondary_regex.match(line)
                if match is not None:
                    parsed = True
                    nb_secondary = int(match.group(1)) + int(match.group(2))
                    nb_reads = nb_lines - nb_secondary
            # Mapped line
            if not parsed:
                match = mapped_regex.match(line)
                if match is not None:
                    parsed = True
                    aln_metrics['mapped'] = int(match.group(1)) + int(match.group(2)) - nb_secondary ;
            # Paired line
            if not parsed:
                match = paired_regex.match(line)
                if match is not None:
                    parsed = True
                    aln_metrics['paired'] = int(match.group(1)) + int(match.group(2))
            # R1 line
            if not parsed:
                match = r1_regex.match(line)
                if match is not None:
                    parsed = True
                    if int(match.group(1)) != 0 and nb_reads != (int(match.group(1)) + int(match.group(2)))*2: # The alignement was processed on paired-end AND single data
                        aln_metrics['nb_r1'] = nb_reads - int(match.group(1)) + int(match.group(2))
                    elif int(match.group(1)) != 0: # The alignement was processed on paired-end data
                        aln_metrics['nb_r1'] = int(match.group(1)) + int(match.group(2))
                    else: # The alignement was processed on single data
                        aln_metrics['nb_r1'] = nb_reads
            # R2 line
            if not parsed:
                match = r2_regex.match(line)
                if match is not None:
                    parsed = True
                    aln_metrics['nb_r2'] = int(match.group(1)) + int(match.group(2))
            # Properly paired line
            if not parsed:
                match = properlyPaired_regex.match(line)
                if match is not None:
                    parsed = True
                    aln_metrics['properly_paired'] = int(match.group(1)) + int(match.group(2))
            # Mate on other chr line
            if not parsed:
                match = mateOtherChr_regex.match(line)
                if match is not None:
                    parsed = True
                    aln_metrics['mate_on_other_chr'] = int(match.group(1)) + int(match.group(2))
    
    return( aln_metrics )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description="Converts the TXT output of samtools flagstat in JSON." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input', required=True, help='Path to the samtools flagstat output (format: TXT).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output', default="flagstat.json", help='Path to the output (format: JSON). [Default: %(default)s]' )
    args = parser.parse_args()

    # Process
    aln_metrics = parseFlagStat(args.input)
    with open(args.output, "w") as FH_out:
        FH_out.write( json.dumps(aln_metrics, default=lambda o: o.__dict__, sort_keys=True ) )
