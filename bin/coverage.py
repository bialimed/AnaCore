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
__email__ = 'escudie.frederiic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import time
import logging
import argparse
import subprocess
from subprocess import Popen, PIPE



########################################################################
#
# FUNCTIONS
#
########################################################################
class TmpFiles:
    """
    @summary: Manager for temporary files.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.dirs = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)


class Cmd:
    """
    @summary : Command wrapper.
    """
    def __init__(self, program, description, exec_parameters, version_parameters=None):
        """
        @param exec_parameters: [str] The parameters to execute the program. Two possibles syntaxes.
                                If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                Otherwise the parameters will be added after the program in command line.
        @param version_parameters: [str] The parameters to get the program version. Two possibles syntaxes.
                                   If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                   Otherwise the parameters will be added after the program in command line.
        """
        self.program = program
        self.description = description
        self.exec_parameters = exec_parameters
        self.version_parameters = version_parameters

    def get_cmd(self):
        """
        @summary : Returns the command line.
        @return : [str] The command line.
        """
        cmd = None
        if '##PROGRAM##' in self.exec_parameters:
            cmd = self.exec_parameters.replace('##PROGRAM##', self.program)
        else:
            cmd = self.program + ' ' + self.exec_parameters
        return cmd

    def get_version(self, location='stderr'):
        """
        @summary : Returns the program version number.
        @param location : [str] If the version command returns the version number on 'stdout' or on 'stderr'.
        @return : [str] version number if this is possible, otherwise this method return 'unknown'.
        """
        if self.version_parameters is None:
            return "unknown"
        else:
            try:
                cmd = self.program + ' ' + self.version_parameters
                if '##PROGRAM##' in self.exec_parameters:
                    cmd = self.version_parameters.replace('##PROGRAM##', self.program)
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                if location == 'stderr':
                    return stderr.strip()
                else:
                    return stdout.strip()
            except:
                raise Exception( "Version cannot be retrieve for the software '" + self.program + "'." )

    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        pass

    def submit(self, log_file=None):
        """
        @summary : Launch command, trace this action in log and parse results.
        @log_file : [str] Path to the sample process log file.
        """
        # Log
        if log_file is not None:
            logger.info( '[' + self.__class__.__name__ + '] DESCRIPTION ' + self.description + ' (' + os.path.basename(self.program) + ' version : ' + self.get_version() + ')' )
            logger.info( '[' + self.__class__.__name__ + '] COMMAND ' + self.get_cmd() )
            logger.info( '[' + self.__class__.__name__ + '] START' )
        # Process
        subprocess.check_output( self.get_cmd(), shell=True )
        # Log
        if log_file is not None:
            logger.info( '[' + self.__class__.__name__ + '] END' )
            # Post-process results
            self.parser(log_file)

"""
class SamtoolsDepthAnalysis(Analysis):
    #fait découpe et iter si thread > 1
    cmd = list()
    for idx_path, aln_map_path in aln_map_paths:
        cmd.append(
            SamtoolsDepthCmd(
                aln_map_path,
                output_path,
                *args,
                **kwargs
            )
        )
    submit(cmd)
"""

class SamtoolsViewCmd(Cmd):
    def __init__(self, aln_map_path, output_path, filter=None):
        # Options
        self.output = output_path
        filter_opt = ""
        if filter is not None:
            filter_opt = "-F " + filter + " "
        # Cmd
        Cmd.__init__( self,
                      "samtools",
                      "Filter alignment",
                      "view -bh " + filter_opt + aln_map_path + " > " + self.output,
                      "--version" )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        res = [line.strip() for line in Cmd.get_version(self, 'stdout').decode().split("\n")]
        return " ".join(res[:-1])

class SamtoolsDepthCmd(Cmd):
    def __init__(self, aln_map_path, output_path, bed=None, skip_zero=False, max_depth=100000):
        # Options
        self.output = output_path
        depth_opt = ""
        depth_opt = "-d " + str(max_depth) + " "
        bed_opt = ""
        if bed is not None:
            bed_opt = "-b " + bed + " "
        skip_zero_opt = ""
        if not skip_zero:
            skip_zero_opt = "-a "
        # Cmd
        Cmd.__init__( self,
                      "samtools",
                      "Coverage depth by position",
                      "depth " + skip_zero_opt + depth_opt + bed_opt + aln_map_path + " > " + self.output,
                      "--version" )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        res = [line.strip() for line in Cmd.get_version(self, 'stdout').decode().split("\n")]
        return " ".join(res[:-1])



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Computes the depth at each position or region.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-bam', required=True, help='The alignment file processed (format: BAM).' )
    group_input.add_argument( '-b', '--input-bed', help='Compute depth at list of positions or regions in specified file (format: BED).' )
    group_output = parser.add_argument_group( 'Outputs' ) # outputs
    group_output.add_argument( '-o', '--output-cov', required=True, help='The coverage file (format: TSV).' )
    args = parser.parse_args()
    ########################/save/fescudie/softwares/samtools-1.3.1/samtools depth -a -d 100000 -b test.bed 17T001312_S17.bam

    # Logger
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s -- [%(name)s][pid:%(process)d][%(levelname)s] -- %(message)s' )
    logger = logging.getLogger( os.path.basename(__file__) )

    # Process coverage
    out_dir = os.path.split( args.output_cov )[0]
    tmpFiles = TmpFiles(out_dir)
    try:
        tmp_bam = tmpFiles.add( "filtered.bam" )
        SamtoolsViewCmd(args.input_bam, tmp_bam, "256").submit( logger )
        SamtoolsDepthCmd(tmp_bam, args.output_cov, args.input_bed).submit( logger )
    finally:
        tmpFiles.deleteAll()
