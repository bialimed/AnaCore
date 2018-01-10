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
import json
import argparse

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from illumina import CompletedJobInfo



########################################################################
#
# FUNCTIONS
#
########################################################################
def writeLog(log_path, log):
    """
    @summary: Writes log information in txt.
    @param log_path: [str] Path to the log file.
    @param log: [CompletedJobInfo] The workflo information.
    """
    with open(log_path, "w") as FH_log:
        FH_log.write(
            "Workflow={}\n".format(log.workflow_name) + \
            "Version={}\n".format(log.version) + \
            "Parameters={}\n".format(json.dumps(log.parameters)) + \
            "Start_time={}\n".format(log.start_datetime) + \
            "End_time={}\n".format(log.end_datetime)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Converts MiSeq Reporter log file in Anapath log file.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input', required=True, help='The path to MiSeq Reporter log file (format: XML).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output', required=True, help='The path to the outputted log file (format: txt).')
    args = parser.parse_args()

    # Process
    log = CompletedJobInfo(args.input)
    writeLog(args.output, log)
