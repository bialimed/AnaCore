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
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class MergeSV(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_a = os.path.join( tmp_folder, unique_id + "_a.tsv")
        self.tmp_in_b = os.path.join( tmp_folder, unique_id + "_b.tsv")
        self.tmp_output = os.path.join( tmp_folder, unique_id + "_out.csv")

        # Exec command
        self.cmd = [
            "mergeSV.py",
            "--new-title", "Gene",
            "--links-titles", "Subject", "GeneName",
            "--out-separator", ";",
            "--input-files", self.tmp_in_a, self.tmp_in_b,
            "--output-file", self.tmp_output
        ]

        # Create input file a
        a_content = """Subject	Count	RPKM
ALK	180	2
ROS1	900	8
unk01	100	1"""
        with open(self.tmp_in_a, "w") as FH_out:
            FH_out.write( a_content )

        # Create input file b
        b_content = """Category	GeneName	Description
Protein Coding	ROS1	ROS Proto-Oncogene 1; Receptor Tyrosine Kinase
Protein Coding	STAT3	Signal Transducer And Activator Of Transcription 3
Protein Coding	ALK	ALK Receptor Tyrosine Kinase
Protein Coding	JAK3	Janus Kinase 3"""
        with open(self.tmp_in_b, "w") as FH_out:
            FH_out.write( b_content )


    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_a, self.tmp_in_b, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = """Gene;Count;RPKM;Category;Description
ALK;180;2;Protein Coding;ALK Receptor Tyrosine Kinase
ROS1;900;8;Protein Coding;ROS Proto-Oncogene 10x3b Receptor Tyrosine Kinase
unk01;100;1;;"""
        observed = ""
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            sorted(expected.split("\n")),
            sorted(observed.split("\n"))
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
