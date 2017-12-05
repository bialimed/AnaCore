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
        self.tmp_R1 = os.path.join(tmp_folder, unique_id + "_R1.fastq.gz")
        self.tmp_R2 = os.path.join(tmp_folder, unique_id + "_R2.fastq.gz")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_combined.fastq.gz")

        # Exec command
        self.cmd = [
            "mergePairs.py",
            "--input-R1", self.tmp_R1,
            "--input-R2", self.tmp_R2,
            "--output-combined", self.tmp_output
        ]

        # Create input file R1
        a_content = """@seq_1	Empty

+

@seq_2	Too short
CACTAAATAGCAGCATTCA
+
HGFEDCBA@?>=<;:9876
@seq_3	Fragment length > reads length and errors 96:A>C,116:A>G
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCACATGCCAAAGTAGCTGATGTGTTGGACGTTCTAAATGAGGTAGATGAATATTCTG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC#BBBAA@@@@@@@@@=<<<<;<<<<;;;;;;;;;;;;;;;;::::::::::::98"""
        with open(self.tmp_in_a, "w") as FH_out:
            FH_out.write( a_content )
		# @seq_4	Fragment length < reads length
        # Create input file R2
        b_content = """@seq_1	Empty

+

@seq_2	Too short
TGAATGCTGCTATTTAGTG
+
FEDCBA@?>=<;:987654
@seq_3	Fragment length > reads length and errors 87:A>T,149:A>G,150:C>G
TTTTACATATTAAAGCCTCATGAGGATCACTGGCCAGTAAGTCTATTTTCTCTGAAGAACCAGAATATTCATCTACCTCATTTAGATCGTCCAATACATCAGCTACTTTGGCATTTGATTCAGACTCCCCATCATGTGAGTCATCAGAGG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=CCCCCCCCBBBBAA@@@@@@@@@=<<<<<<<<<;;;;;;;;;;;;;;;;::::::::::::98"""
        with open(self.tmp_in_b, "w") as FH_out:
            FH_out.write( b_content )


    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_R1, self.tmp_R2, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = """@seq_3	Fragment length > reads length
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGTTCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCTCATGAGGCTTTAATATGTAAAA
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
