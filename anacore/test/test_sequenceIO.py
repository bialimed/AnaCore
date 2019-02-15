#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.dirname(os.path.dirname(CURRENT_DIR)))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastaIO, FastqIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFastaIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_multi_line = os.path.join(tmp_folder, unique_id + ".fasta")
        self.tmp_mono_line = os.path.join(tmp_folder, unique_id + ".fasta")

        # Create multi line
        content = """>seq1 test description
ATAGATAGCATCCCCCCNATAC
ATAGATAGCATCCCCCCNATAC
ATAGATAGCATCCCCCCNATAC
>seq2 test description 2
ATGAAAAAAAAAAAAANTG
ATGAAAAAAAAAAAAANTG"""
        with open(self.tmp_multi_line, "w") as FH_out:
            FH_out.write(content)

        # Create mono line
        content = """>seq1 test description
ATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATAC
>seq2 test description 2
ATGAAAAAAAAAAAAANTGATGAAAAAAAAAAAAANTG"""
        with open(self.tmp_mono_line, "w") as FH_out:
            FH_out.write(content)

    def testNbSeq(self):
        nb_seq = FastaIO.nbSeq(self.tmp_mono_line)
        self.assertEqual(2, nb_seq)
        nb_seq = FastaIO.nbSeq(self.tmp_multi_line)
        self.assertEqual(2, nb_seq)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_mono_line, self.tmp_multi_line]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


class TestFastqIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_seq = os.path.join(tmp_folder, unique_id + ".fastq")

        # Create file
        content = """@seq1 test description
ATAGATAGCATCCCCCCNATAC
+
@?>=<;:9876543210##-,+
>seq2 test description 2
ATGAAAAAAAAAAAAANTG
+
@?>=<;:987654321#/."""
        with open(self.tmp_seq, "w") as FH_out:
            FH_out.write(content)

    def testNbSeq(self):
        nb_seq = FastqIO.nbSeq(self.tmp_seq)
        self.assertEqual(2, nb_seq)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_seq]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
