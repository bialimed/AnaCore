#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import gzip
import os
import shutil
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.abstractFile import isEmpty, isGzip, checksum


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestChecksum(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.txt")

    def tearDown(self):
        for curr_file in [self.tmp_in]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        data = [
            {"content": "", "expected": "d41d8cd98f00b204e9800998ecf8427e"},
            {"content": "\n", "expected": "68b329da9893e34099c7d8ad5cb9c940"},
            {"content": "abc", "expected": "900150983cd24fb0d6963f7d28e17f72"}
        ]
        observed = list()
        for curr in data:
            with open(self.tmp_in, "w") as writer:
                writer.write(curr["content"])
            observed.append(checksum(self.tmp_in))
        self.assertEqual(observed, [elt["expected"] for elt in data])


class TestIsEmpty(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.txt")
        self.tmp_in_gz = os.path.join(tmp_folder, unique_id + "_in.txt.gz")

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_in_gz]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        data = [
            {"content": "", "expected": True},
            {"content": "\n", "expected": False},
            {"content": "abc", "expected": False}
        ]
        observed = list()
        expected = list()
        for curr in data:
            # raw
            with open(self.tmp_in, "w") as writer:
                writer.write(curr["content"])
            observed.append(isEmpty(self.tmp_in))
            expected.append(curr["expected"])
            # gzip
            with open(self.tmp_in, 'rb') as reader:
                with gzip.open(self.tmp_in_gz, 'wb') as writer:
                    shutil.copyfileobj(reader, writer)
            observed.append(isEmpty(self.tmp_in_gz))
            expected.append(curr["expected"])
        self.assertEqual(observed, expected)


class TestIsGzip(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.txt")
        self.tmp_in_gz = os.path.join(tmp_folder, unique_id + "_in2.txt")

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_in_gz]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        data = [
            "",
            "\n",
            "abc"
        ]
        observed = list()
        expected = list()
        for curr in data:
            # raw
            with open(self.tmp_in, "w") as writer:
                writer.write(curr)
            observed.append(isGzip(self.tmp_in))
            expected.append(False)
            # gzip
            with open(self.tmp_in, 'rb') as reader:
                with gzip.open(self.tmp_in_gz, 'wb') as writer:
                    shutil.copyfileobj(reader, writer)
            observed.append(isGzip(self.tmp_in_gz))
            expected.append(True)
        self.assertEqual(observed, expected)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
