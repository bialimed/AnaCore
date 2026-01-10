#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'


import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.sequel.samplesheet import SampleSheet


class TestSampleSheet(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_SampleSheet.csv")
        self.test_cases = [
            {
                "name": "single",
                "content": """Barcodes,Bio Sample
bc1001,Alfred
bc1003,Berthold"""
            },
            {
                "name": "dual_annot",
                "content": """Barcode Name,Bio Sample Name,CT,Design
M13_bc1011--M13_bc1062,splA,100,Cov
M13_bc1011--M13_bc1063,splB,14,"""
            }
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testParse(self):
        expected = {
            "single": [
                {
                    'barcodes': {"index": "bc1001"},
                    'barcode_name': "bc1001",
                    'basename': 'Alfred',
                    'description': None,
                    'id': 'Alfred',
                    'metadata': {}
                },
                {
                    'barcodes': {"index": "bc1003"},
                    'barcode_name': "bc1003",
                    'basename': 'Berthold',
                    'description': None,
                    'id': 'Berthold',
                    'metadata': {}
                }
            ],
            "dual_annot": [
                {
                    'barcodes': {"index": "M13_bc1011", "index2": "M13_bc1062"},
                    'barcode_name': "M13_bc1011--M13_bc1062",
                    'basename': 'splA',
                    'description': None,
                    'id': 'splA',
                    'metadata': {"CT": "100", "Design": "Cov"}
                },
                {
                    'barcodes': {"index": "M13_bc1011", "index2": "M13_bc1063"},
                    'barcode_name': "M13_bc1011--M13_bc1063",
                    'basename': 'splB',
                    'description': None,
                    'id': 'splB',
                    'metadata': {"CT": "14", "Design": ""}
                }
            ]
        }
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            observed = SampleSheet(self.tmp_file)
            observed = [spl.toDict() for spl in observed.samples]
            self.assertEqual(observed, expected[curr_test["name"]])


if __name__ == "__main__":
    unittest.main()
