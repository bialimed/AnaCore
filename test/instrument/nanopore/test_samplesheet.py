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

from anacore.instrument.nanopore.samplesheet import SampleSheet


class TestSampleSheet(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_SampleSheet.csv")
        self.test_cases = [
            {
                "name": "no_sample",
                "content": """protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit
500ea8b3-e682-451f-a47e-e2e845f55c88,MN42938,FAX28220,no_sample,RUN4_VM,FLO-MIN106,"""
            }, {
                "name": "single_barcode_with_alias",
                "content": """position_id,sample_id,experiment_id,flow_cell_product_code,kit,barcode,alias
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode01,Test1
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode02,Test2
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode03,Test3"""
            }, {
                "name": "single_barcode_without_alias",
                "content": """position_id,sample_id,experiment_id,flow_cell_product_code,kit,barcode
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode01
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode02
X1,Test,Sample_sheet,FLO-MIN106,SQK-LSK109 EXP-NBD196,barcode03"""
            },
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testIsValid(self):
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            self.assertTrue(
                SampleSheet.isValid(self.tmp_file)
            )
        test_cases_invalid = [
            {
                "name": "empty",
                "content": """protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit"""
            },
            {
                "name": "empty_cr",
                "content": """protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit
"""
            },
            {
                "name": "missing_kit",
                "content": """protocol_run_id,position_id,flow_cell_id,sample_id,flow_cell_product_code
500ea8b3-e682-451f-a47e-e2e845f55c88,MN42938,FAX28220,no_sample,FLO-MIN106"""
            },
            {
                "name": "nb_column_change",
                "content": """protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit
500ea8b3-e682-451f-a47e-e2e845f55c88,MN42938,FAX28220,no_sample,RUN4_VM,FLO-MIN106,,"""
            },
        ]
        for curr_test in test_cases_invalid:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            self.assertFalse(
                SampleSheet.isValid(self.tmp_file)
            )

    def testParse(self):
        expected = {
            "no_sample": [
                {
                    'barcodes': dict(),
                    'basename': 'no_sample',
                    'description': None,
                    'id': 'no_sample',
                    'library_basename': 'no_sample',
                    'library_name': 'no_sample',
                    'metadata': {
                        "protocol_run_id": "500ea8b3-e682-451f-a47e-e2e845f55c88",
                        "position_id": "MN42938",
                        "flow_cell_id": "FAX28220",
                        "experiment_id": "RUN4_VM",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "",
                        "sample_id": "no_sample"
                    }
                }
            ],
            "single_barcode_with_alias": [
                {
                    'barcodes': {"index": "barcode01"},
                    'basename': 'X1_pass_Test1',
                    'description': None,
                    'id': 'Test1',
                    'library_basename': 'X1_pass_Test1',
                    'library_name': 'Test1',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
                },
                {
                    'barcodes': {"index": "barcode02"},
                    'basename': 'X1_pass_Test2',
                    'description': None,
                    'id': 'Test2',
                    'library_basename': 'X1_pass_Test2',
                    'library_name': 'Test2',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
                },
                {
                    'barcodes': {"index": "barcode03"},
                    'basename': 'X1_pass_Test3',
                    'description': None,
                    'id': 'Test3',
                    'library_basename': 'X1_pass_Test3',
                    'library_name': 'Test3',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
                }
            ],
            "single_barcode_without_alias": [
                {
                    'barcodes': {"index": "barcode01"},
                    'basename': 'X1_pass_barcode01',
                    'description': None,
                    'id': 'barcode01',
                    'library_basename': 'X1_pass_barcode01',
                    'library_name': 'barcode01',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
                },
                {
                    'barcodes': {"index": "barcode02"},
                    'basename': 'X1_pass_barcode02',
                    'description': None,
                    'id': 'barcode02',
                    'library_basename': 'X1_pass_barcode02',
                    'library_name': 'barcode02',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
                },
                {
                    'barcodes': {"index": "barcode03"},
                    'basename': 'X1_pass_barcode03',
                    'description': None,
                    'id': 'barcode03',
                    'library_basename': 'X1_pass_barcode03',
                    'library_name': 'barcode03',
                    'metadata': {
                        "position_id": "X1",
                        "experiment_id": "Sample_sheet",
                        "flow_cell_product_code": "FLO-MIN106",
                        "kit": "SQK-LSK109 EXP-NBD196",
                        "sample_id": "Test"
                    }
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
