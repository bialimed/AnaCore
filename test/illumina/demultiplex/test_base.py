#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
import shutil
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina.demultiplex.base import _getNearestBarcode

class TestGetNearestbarcode(unittest.TestCase):
    def test(self):
        db = ["GTGGAT+GTGGAT", "TTCA", "ATTGCT+ATT"]
        self.assertEqual(  # Same barcode
            {"bc": "ATTGCT+ATT", "dist": 0},
            _getNearestBarcode("ATTGCT+ATTNNN", db)
        )  
        self.assertEqual(  # Same barcode with query extension
            {"bc": "ATTGCT+ATT", "dist": 0},
            _getNearestBarcode("ATTGCT+ATTGCC", db)
        )
        self.assertEqual(  # Barcode with one difference
            {"bc": "ATTGCT+ATT", "dist": 1},
            _getNearestBarcode("ATCGCT+ATTNNN", db)
        ) 
        self.assertEqual(  # Barcode with only one index and no difference
            {"bc": "TTCA", "dist": 0},
            _getNearestBarcode("TTCAGG+CCGTCC", db)
        ) 
        self.assertEqual(  # No corresponding
            {"bc": "TTCA", "dist": 4},
            _getNearestBarcode("NNNNNN+NNNNNN", db)
        )
        with self.assertRaises(Exception) as context:  # Insufficient number of phases
            _getNearestBarcode("NNNNNN", db)
        with self.assertRaises(Exception) as context:  # Insufficient number of cycles
            _getNearestBarcode("NNNNNN+N", db)


if __name__ == "__main__":
    unittest.main()
