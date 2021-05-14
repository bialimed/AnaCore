#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.db.homo_sapiens.accession import ChrAccession


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestChrAccession(unittest.TestCase):
    def testFromDict(self):
        data = [
            {"in": "chr1", "expected": "1"},
            {"in": "chrM", "expected": "MT"},
            {"in": "chrMT", "expected": "MT"},
            {"in": "CHR1", "expected": "1"},
            {"in": "1", "expected": "1"},
            {"in": "M", "expected": "MT"},
            {"in": "MT", "expected": "MT"},
            {"in": "CM000663", "expected": "1"},
            {"in": "CM000663.2", "expected": "1"},
            {"in": "J01415", "expected": "MT"},
            {"in": "NC_000001", "expected": "1"},
            {"in": "NC_000001.12", "expected": "1"},
            {"in": "NC_0018007", "expected": "MT"},
            {"in": "NC_012920", "expected": "MT"},
        ]
        for curr in data:
            self.assertEqual(ChrAccession.toHumanName(curr["in"]), curr["expected"])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
