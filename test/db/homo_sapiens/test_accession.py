#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.db.homo_sapiens.accession import AssemblyAccession, ChrAccession


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestAssemblyAccession(unittest.TestCase):
    def testFromDict(self):
        data = [
            {"in": "GRCb36", "expected": "GRCh36"},
            {"in": "Hg18", "expected": "GRCh36"},
            {"in": "GRCb37", "expected": "GRCh37"},
            {"in": "GRCh37", "expected": "GRCh37"},
            {"in": "Hg19", "expected": "GRCh37"},
            {"in": "GRCh38", "expected": "GRCh38"},
            {"in": "hg38", "expected": "GRCh38"},
            {"in": "GCF_000001405.12", "expected": "GRCh36"},
            {"in": "GCF_000001405.13", "expected": "GRCh37"},
            {"in": "GCF_000001405.25", "expected": "GRCh37"},
            {"in": "GCF_000001405.26", "expected": "GRCh38"},
            {"in": "GCF_000001405.39", "expected": "GRCh38"},
            {"in": "GCA_000001405.1", "expected": "GRCh37"},
            {"in": "GCA_000001405.14", "expected": "GRCh37"},
            {"in": "GCA_000001405.15", "expected": "GRCh38"},
            {"in": "GCA_000001405.28", "expected": "GRCh38"}
        ]
        for curr in data:
            self.assertEqual(AssemblyAccession.toHumanName(curr["in"]), curr["expected"])


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
