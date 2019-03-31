#!/usr/bin/env python3
#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.join(APP_DIR, "lib")
sys.path.append(BIN_DIR)
sys.path.append(LIB_DIR)

from anacore.vcf import VCFRecord
from mergeCoOccurVar import mergedRecord

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class MergeCoOccurVar(unittest.TestCase):
    def setUp(self):
        self.ref_seq = "ACGCAAATCTCGGCATGCCGATT"
        #               | | | | | |  |  |  |  |
        #               1 3 5 7 9 11 14 17 20 23
        self.variant_1 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            "artificial_1",  # id
            None,  # ref
            None,  # alt
            10,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.05]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [10], "DP": 100},
                "splB": {"AD": [40], "DP": 4900},
            }
        )
        self.variant_2 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            30,  # qual
            ["PASS"],  # filter
            {"AF": [0.06]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )
        self.expected_merge = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            20,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=A/T", "chr1:20=G/C"]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )

    def testMergedRecord_1_substit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "A"
        self.variant_1.alt = ["T"]
        # Variant 2
        self.variant_2.pos = 20
        self.variant_2.ref = "G"
        self.variant_2.alt = ["C"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCTCGGCATGCCG"
        self.expected_merge.alt = ["TAATCTCGGCATGCCC"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=A/T", "chr1:20=G/C"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_2_largeSubstit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["TGCA"]
        # Variant 2
        self.variant_2.pos = 10
        self.variant_2.ref = "TC"
        self.variant_2.alt = ["GG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCTC"
        self.expected_merge.alt = ["TGCACGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/TGCA", "chr1:10=TC/GG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_3_largeCloseSubstit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["TGCA"]
        # Variant 2
        self.variant_2.pos = 9
        self.variant_2.ref = "CT"
        self.variant_2.alt = ["GG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCT"
        self.expected_merge.alt = ["TGCAGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/TGCA", "chr1:9=CT/GG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_4_delIns(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["-"]
        # Variant 2
        self.variant_2.pos = 10
        self.variant_2.ref = "-"
        self.variant_2.alt = ["GGCATCT"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATC"
        self.expected_merge.alt = ["CGGCATCT"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/-", "chr1:10=-/GGCATCT"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_5_coDelIns(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["-"]
        # Variant 2
        self.variant_2.pos = 9
        self.variant_2.ref = "-"
        self.variant_2.alt = ["AGG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAAT"
        self.expected_merge.alt = ["AGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/-", "chr1:9=-/AGG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_6_insDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 7
        self.variant_2.ref = "ATC"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATC"
        self.expected_merge.alt = ["GTGTGAA"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:7=ATC/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_7_closeInsDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 6
        self.variant_2.ref = "AA"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAA"
        self.expected_merge.alt = ["GTGTGA"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:6=AA/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_8_coInsDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 5
        self.variant_2.ref = "AA"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AA"
        self.expected_merge.alt = ["GTGTG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:5=AA/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    # def testMergedRecord_9_delInsOtherSyntax(self):
    #     # Variant 1
    #     self.variant_1.pos = 4
    #     self.variant_1.ref = "CAAAT"
    #     self.variant_1.alt = ["C"]
    #     # Variant 2
    #     self.variant_2.pos = 9
    #     self.variant_2.ref = "C"
    #     self.variant_2.alt = ["CGCATCT"]
    #     # Expected merge
    #     self.expected_merge.pos = 5
    #     self.expected_merge.ref = "AAATC"
    #     self.expected_merge.alt = ["CGGCATCT"]
    #     self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:4=CAAAT/-", "chr1:9=C/CGGCATCT"]}
    #     # Eval
    #     observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
    #     self.assertEqual(
    #         strVariant(observed_merge),
    #         strVariant(self.expected_merge)
    #     )


def strVariant(var):
    info = []
    for info_key, info_val in sorted(var.info.items()):
        info.append(
            "{}:{}".format(info_key, info_val)
        )
    samples = []
    for spl_name, spl_val in sorted(var.samples.items()):
        for key, val in sorted(spl_val.items()):
            samples.append(
                "{}:{}:{}".format(spl_name, key, val)
            )
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        var.getName(),
        var.id,
        var.qual,
        sorted(var.filter),
        info,
        sorted(var.format),
        samples
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
