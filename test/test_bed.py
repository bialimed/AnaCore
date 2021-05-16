#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.bed import BEDIO, BEDRecord, getSortedAreasByChr
from anacore.region import Region


########################################################################
#
# FUNCTIONS
#
########################################################################
def diffRecord(record_a, record_b):
    are_diff = False
    attr_list = [
        "start", "end", "chrom", "score", "thickStart", "thickEnd",
        "itemRgb", "blockCount", "blockSizes", "blockStarts"
    ]
    for attr_name in attr_list:
        if str(getattr(record_a, attr_name)) != str(getattr(record_b, attr_name)):
            are_diff = True
    return are_diff


class TestBEDRecord(unittest.TestCase):
    def testSetAttrChrom(self):
        record = BEDRecord("7", 127471197, 127472363)
        record.chrom = "chr7"
        self.assertEqual(record.chrom, "chr7")
        self.assertEqual(record.reference.name, "chr7")


class TestBEDIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.bed")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.bed")

        self.expected_records = [
            BEDRecord(
                "chr7", 127471197, 127472363, "Pos1", 0, "+",
                127471197, 127472363, [255, 0, 0], None, None, None
            ),
            BEDRecord(
                "chr1", 137475865, 137477031, "Pos2", 0, "+",
                137475865, 137477031, None, None, None, None
            ),
            BEDRecord(
                "chr1", 127475865, 127477031, "Neg1", 0, "-",
                127475865, 127477031, [0, 0, 255], None, None, None
            ),
            BEDRecord(
                "chr22", 1001, 5000, "cloneA", 960, "+",
                1001, 5000, None, 2, [567, 488], [0, 3512]
            ),
            BEDRecord(
                "chr22", 2001, 6000, "cloneB", 900, "-",
                2001, 6000, None, 2, [433, 399], [0, 3601]
            )
        ]
        self.expected_content = """browser hide all
track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
chr7	127471196	127472363	Pos1	0	+	127471196	127472363	255,0,0
chr1	137475864	137477031	Pos2	0	+	137475864	137477031
chr1	127475864	127477031	Neg1	0	-	127475864	127477031	0,0,255
track name=pairedReads description="Clone Paired Reads" useScore=1
chr22	1000	5000	cloneA	960	+	1000	5000	0	2	567,488	0,3512
chr22	2000	6000	cloneB	900	-	2000	6000	0	2	433,399	0,3601"""
        with open(self.tmp_in, "w") as writer:
            writer.write(self.expected_content)

    def testGetMaxNbCol(self):
        self.assertEqual(
            BEDIO.getMaxNbCol(self.tmp_in),
            12
        )

    def testGetSortedAreasByChr(self):
        # Read
        observed = getSortedAreasByChr(self.tmp_in)
        # Expected
        expected = {
            "chr1": [
                BEDRecord(
                    "chr1", 127475865, 127477031, "Neg1", 0, "-",
                    127475865, 127477031, [0, 0, 255]
                ),
                BEDRecord(
                    "chr1", 137475865, 137477031, "Pos2", 0, "+",
                    137475865, 137477031, None
                )
            ],
            "chr7": [
                BEDRecord(
                    "chr7", 127471197, 127472363, "Pos1", 0, "+",
                    127471197, 127472363, [255, 0, 0]
                )
            ],
            "chr22": [
                BEDRecord(
                    "chr22", 1001, 5000, "cloneA", 960, "+",
                    1001, 5000, None, 2, [567, 488], [0, 3512]
                ),
                BEDRecord(
                    "chr22", 2001, 6000, "cloneB", 900, "-",
                    2001, 6000, None, 2, [433, 399], [0, 3601]
                )
            ]
        }
        # Assert
        self.assertTrue(len(expected) == len(observed))
        for chrom in expected:
            for record_a, record_b in zip(expected[chrom], observed[chrom]):
                self.assertFalse(diffRecord(record_a, record_b))

    def testIsValid(self):
        # Valid
        self.assertTrue(
            BEDIO.isValid(self.tmp_in)
        )
        # Invalid pos type
        with open(self.tmp_out, "w") as writer:
            writer.write(self.expected_content)
            writer.write("\nchr1	137475864	test	Neg2	0	-")
        self.assertFalse(
            BEDIO.isValid(self.tmp_out)
        )
        # Invalid strand value
        with open(self.tmp_out, "w") as writer:
            writer.write(self.expected_content)
            writer.write("\nchr1	137475864	137475868	Neg2	0	0")
        self.assertFalse(
            BEDIO.isValid(self.tmp_out)
        )
        # Empty RGB value
        with open(self.tmp_out, "w") as writer:
            writer.write(self.expected_content)
            writer.write("\nchr7	127471196	127472363	testRGBNone	0	+	127471196	127472363	.")
        self.assertTrue(
            BEDIO.isValid(self.tmp_out)
        )
        # Invalid RGB value
        with open(self.tmp_out, "w") as writer:
            writer.write(self.expected_content)
            writer.write("\nchr7	127471196	127472363	Pos1	0	+	127471196	127472363	255")
        self.assertFalse(
            BEDIO.isValid(self.tmp_out)
        )

    def testRead(self):
        # Read
        observed_records = []
        with BEDIO(self.tmp_in) as reader:
            for record in reader:
                observed_records.append(record)
        # Assert
        self.assertTrue(len(self.expected_records) == len(observed_records))
        for record_a, record_b in zip(self.expected_records, observed_records):
            self.assertFalse(diffRecord(record_a, record_b))

    def testWriteFromBEDRecord(self):
        expected = """chr7	127471196	127472363	Pos1	0	+	127471196	127472363	255,0,0
chr1	137475864	137477031	Pos2	0	+	137475864	137477031	.
chr1	127475864	127477031	Neg1	0	-	127475864	127477031	0,0,255"""
        with BEDIO(self.tmp_out, "w", 9) as writer:
            for curr_region in self.expected_records[:3]:
                writer.write(curr_region)
        observed = None
        with open(self.tmp_out) as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)

    def testWriteFromRegions(self):
        regions = [
            Region(158, 367, "-", "chr1", "target_01"),
            Region(54899, 54999, "+", "chr1", "target_02"),
        ]

        # From regions and over columns
        expected = """chr1\t157\t367\ttarget_01\t.\t-\t.\t.\t.
chr1\t54898\t54999\ttarget_02\t.\t+\t.\t.\t."""
        with BEDIO(self.tmp_out, "w", 9) as writer:
            for curr_region in regions:
                writer.write(curr_region)
        observed = None
        with open(self.tmp_out) as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)

        # From regions and limit columns
        expected = """chr1\t157\t367
chr1\t54898\t54999"""
        with BEDIO(self.tmp_out, "w", 3) as writer:
            for curr_region in regions:
                writer.write(curr_region)
        observed = None
        with open(self.tmp_out) as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
