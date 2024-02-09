#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import gzip
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.region import Region
from anacore.wig import getInfoFromFixedDef, getInfoFromVariableDef, WigIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def diffRecord(record_a, record_b):
    are_diff = True
    if record_a.reference.name == record_b.reference.name:
        if record_a.start == record_b.start:
            if record_a.end == record_b.end:
                if record_a.annot["value"] == record_b.annot["value"]:
                    are_diff = False
    return are_diff


class TestWigIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.wig")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.wig")

    def testGetInfoFromFixedDef(self):
        test_cases = [
            {
                "def": "fixedStep chrom=chr1 start=12 step=50",
                "expected": ("chr1", 12, 50, 1)
            },
            {
                "def": "fixedStep chrom=chr2 start=200 step=50 span=1",
                "expected": ("chr2", 200, 50, 1)
            },
            {
                "def": "fixedStep chrom=chr3 start=300 step=1 span=60",
                "expected": ("chr3", 300, 1, 60)
            }
        ]
        observed = [getInfoFromFixedDef(curr["def"]) for curr in test_cases]
        self.assertEqual(observed, [curr["expected"] for curr in test_cases])

    def testGetInfoFromVariableDef(self):
        test_cases = [
            {
                "def": "variableStep chrom=chr1",
                "expected": ("chr1", None, None, 1)
            },
            {
                "def": "variableStep chrom=chr2 span=1",
                "expected": ("chr2", None, None, 1)
            },
            {
                "def": "variableStep chrom=chr3 span=60",
                "expected": ("chr3", None, None, 60)
            }
        ]
        observed = [getInfoFromVariableDef(curr["def"]) for curr in test_cases]
        self.assertEqual(observed, [curr["expected"] for curr in test_cases])

    def testIsValid(self):
        # Standard
        with open(self.tmp_in, "w") as writer:
            writer.write("""variableStep chrom=chr1
300701 12.5
300702 12.5
variableStep chrom=chr2 span=5
300701 13
300710 18
fixedStep chrom=chr4 start=400601 step=100
11
22
33
fixedStep chrom=chr5 start=400601 step=100 span=20
11
22
33""")
        self.assertTrue(WigIO.isValid(self.tmp_in))
        # Short
        with open(self.tmp_in, "w") as writer:
            writer.write("""variableStep chrom=chr1
300701 12.5
300702 12.5
variableStep chrom=chr2 span=5""")
        self.assertTrue(WigIO.isValid(self.tmp_in))
        # gzip
        with gzip.open(self.tmp_in + ".gz", "wt") as writer:
            writer.write("""variableStep chrom=chr1
300701 12.5
300702 12.5
variableStep chrom=chr2 span=5""")
        self.assertTrue(WigIO.isValid(self.tmp_in + ".gz"))
        # Empty
        with open(self.tmp_in, "w") as writer:
            writer.write("")
        self.assertTrue(WigIO.isValid(self.tmp_in))
        # Invalid: missing def
        with open(self.tmp_in, "w") as writer:
            writer.write("300701 12.5")
        self.assertFalse(WigIO.isValid(self.tmp_in))
        # gzip invalid: missing def
        with open(self.tmp_in + ".gz", "wt") as writer:
            writer.write("300701 12.5")
        self.assertFalse(WigIO.isValid(self.tmp_in + ".gz"))
        # Invalid: error in record
        with open(self.tmp_in, "wt") as writer:
            writer.write("variableStep chrom=chr1\n300701 error")
        self.assertFalse(WigIO.isValid(self.tmp_in))
        # Invalid: error in def
        with open(self.tmp_in, "wt") as writer:
            writer.write("fixedStep chrom=chr1\n300701 error")
        self.assertFalse(WigIO.isValid(self.tmp_in))

    def testRead(self):
        with open(self.tmp_in, "w") as writer:
            writer.write("""browser hide all
variableStep chrom=chr1
300701 12.5
300702 12.5
track type=wiggle_0 name="variableStep" description="variableStep format" visibility=full autoScale=off
variableStep chrom=chr2 span=5
300701 13
300710 18
# comment
fixedStep chrom=chr4 start=400601 step=100
11
22
33
fixedStep chrom=chr5 start=400601 step=100 span=20
11
22
33""")
        expected_records = [
            Region(300701, 300701, None, "chr1", annot={"value": 12.5}),
            Region(300702, 300702, None, "chr1", annot={"value": 12.5}),
            Region(300701, 300705, None, "chr2", annot={"value": 13}),
            Region(300710, 300714, None, "chr2", annot={"value": 18}),
            Region(400601, 400601, None, "chr4", annot={"value": 11}),
            Region(400701, 400701, None, "chr4", annot={"value": 22}),
            Region(400801, 400801, None, "chr4", annot={"value": 33}),
            Region(400601, 400620, None, "chr5", annot={"value": 11}),
            Region(400701, 400720, None, "chr5", annot={"value": 22}),
            Region(400801, 400820, None, "chr5", annot={"value": 33})
        ]
        # Read
        with WigIO(self.tmp_in) as reader:
            observed_records = [record for record in reader]
        # Assert
        self.assertEqual(len(expected_records), len(observed_records))
        for record_a, record_b in zip(expected_records, observed_records):
            self.assertFalse(diffRecord(record_a, record_b))

    def testReadStrange(self):
        with open(self.tmp_in, "w") as writer:
            writer.write("""fixedStep chrom=chr0 start=200 step=100
variableStep chrom=chr1
300701 12.5
300702 12.5
variableStep chrom=chr2 span=5
fixedStep chrom=chr4 start=400601 step=100
11
22
33
variableStep chrom=chr6 span=5
fixedStep chrom=chr6 start=400601 step=100 span=5""")
        expected_records = [
            Region(300701, 300701, None, "chr1", annot={"value": 12.5}),
            Region(300702, 300702, None, "chr1", annot={"value": 12.5}),
            Region(400601, 400601, None, "chr4", annot={"value": 11}),
            Region(400701, 400701, None, "chr4", annot={"value": 22}),
            Region(400801, 400801, None, "chr4", annot={"value": 33})
        ]
        # Read
        with WigIO(self.tmp_in) as reader:
            observed_records = [record for record in reader]
        # Assert
        self.assertTrue(len(expected_records) == len(observed_records))
        for record_a, record_b in zip(expected_records, observed_records):
            self.assertFalse(diffRecord(record_a, record_b))

    def testWriteFromRecord(self):
        expected = """variableStep chrom=chr1
300701 12.5
300702 12.5
variableStep chrom=chr1 span=5
300801 13
300810 18
variableStep chrom=chr4 span=5
300 14
310 16
fixedStep chrom=chr4 start=400601 step=100 span=5
11
22
33
fixedStep chrom=chr5 start=400601 step=100 span=20
11
22
33
fixedStep chrom=chr5 start=405101 step=100 span=20
44
"""
        # Write
        with WigIO(self.tmp_out, "w") as writer:
            writer.write("chr1", 300701, 12.5)
            writer.write("chr1", 300702, 12.5)
            writer.span = 5  # change span => new def
            writer.write("chr1", 300801, 13)
            writer.write("chr1", 300810, 18)
            writer.write("chr4", 300, 14)  # change chrom => new def
            writer.write("chr4", 310, 16)
            writer.step = 100  # change step => new type
            writer.write("chr4", 400601, 11)
            writer.write("chr4", 400701, 22)
            writer.write("chr4", 400801, 33)
            writer.span = 20  # change span => new def
            writer.write("chr5", 400601, 11)
            writer.write("chr5", 400701, 22)
            writer.write("chr5", 400801, 33)
            writer.write("chr5", 405101, 44)  # start > step => new def
        # Observed
        with open(self.tmp_out) as reader:
            observed = reader.read()
        # Assert
        self.assertEqual(observed, expected)

    def testWriteFromRecordAppend(self):
        # Test append in same def
        with open(self.tmp_out, "w") as writer:
            writer.write("""variableStep chrom=chr1 span=5
300801 13
300810 18
fixedStep chrom=chr4 start=400601 step=100 span=5
11
22
33""")
        with WigIO(self.tmp_out, "a") as writer:
            writer.write("chr4", 400901, 44)
        expected = """variableStep chrom=chr1 span=5
300801 13
300810 18
fixedStep chrom=chr4 start=400601 step=100 span=5
11
22
33
44
"""
        with open(self.tmp_out) as reader:
            observed = reader.read()
        self.assertEqual(observed, expected)
        # Test append with change def
        with open(self.tmp_out, "w") as writer:
            writer.write("""variableStep chrom=chr1 span=5
300801 13
300810 18
fixedStep chrom=chr4 start=400601 step=100 span=5
11
22
33""")
        with WigIO(self.tmp_out, "a") as writer:
            writer.write("chr4", 401001, 44)
        expected = """variableStep chrom=chr1 span=5
300801 13
300810 18
fixedStep chrom=chr4 start=400601 step=100 span=5
11
22
33
fixedStep chrom=chr4 start=401001 step=100 span=5
44
"""
        with open(self.tmp_out) as reader:
            observed = reader.read()
        self.assertEqual(observed, expected)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in, self.tmp_in + ".gz", self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
