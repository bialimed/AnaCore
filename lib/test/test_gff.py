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
import uuid
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.gff import GFF3IO, GFF3Record


########################################################################
#
# FUNCTIONS
#
########################################################################
def diffRecord(record_a, record_b):
    are_diff = False
    for attr_name in ["seq_id", "source", "type", "start", "end", "score", "strand", "phase"]:
        if getattr(record_a, attr_name) != getattr(record_b, attr_name):
            are_diff = True
    if sorted(record_a.annot.keys()) != sorted(record_b.annot.keys()):
        are_diff = True
    else:
        for key in sorted(record_a.annot.keys()):
            if record_a.annot[key] != record_b.annot[key]:
                are_diff = True
    return are_diff


class TestGTFIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_gff = os.path.join(tmp_folder, unique_id + "_in.gff3")
        self.tmp_out_gff = os.path.join(tmp_folder, unique_id + "_out.gff3")

        self.expected_records = [
            GFF3Record("ctg123", None, "operon", 1300, 15000, None, "+", None, {"ID": "operon001", "Name": "superOperon"}),
            GFF3Record("ctg123", None, "mRNA", 1300, 9000, None, "+", None, {"ID": "mrna0001", "Name": "sonichedgehog", "Parent": "operon001"}),
            GFF3Record("ctg123", None, "exon", 1300, 1500, None, "+", None, {"Parent": "mrna0001"}),
            GFF3Record("ctg123", None, "exon", 1050, 1500, None, "+", None, {"Parent": "mrna0001"}),
            GFF3Record("ctg123", None, "exon", 3000, 3902, None, "+", None, {"Parent": "mrna0001"}),
            GFF3Record("ctg123", None, "exon", 5000, 5500, None, "+", None, {"Parent": "mrna0001"}),
            GFF3Record("ctg123", None, "exon", 7000, 9000, None, "+", None, {"Parent": "mrna0001"}),
            GFF3Record("ctg123", None, "mRNA", 10000, 15000, None, "+", None, {"ID": "mrna0002", "Name": "subsonicsquirrel", "Parent": "operon001", "note": "test	encode and decode (cleanning characters ;, =)."}),
            GFF3Record("ctg123", None, "exon", 10000, 12000, None, "+", None, {"Parent": "mrna0002"}),
            GFF3Record("ctg123", None, "exon", 14000, 15000, None, "+", None, {"Parent": "mrna0002"})
        ]
        self.expected_content = """##gff-version 3
ctg123	.	operon	1300	15000	.	+	.	Name=superOperon;ID=operon001
ctg123	.	mRNA	1300	9000	.	+	.	Name=sonichedgehog;ID=mrna0001;Parent=operon001
ctg123	.	exon	1300	1500	.	+	.	Parent=mrna0001
ctg123	.	exon	1050	1500	.	+	.	Parent=mrna0001
ctg123	.	exon	3000	3902	.	+	.	Parent=mrna0001
ctg123	.	exon	5000	5500	.	+	.	Parent=mrna0001
ctg123	.	exon	7000	9000	.	+	.	Parent=mrna0001
ctg123	.	mRNA	10000	15000	.	+	.	Name=subsonicsquirrel;ID=mrna0002;Parent=operon001;note=test	encode and decode (cleanning characters %3B%2C %3D).
ctg123	.	exon	10000	12000	.	+	.	Parent=mrna0002
ctg123	.	exon	14000	15000	.	+	.	Parent=mrna0002"""
        with open(self.tmp_in_gff, "w") as FH_gtf:
            FH_gtf.write(self.expected_content)

    def testWrite(self):
        # Write
        with GFF3IO(self.tmp_out_gff, "w") as FH_out:
            for record in self.expected_records:
                FH_out.write(record)
        # Read
        observed_content = None
        with open(self.tmp_out_gff) as FH_out:
            observed_content = "".join(FH_out.readlines()).strip()
        # Assert
        self.assertEqual(
            self.expected_content.replace("test	encode and decode", "test encode and decode"),
            observed_content
        )

    def testRead(self):
        # Read
        observed_records = []
        with GFF3IO(self.tmp_in_gff) as FH_in:
            for record in FH_in:
                observed_records.append(record)
        # Assert
        self.assertTrue(len(self.expected_records) == len(observed_records))
        for record_a, record_b in zip(self.expected_records, observed_records):
            self.assertFalse(diffRecord(record_a, record_b))

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_gff, self.tmp_out_gff]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
