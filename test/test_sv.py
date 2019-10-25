#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT-O
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
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import gzip
import uuid
import tempfile
import unittest
from copy import deepcopy

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.sv import SVIO, HashedSVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestSVIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_tsv = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_in_csv = os.path.join(tmp_folder, unique_id + "_in.csv")
        self.tmp_out_dsv = os.path.join(tmp_folder, unique_id + "_out.dsv")
        self.tmp_out_dsv_gz = os.path.join(tmp_folder, unique_id + "_out.dsv.gz")

        # Data
        self.data = {
            "titles": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
            "rows": [
                ["20", "14370", "rs6054257", "G", "A", "29", "PASS"],
                ["20", "17330", "", "", "A", "3", "q10 low"],
                ["20", "1110696", "rs6040355", "A", "G,T", "67", "PASS"],
                ["", "1230237", "", "T", "", "47", ""],
                ["20", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
                ["", "", "", "", "", "", ""]
            ]
        }

        # Create TSV
        tsv_content = """#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER
20\t14370\trs6054257\tG\tA\t29\tPASS
20\t17330\t   \t\tA\t3\tq10 low
20\t1110696\trs6040355\tA\tG,T\t67\tPASS
\t1230237\t\tT\t\t47\t
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\t
\t\t\t\t\t\t"""
        with open(self.tmp_in_tsv, "w") as FH_tsv:
            FH_tsv.write(tsv_content)

        # Create CSV
        csv_content = tsv_content[1:].replace("\t", ";")
        with open(self.tmp_in_csv, "w") as FH_csv:
            FH_csv.write(csv_content)

    def testIsValidTrue(self):
        # Create one valid DSV file
        rows = deepcopy(self.data["rows"])
        rows.insert(0, self.data["titles"])
        rows.extend([
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""]
        ])
        with open(self.tmp_out_dsv, "w") as FH_out:
            for curr_row in rows:
                FH_out.write(".".join(curr_row) + "\n")
        # Check isValid
        self.assertEqual(True, SVIO.isValid(self.tmp_out_dsv, "."))
        # Check if fie pointer is ok in reopen
        with open(self.tmp_out_dsv) as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                expected_row = ".".join(rows[row_idx]) + "\n"
                self.assertEqual(expected_row, readed_row)

    def testIsValidFalse(self):
        # Create one invalid TSV file
        rows = deepcopy(self.data["rows"])
        rows.insert(0, self.data["titles"])
        rows.extend([
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""],
            ["1", "1234567", "microsat1", "GTC", "G,GTCT", "50", ""]
        ])
        with open(self.tmp_out_dsv, "w") as FH_out:
            for curr_row in rows:
                FH_out.write(".".join(curr_row) + "\n")
        # Check isValid
        self.assertEqual(False, SVIO.isValid(self.tmp_out_dsv, "."))
        # Check if file pointer is ok in reopen
        with open(self.tmp_out_dsv) as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                expected_row = ".".join(rows[row_idx]) + "\n"
                self.assertEqual(expected_row, readed_row)

    def testIterCSV(self):
        with SVIO(self.tmp_in_csv, separator=";", title_starter=None) as FH_in:
            # Header
            self.assertEqual(FH_in.titles, self.data["titles"])
            # Records
            for record_idx, readed_records in enumerate(FH_in):
                self.assertEqual(self.data["rows"][record_idx], readed_records)

    def testIterTSV(self):
        with SVIO(self.tmp_in_tsv, title_starter="#") as FH_in:
            # Header
            self.assertEqual(FH_in.titles, self.data["titles"])
            # Records
            for record_idx, readed_records in enumerate(FH_in):
                self.assertEqual(self.data["rows"][record_idx], readed_records)

    def testWriteDSV(self):
        # Write file
        with SVIO(self.tmp_out_dsv, "w", separator=".", title_starter="#") as FH_out:
            FH_out.titles = self.data["titles"]
            for row in self.data["rows"]:
                FH_out.write(row)
        # Assert result
        with open(self.tmp_out_dsv) as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)

    def testAppendDSV(self):
        split_limit = 3  # Numbers of row in first write
        # 1st write step
        with SVIO(self.tmp_out_dsv, "a", separator=".", title_starter="#") as FH_out:
            FH_out.titles = self.data["titles"]
            for row in self.data["rows"][:split_limit]:
                FH_out.write(row)
        # Assert result
        with open(self.tmp_out_dsv) as FH_in:
            nb_rows = -1
            for row_idx, readed_row in enumerate(FH_in):
                nb_rows += 1
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)
            self.assertEqual(split_limit, nb_rows)
        # 2nd write step
        with SVIO(self.tmp_out_dsv, "a", separator=".", title_starter="#") as FH_out:
            self.assertEqual(self.data["titles"], FH_out.titles)  # Assert titles retrieval
            for row in self.data["rows"][split_limit:]:
                FH_out.write(row)
        # Assert result
        with open(self.tmp_out_dsv) as FH_in:
            nb_rows = -1
            for row_idx, readed_row in enumerate(FH_in):
                nb_rows += 1
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)
            self.assertEqual(len(self.data["rows"]), nb_rows)

    def testAppendEmptyDSV(self):
        # Write file
        with SVIO(self.tmp_out_dsv, "a", separator=".", title_starter="#") as FH_out:
            FH_out.titles = self.data["titles"]
            for row in self.data["rows"]:
                FH_out.write(row)
        # Assert result
        with open(self.tmp_out_dsv) as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)

    def testAppendCompressedEmptyDSV(self):
        # Cretes empty file
        with gzip.open(self.tmp_out_dsv_gz, "wt") as FH_out:
            pass
        # Write file
        with SVIO(self.tmp_out_dsv_gz, "a", separator=".", title_starter="#") as FH_out:
            FH_out.titles = self.data["titles"]
            for row in self.data["rows"]:
                FH_out.write(row)
        # Assert result
        nb_rows = -1
        with gzip.open(self.tmp_out_dsv_gz, "rt") as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                nb_rows += 1
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)
        self.assertEqual(len(self.data["rows"]), nb_rows)

    def testAppendCompressedDSV(self):
        split_limit = 3  # Numbers of row in first write
        # 1st write step
        with SVIO(self.tmp_out_dsv_gz, "a", separator=".", title_starter="#") as FH_out:
            FH_out.titles = self.data["titles"]
            for row in self.data["rows"][:split_limit]:
                FH_out.write(row)
        # Assert result
        nb_rows = -1
        with gzip.open(self.tmp_out_dsv_gz, "rt") as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                nb_rows += 1
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)
        self.assertEqual(split_limit, nb_rows)
        # 2st write step
        with SVIO(self.tmp_out_dsv_gz, "a", separator=".", title_starter="#") as FH_out:
            self.assertEqual(self.data["titles"], FH_out.titles)  # Assert titles retrieval
            for row in self.data["rows"][split_limit:]:
                FH_out.write(row)
        # Assert result
        nb_rows = -1
        with gzip.open(self.tmp_out_dsv_gz, "rt") as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                nb_rows += 1
                expected_row = ""
                if row_idx == 0:
                    expected_row = "#" + ".".join(self.data["titles"]) + "\n"
                else:
                    expected_row = ".".join(self.data["rows"][row_idx - 1]) + "\n"
                self.assertEqual(expected_row, readed_row)
        self.assertEqual(len(self.data["rows"]), nb_rows)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_tsv, self.tmp_in_csv, self.tmp_out_dsv, self.tmp_out_dsv_gz]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
