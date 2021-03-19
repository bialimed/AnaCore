#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from copy import deepcopy
import gzip
import os
import sys
import tempfile
import unittest
import uuid

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
        self.tmp_in_metadata_tsv = os.path.join(tmp_folder, unique_id + "_in_withMeta.tsv")
        self.tmp_in_csv = os.path.join(tmp_folder, unique_id + "_in.csv")
        self.tmp_out_tsv = os.path.join(tmp_folder, unique_id + "_out.tsv")
        self.tmp_out_dsv = os.path.join(tmp_folder, unique_id + "_out.dsv")
        self.tmp_out_dsv_gz = os.path.join(tmp_folder, unique_id + "_out.dsv.gz")

        # Data
        self.data = {
            "metadata": ["release: 17 ; date: 2021-03-01", "author: anapath"],
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

        # Create metadata TSV
        metadata_tsv_content = """##release: 17 ; date: 2021-03-01
##author: anapath"""
        with open(self.tmp_in_metadata_tsv, "w") as FH_tsv:
            FH_tsv.write(metadata_tsv_content)

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

    def testIsValidWithMetadata(self):
        # Create one valid DSV file
        content = """@title: test
#CHROM.POS.ID.REF.ALT.QUAL.FILTER
20.14370.rs6054257.G.A.29.PASS
20.17330.   ..A.3.q10 low

20.1110696.rs6040355.A.G,T.67.PASS
.1230237..T..47.
20.1234567.microsat1.GTC.G,GTCT.50.
......"""
        with open(self.tmp_out_dsv, "w") as FH_out:
            FH_out.write(content)
        # Check isValid
        self.assertTrue(not SVIO.isValid(self.tmp_out_dsv, ".", "@"))
        # Create one valid DSV file
        content = """@title: test
#CHROM.POS.ID.REF.ALT.QUAL.FILTER
20.14370.rs6054257.G.A.29.PASS
20.17330.   ..A.3.q10 low
20.1110696.rs6040355.A.G,T.67.PASS
.1230237..T..47.
20.1234567.microsat1.GTC.G,GTCT.50.
......"""
        with open(self.tmp_out_dsv, "w") as FH_out:
            FH_out.write(content)
        # Check isValid
        self.assertTrue(SVIO.isValid(self.tmp_out_dsv, ".", "@"))
        # Check if file pointer is ok in reopen
        observed_rows = []
        with open(self.tmp_out_dsv) as FH_in:
            for row_idx, readed_row in enumerate(FH_in):
                observed_rows.append(readed_row)
        self.assertEqual(content, "".join(observed_rows))

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
            self.assertEqual(split_limit + 1, FH_out.current_line_nb)
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

    def testWriteHeader(self):
        # With titles and metadata
        with SVIO(self.tmp_out_tsv, "w") as writer:
            writer.titles = self.data["titles"]
            writer.metadata = self.data["metadata"]
            writer.writeHeader()
            writer.write(self.data["rows"][0])
        with SVIO(self.tmp_out_tsv, "r") as reader:
            self.assertEqual(reader.titles, self.data["titles"])
            self.assertEqual(reader.metadata, self.data["metadata"])
            self.assertEqual(reader.read(), [self.data["rows"][0]])
        # Without metadata
        with SVIO(self.tmp_out_tsv, "w") as writer:
            writer.titles = self.data["titles"]
            writer.writeHeader()
            writer.write(self.data["rows"][0])
        with SVIO(self.tmp_out_tsv, "r") as reader:
            self.assertEqual(reader.titles, self.data["titles"])
            self.assertEqual(reader.metadata, [])
            self.assertEqual(reader.read(), [self.data["rows"][0]])
        # Without titles
        with SVIO(self.tmp_out_tsv, "w", has_title=False) as writer:
            writer.metadata = self.data["metadata"]
            writer.writeHeader()
            writer.write(self.data["rows"][0])
        with SVIO(self.tmp_out_tsv, "r", has_title=False) as reader:
            self.assertEqual(reader.titles, None)
            self.assertEqual(reader.metadata, self.data["metadata"])
            self.assertEqual(reader.read(), [self.data["rows"][0]])
        # Without titles and metadata
        with SVIO(self.tmp_out_tsv, "w", has_title=False) as writer:
            writer.writeHeader()
            writer.write(self.data["rows"][0])
        with SVIO(self.tmp_out_tsv, "r", has_title=False) as reader:
            self.assertEqual(reader.titles, None)
            self.assertEqual(reader.metadata, [])
            self.assertEqual(reader.read(), [self.data["rows"][0]])
        # Only header with titles and metadata
        with SVIO(self.tmp_out_tsv, "w") as writer:
            writer.titles = self.data["titles"]
            writer.metadata = self.data["metadata"]
            writer.writeHeader()
        with SVIO(self.tmp_out_tsv, "r") as reader:
            self.assertEqual(reader.titles, self.data["titles"])
            self.assertEqual(reader.metadata, self.data["metadata"])
        # Only header without metadata
        with SVIO(self.tmp_out_tsv, "w") as writer:
            writer.titles = self.data["titles"]
            writer.writeHeader()
        with SVIO(self.tmp_out_tsv, "r") as reader:
            self.assertEqual(reader.titles, self.data["titles"])
            self.assertEqual(reader.metadata, [])
        # Only header without titles
        with SVIO(self.tmp_out_tsv, "w", has_title=False) as writer:
            writer.metadata = self.data["metadata"]
            writer.writeHeader()
        with SVIO(self.tmp_out_tsv, "r", has_title=False) as reader:
            self.assertEqual(reader.titles, None)
            self.assertEqual(reader.metadata, self.data["metadata"])
        # Empty without titles and metadata
        with SVIO(self.tmp_out_tsv, "w", has_title=False) as writer:
            writer.writeHeader()
        with SVIO(self.tmp_out_tsv, "r", has_title=False) as reader:
            self.assertEqual(reader.titles, None)
            self.assertEqual(reader.metadata, [])

    def testAppendEmptyDSV(self):
        # Write file
        with SVIO(self.tmp_out_dsv, "a", separator=".", title_starter="#") as FH_out:
            self.assertEqual(0, FH_out.current_line_nb)
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
        for curr_file in [self.tmp_in_tsv, self.tmp_in_metadata_tsv, self.tmp_in_csv, self.tmp_out_tsv, self.tmp_out_dsv, self.tmp_out_dsv_gz]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
