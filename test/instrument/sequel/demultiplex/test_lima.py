#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
from pathlib import Path
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.sequel.demultiplex.lima import filterSamples, samplesheetToLima


class TestFilterSamples(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_out_samplesheet = os.path.join(self.tmp_folder, unique_id + "_out_samplesheet.csv")
        self.tmp_in_samplesheet = os.path.join(self.tmp_folder, unique_id + "_in_samplesheet.csv")
        with open(self.tmp_in_samplesheet, "w") as writer:
            writer.write("""Barcode Name,CtN,Bio Sample Name
bc1011--bc1062,4.21,splA
bc1011--bc1063,4.20,splB
bc1012--bc1064,4.81,splC
bc1012--bc1065,4.54,splD""")
        self.tmp_demult_folder = os.path.join(self.tmp_folder, unique_id + "_demult")
        os.makedirs(self.tmp_demult_folder)
        Path(os.path.join(self.tmp_demult_folder, "lima_bcless.bc1011--bc1062.bam")).touch()
        Path(os.path.join(self.tmp_demult_folder, "lima_bcless.bc1012--bc1064.bam")).touch()
        Path(os.path.join(self.tmp_demult_folder, "lima_bcless.bc1012--bc1065.bam")).touch()
        Path(os.path.join(self.tmp_demult_folder, "lima_bcless.bc1015--bc1083.bam")).touch()

    def tearDown(self):
        if os.path.exists(self.tmp_in_samplesheet):
            os.remove(self.tmp_in_samplesheet)
        if os.path.exists(self.tmp_out_samplesheet):
            os.remove(self.tmp_out_samplesheet)
        if os.path.exists(self.tmp_demult_folder):
            for filename in os.listdir(self.tmp_demult_folder):
                os.remove(os.path.join(self.tmp_demult_folder, filename))
            os.rmdir(self.tmp_demult_folder)

    def test(self):
        expected = """Barcode Name,Bio Sample Name,CtN
bc1011--bc1062,splA,4.21
bc1012--bc1065,splD,4.54"""
        filterSamples(
            self.tmp_in_samplesheet, self.tmp_out_samplesheet, self.tmp_demult_folder,
            ["splA", "splB", "splD"]
        )
        with open(self.tmp_out_samplesheet) as reader:
            observed = "\n".join([line.strip() for line in reader])
        self.assertEqual(observed, expected)


class TestSamplesheetToLima(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in_samplesheet = os.path.join(self.tmp_folder, unique_id + "_in_samplesheet.csv")
        with open(self.tmp_in_samplesheet, "w") as writer:
            writer.write("""Barcode Name,CtN,Bio Sample Name
bc1011--bc1062,4.21,splA
bc1011--bc1063,4.20,splB
bc1012--bc1064,4.81,splC
bc1012--bc1065,4.54,splD""")
        self.tmp_out_samplesheet = os.path.join(self.tmp_folder, unique_id + "_out_samplesheet.csv")

    def tearDown(self):
        if os.path.exists(self.tmp_in_samplesheet):
            os.remove(self.tmp_in_samplesheet)
        if os.path.exists(self.tmp_out_samplesheet):
            os.remove(self.tmp_out_samplesheet)

    def test(self):
        expected = """Barcodes,Bio Sample
bc1011--bc1062,splA
bc1011--bc1063,splB
bc1012--bc1064,splC
bc1012--bc1065,splD"""
        samplesheetToLima(self.tmp_in_samplesheet, self.tmp_out_samplesheet)
        with open(self.tmp_out_samplesheet) as reader:
            observed = "\n".join([line.strip() for line in reader])
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
