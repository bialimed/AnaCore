#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import datetime
import csv
import os
import shutil
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina.demultiplex.bclconvert import DemultLog, DemultStat


class TestBclConvertLog(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in_folder = os.path.join(tmp_folder, unique_id + "_Logs")
        self.test_cases = [
            {  # Complete
                "info_content": """2023-01-26T20:37:12Z thread 1008163 Sample sheet being processed by common lib? No
2023-01-26T20:37:12Z thread 1008163 SampleSheet Settings:
2023-01-26T20:37:12Z thread 1008163   <none>
2023-01-26T20:37:12Z thread 1008163
2023-01-26T20:37:12Z thread 1008163 Sample Sheet Warning: 'Sample_Name' column is present in the sample sheet, but not enabled on the command line.
2023-01-26T20:37:12Z thread 1008163 Sample Sheet Warning: 'Sample_Project' column is present in the sample sheet, but not enabled on the command line.
2023-01-26T20:37:12Z thread 1008163 shared-thread-linux-native-asio output is disabled
2023-01-26T20:37:12Z thread 1008163 bcl-convert Version 00.000.000.4.0.3
2023-01-26T20:37:12Z thread 1008163 Copyright (c) 2014-2022 Illumina, Inc.
2023-01-26T20:37:12Z thread 1008163 Command Line: --bcl-input-directory /seq/230111_NDX550421_RUO_0378_AH7NLTAFX5 --output-directory out --no-lane-splitting true --bcl-num-conversion-threads 28 --bcl-num-compression-threads 2 --bcl-num-decompression-threads 2
2023-01-26T20:37:12Z thread 1008163 Conversion Begins.
2023-01-26T20:37:12Z thread 1008163 # CPU hw threads available: 64
2023-01-26T20:37:12Z thread 1008163 Parallel Tiles: 1. Threads Per Tile: 28
2023-01-26T20:37:12Z thread 1008163 SW compressors: 2
2023-01-26T20:37:12Z thread 1008163 SW decompressors: 2
2023-01-26T20:37:12Z thread 1008163 SW FASTQ compression level: 1
2023-01-26T20:43:39Z thread 1008163 Conversion Complete.""",
                "error_content": """""",
                "expected": {
                    "command": "bcl-convert --bcl-input-directory /seq/230111_NDX550421_RUO_0378_AH7NLTAFX5 --output-directory out --no-lane-splitting true --bcl-num-conversion-threads 28 --bcl-num-compression-threads 2 --bcl-num-decompression-threads 2",
                    "version": "4.0.3",
                    "start_time": datetime.datetime(2023, 1, 26, 20, 37, 12, tzinfo=datetime.timezone.utc),
                    "status": "COMPLETED",
                    "end_time": datetime.datetime(2023, 1, 26, 20, 43, 39, tzinfo=datetime.timezone.utc)
                }
            }, {  # Incomplete
                "info_content": """2023-01-26T20:37:12Z thread 1008163 Sample sheet being processed by common lib? No
2023-01-26T20:37:12Z thread 1008163 SampleSheet Settings:
2023-01-26T20:37:12Z thread 1008163   <none>
2023-01-26T20:37:12Z thread 1008163
2023-01-26T20:37:12Z thread 1008163 Sample Sheet Warning: 'Sample_Name' column is present in the sample sheet, but not enabled on the command line.
2023-01-26T20:37:12Z thread 1008163 Sample Sheet Warning: 'Sample_Project' column is present in the sample sheet, but not enabled on the command line.
2023-01-26T20:37:12Z thread 1008163 shared-thread-linux-native-asio output is disabled
2023-01-26T20:37:12Z thread 1008163 bcl-convert Version 00.000.000.4.0.3
2023-01-26T20:37:12Z thread 1008163 Copyright (c) 2014-2022 Illumina, Inc.
2023-01-26T20:37:12Z thread 1008163 Command Line: --bcl-input-directory /seq/230111_NDX550421_RUO_0378_AH7NLTAFX5 --output-directory out --no-lane-splitting true --bcl-num-conversion-threads 28 --bcl-num-compression-threads 2 --bcl-num-decompression-threads 2
2023-01-26T20:37:12Z thread 1008163 Conversion Begins.
2023-01-26T20:37:12Z thread 1008163 # CPU hw threads available: 64
2023-01-26T20:37:12Z thread 1008163 Parallel Tiles: 1. Threads Per Tile: 28
2023-01-26T20:37:12Z thread 1008163 SW compressors: 2
2023-01-26T20:37:12Z thread 1008163 SW decompressors: 2
2023-01-26T20:37:12Z thread 1008163 SW FASTQ compression level: 1""",
                "error_content": """""",
                "expected": {
                    "command": "bcl-convert --bcl-input-directory /seq/230111_NDX550421_RUO_0378_AH7NLTAFX5 --output-directory out --no-lane-splitting true --bcl-num-conversion-threads 28 --bcl-num-compression-threads 2 --bcl-num-decompression-threads 2",
                    "version": "4.0.3",
                    "start_time": datetime.datetime(2023, 1, 26, 20, 37, 12, tzinfo=datetime.timezone.utc),
                    "status": "RUNNING",
                    "end_time": None
                }
            },
            {  # Error
                "info_content": """""",
                "error_content": """2023-01-31T14:54:01Z thread 3238972 ERROR: RunInfo.xml does not exist at /seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/RunInfo.xml""",
                "expected": {
                    "command": None,
                    "version": None,
                    "start_time": datetime.datetime(2023, 1, 31, 14, 54, 1, tzinfo=datetime.timezone.utc),
                    "status": "FAILED",
                    "end_time": datetime.datetime(2023, 1, 31, 14, 54, 1, tzinfo=datetime.timezone.utc)
                }
            }
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_in_folder):
            shutil.rmtree(self.tmp_in_folder)

    def testParse(self):
        expected = [curr_test["expected"] for curr_test in self.test_cases]
        observed = []
        os.makedirs(self.tmp_in_folder)
        for curr_test in self.test_cases:
            with open(os.path.join(self.tmp_in_folder, "Info.log"), "w") as handle:
                handle.write(curr_test["info_content"])
            with open(os.path.join(self.tmp_in_folder, "Errors.log"), "w") as handle:
                handle.write(curr_test["error_content"])
            res = DemultLog(self.tmp_in_folder)
            observed.append({
                "command": res.command,
                "version": res.version,
                "start_time": res.start_time,
                "status": res.status,
                "end_time": res.end_time
            })
        self.assertEqual(expected, observed)


class TestDemultStatBclConvert(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_stat_file = os.path.join(tmp_folder, unique_id + "_Demultiplex_Stats.csv")
        with open(self.tmp_stat_file, "w") as writer:
            writer.write("""Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads
1,splA,CTGATCGT-GCGCATAT,1546203,1511172,35031,0,0.135021,0.176742,0.181226,0.000000
1,splB,ACTCTCGA-CTGTACCA,2585489,2522498,62991,0,0.225775,0.295023,0.325872,0.000000
2,splA,CGCATGAT-AAGCCTGA,2293141,2246343,46798,0,0.200246,0.262725,0.242100,0.000000
2,splB,ACGGAACA-ACGAGAAC,2318638,2270158,48480,0,0.202473,0.265510,0.250802,0.000000""")

        self.tmp_undet_file = os.path.join(tmp_folder, unique_id + "_Top_unknown_barcodes.csv")
        with open(self.tmp_undet_file, "w") as writer:
            writer.write("""Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads
1,GGGGGGGG,AGATCTCG,2610292,0.963872,0.227941
1,TGGACTCT,TCTTACGG,37763,0.013944,0.003298
2,TTGCATTC,AAGGCGTA,31132,0.011496,0.002719
2,TGGACTCT,TCTTACGG,28943,0.010687,0.002527""")

    def tearDown(self):
        for path in [self.tmp_stat_file, self.tmp_undet_file]:
            if os.path.exists(path):
                os.remove(path)

    def testSamplesCounts(self):
        expected = {
            "splA": 1546203 + 2293141,
            "splB": 2585489 + 2318638
        }
        stats = DemultStat(self.tmp_stat_file)
        self.assertEqual(stats.samplesCounts(), expected)

    def testUndeterminedCounts(self):
        expected = {
            "GGGGGGGG+AGATCTCG": 2610292,
            "TGGACTCT+TCTTACGG": 37763 + 28943,
            "TTGCATTC+AAGGCGTA": 31132
        }
        stats = DemultStat(self.tmp_stat_file)
        with self.assertRaises(Exception):
            self.assertEqual(stats.undeterminedCounts(), expected)
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.undeterminedCounts(), expected)

    def testUnexpectedBarcodes(self):
        expected = []
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # Increase count of undetermined
        data = []
        with open(self.tmp_undet_file) as handle:
            reader = csv.DictReader(handle, delimiter=',')
            for row in reader:
                row["# Reads"] = int(row["# Reads"]) * 60
                data.append(row)
        fieldnames = ["Lane", "index", "index2", "# Reads", "% of Unknown Barcodes", "% of All Reads"]
        with open(self.tmp_undet_file, "w") as handle:
            writer = csv.DictWriter(handle, fieldnames, delimiter=',')
            writer.writeheader()
            for row in data:
                writer.writerow(row)
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        # RNA under
        expected = [  # "GGGGGGGG+AGATCTCG" is skipped because it corresponds to without UDI
            {"spl": "TGGACTCT+TCTTACGG", "ct": (37763 + 28943) * 60},
        ]
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # RNA under with skipped
        expected = []
        self.assertEqual(stats.unexpectedBarcodes({"splA"}), expected)


if __name__ == "__main__":
    unittest.main()
