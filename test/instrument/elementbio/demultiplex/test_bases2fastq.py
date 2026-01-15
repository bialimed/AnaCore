#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import datetime
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.elementbio.demultiplex.bases2fastq import DemultLog, DemultStat


class TestDemultLogLog(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_bases2fastq.log")
        self.test_cases = [
            {  # Complete with warnings
                "content": """<2025-12-02 09:32:33> [info]: bases2fastq version: 2.3.0.2116803307, use subject to license available at elementbiosciences.com
<2025-12-02 09:32:33> [info]: Run command: /soft/bases2fastq/v2.3.0/bases2fastq /raw/aviti_run4 /work/aviti_run4 --num-threads 20 --group-fastq --no-projects 
<2025-12-02 09:32:33> [info]: Parsing run manifest /raw/aviti_run4/RunManifest.csv
<2025-12-02 09:33:04> [info]: Maintaining orientation of I1 and I2 sequence(s) in lane 1
<2025-12-02 09:33:05> [info]: Maintaining orientation of I1 and I2 sequence(s) in lane 2
<2025-12-02 09:33:05> [info]: Skipping adapter sequence detection
<2025-12-02 09:33:05> [info]: Starting FASTQ generation for SIDEB_ADA_2510541191
<2025-12-02 09:33:05> [info]: Processing tile L1R01C02S1
<2025-12-02 09:33:05> [info]: Processing tile L1R01C02S2
<2025-12-02 09:33:05> [info]: Processing tile L1R01C03S1
<2025-12-02 09:33:05> [info]: Processing tile L1R01C03S2
<2025-12-02 09:33:05> [info]: Processing tile L1R01C04S1
<2025-12-02 09:33:05> [info]: Processing tile L1R01C04S2
<2025-12-02 09:33:05> [info]: Processing tile L1R02C01S1
<2025-12-02 09:33:05> [info]: Processing tile L1R02C01S2
<2025-12-02 09:33:05> [info]: Processing tile L1R02C02S1
<2025-12-02 09:33:05> [info]: Processing tile L1R02C02S2
<2025-12-02 10:39:56> [info]: Aggregating run stats
<2025-12-02 10:40:01> [info]: [SIDEB_ADA_2510541191] Generating HTML QC-report
<2025-12-02 10:40:01> [warning]: Failed to generate HTML QC report. To suppress this warning, use option --skip-qc-report. Errors from QC report generation can be found in /work/aviti_run4/info/SIDEB_ADA_2510541191_QC.errors
<2025-12-02 10:40:01> [info]: [SIDEB_ADA_2510541191] Generating MultiQC report
<2025-12-02 10:40:01> [warning]: Failed to run MultiQC. See log for details: /work/aviti_run4/info/MultiQCReport.log
<2025-12-02 10:40:02> [info]: ========= Stats Summary =========
<2025-12-02 10:40:02> [info]: Polony count:      1,172,877,956
<2025-12-02 10:40:02> [info]: Reads assigned:    96.068%
<2025-12-02 10:40:02> [info]: Mean Q score:      41.693
<2025-12-02 10:40:02> [info]: Percent Q30:       94.226%
<2025-12-02 10:40:02> [info]: Percent Q40:       82.863%
<2025-12-02 10:40:02> [info]: Assigned yield:    338.008 Gb
<2025-12-02 10:40:02> [info]: =================================
<2025-12-02 10:40:02> [info]: ============ Timing =============
<2025-12-02 10:40:02> [info]: FASTQ generation:    1h7m22.974s
<2025-12-02 10:40:02> [info]: Stats reports:       5.804s
<2025-12-02 10:40:02> [info]: Total elapsed:       1h7m28.778s
<2025-12-02 10:40:02> [info]: =================================
<2025-12-02 10:40:02> [info]: Output stored in /work/aviti_run4""",
                "expected": {
                    "command": "/soft/bases2fastq/v2.3.0/bases2fastq /raw/aviti_run4 /work/aviti_run4 --num-threads 20 --group-fastq --no-projects",
                    "version": "2.3.0.2116803307",
                    "start_time": datetime.datetime(2025, 12, 2, 9, 32, 33),
                    "status": "COMPLETED",
                    "end_time": datetime.datetime(2025, 12, 2, 10, 40, 2)
                }
            },
            {  # Incomplete
                "content": """<2025-12-02 09:32:33> [info]: bases2fastq version: 2.3.0.2116803307, use subject to license available at elementbiosciences.com
<2025-12-02 09:32:33> [info]: Run command: /soft/bases2fastq/v2.3.0/bases2fastq /raw/aviti_run4 /work/aviti_run4 --num-threads 20 --group-fastq --no-projects 
<2025-12-02 09:32:33> [info]: Parsing run manifest /raw/aviti_run4/RunManifest.csv
<2025-12-02 09:33:04> [info]: Maintaining orientation of I1 and I2 sequence(s) in lane 1
<2025-12-02 09:33:05> [info]: Maintaining orientation of I1 and I2 sequence(s) in lane 2
<2025-12-02 09:33:05> [info]: Skipping adapter sequence detection
<2025-12-02 09:33:05> [info]: Starting FASTQ generation for SIDEB_ADA_2510541191
<2025-12-02 09:33:05> [info]: Processing tile L1R01C02S1
<2025-12-02 09:33:05> [info]: Processing tile L1R01C02S2
<2025-12-02 09:33:05> [info]: Processing tile L1R01C03S1
<2025-12-02 09:33:05> [info]: Processing tile L1R01C03S2""",
                "expected": {
                    "command": "/soft/bases2fastq/v2.3.0/bases2fastq /raw/aviti_run4 /work/aviti_run4 --num-threads 20 --group-fastq --no-projects",
                    "version": "2.3.0.2116803307",
                    "start_time": datetime.datetime(2025, 12, 2, 9, 32, 33),
                    "status": "RUNNING",
                    "end_time": None
                }
            },
            {  # Error
                "content": """<2026-01-15 10:30:30> [info]: bases2fastq version: 2.3.0.2116803307, use subject to license available at elementbiosciences.com
<2026-01-15 10:30:30> [info]: Run command: /app/bases2fastq /work/aviti_run4 /work/aviti_run4/demultiplex --group-fastq --no-projects --num-threads 20 --num-unassigned 10000
<2026-01-15 10:30:30> [info]: Parsing run manifest /work/aviti_run4/RunManifest.csv
Error: error on line(s) [], Index1 length (0) for sample (TEM-HD730-S45-2025-R2) is invalid. Based on the I1Mask and I1Cycles, the software expects the Index1 length to be (10).: error on line(s) [], Index1 length (0) for sample (TEM-HD730-S45-2025-R2) is invalid. Based on the I1Mask and I1Cycles, the software expects the Index1 length to be (10).: error on line(s) [12], Index2 length (0) is invalid for sample (25T051952-RE1). Based on the I2Mask and I2Cycles, the software expects the Index2 length to be (10).: error on line(s) [12], Index2 length (0) is invalid for sample (25T051952-RE1). Based on the I2Mask and I2Cycles, the software expects the Index2 length to be (10).
Usage:
  elemctl sequence import [file name] [flags]

Flags:
  -c, --chemistryVersion string      set's the Chemistry Version
  -e, --errorFile string             overrides the default name of the warnings and errors file (default "errors.json")
  -f, --file                         if false process and validates the file but won't write the file, errors and manifest will be shown in the console (default true)
  -h, --help                         help for import
  -i, --i1Cycles int                 set's I1 Cycles (default -1)
  -j, --i2Cycles int                 set's I2 Cycles (default -1)
  -d, --kit-info                     includes the kit configuration detail in the output file
  -k, --kitConfiguration string      set's the kit configuration
  -l, --logLevel string              set's the console output level between debug/error/quiet (default "debug")
  -z, --lowDiversity                 set's the low diversity option
  -o, --output string                indicates output file for the manifest
  -p, --preparationWorkflow string   set's the preparation workflow
  -m, --r1Cycles int                 set's R1 Cycles (default -1)
  -n, --r2Cycles int                 set's R2 Cycles (default -1)
  -x, --runParamsFile string         uses the indicated file (json) to complete the RunParameters section
  -r, --runValues stringArray        overwrites any Run Value from the manifest key,value
  -s, --settings stringArray         overwrites any Setting from the manifest key,value,lane
  -v, --verbose                      prints detail process steps and errors

<2026-01-15 10:30:30> [error]: Run manifest error: Index1 length (0) for sample (TEM-HD730-S45-2025-R2) is invalid. Based on the I1Mask and I1Cycles, the software expects the Index1 length to be (10).
<2026-01-15 10:30:30> [error]: Run manifest error: Index1 length (0) for sample (TEM-HD730-S45-2025-R2) is invalid. Based on the I1Mask and I1Cycles, the software expects the Index1 length to be (10).
<2026-01-15 10:30:30> [error]: Run manifest error: Index2 length (0) is invalid for sample (25T051952-RE1). Based on the I2Mask and I2Cycles, the software expects the Index2 length to be (10).
<2026-01-15 10:30:30> [error]: Run manifest error: Index2 length (0) is invalid for sample (25T051952-RE1). Based on the I2Mask and I2Cycles, the software expects the Index2 length to be (10).
<2026-01-15 10:30:30> [error]: Failed to parse run manifest. manifest-parser exited with status -1. See details in /work/aviti_run4/demultiplex/info/RunManifestErrors.json
terminate called after throwing an instance of 'std::runtime_error'
  what():  Failed to parse run manifest. manifest-parser exited with status -1. See details in /work/aviti_run4/demultiplex/info/RunManifestErrors.json
Aborted (core dumped)""",
                "expected": {
                    "command": "/app/bases2fastq /work/aviti_run4 /work/aviti_run4/demultiplex --group-fastq --no-projects --num-threads 20 --num-unassigned 10000",
                    "version": "2.3.0.2116803307",
                    "start_time": datetime.datetime(2026, 1, 15, 10, 30, 30),
                    "status": "FAILED",
                    "end_time": datetime.datetime(2026, 1, 15, 10, 30, 30)
                }
            }
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testParse(self):
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as handle:
                handle.write(curr_test["content"] + "\n")
            res = DemultLog(self.tmp_file)
            observed = {
                "command": res.command,
                "version": res.version,
                "start_time": res.start_time,
                "status": res.status,
                "end_time": res.end_time
            }
            self.assertEqual(observed, curr_test["expected"])


class TestDemultStat(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_stat_file = os.path.join(tmp_folder, unique_id + "_IndexAssignment.csv")
        with open(self.tmp_stat_file, "w") as writer:
            writer.write("""SampleNumber,SampleName,I1,I2,NumPoloniesAssigned,PercentPoloniesAssigned,Yield(Gb),Lane
64,splh91,AACTCTTG,AAGCCATC,80175404,11.375452,24.052455,1+2
65,splh92,GTAGTCAT,AACTCTTG,79627771,11.297752,23.888171,1+2
63,splh90,AAGCCATC,GAGATTCA,45058050,6.392929,13.517332,1+2
2,splVo1,CAAGGTGA,ATGGTTAG,25986110,3.686963,7.793348,1+2
3,splVo2,GCAGGTTC,AGAAGCAA,23901913,3.391253,7.168213,1+2
4,splVo3,AAGTGTCT,GCAGGTTC,22121276,3.138612,6.634148,1+2
1,PhiX,TGTGTCGA,TGTCTGAC,2059054,0.292143,0.617575,1+2
1,PhiX,GCACATAG,GACTACTA,2030740,0.288126,0.609078,1+2
1,PhiX,CACAGATC,ACGAGAGT,1882306,0.267065,0.564561,1+2
1,PhiX,ATGTCGCT,CTAGCTCG,1865140,0.264630,0.559412,1+2
6,splh91,AACTCTTG,AAGCCATC,80175404,11.375452,24.052455,1
7,splh92,GTAGTCAT,AACTCTTG,79627771,11.297752,23.888171,1
5,splh90,AAGCCATC,GAGATTCA,45058050,6.392929,13.517332,1
1,PhiX,TGTGTCGA,TGTCTGAC,755464,0.107187,0.226638,1
1,PhiX,GCACATAG,GACTACTA,699669,0.099270,0.209900,1
1,PhiX,CACAGATC,ACGAGAGT,684293,0.097089,0.205287,1
1,PhiX,ATGTCGCT,CTAGCTCG,664079,0.094221,0.199223,1
2,splVo1,CAAGGTGA,ATGGTTAG,22634514,3.211432,6.787882,2
3,splVo2,GCAGGTTC,AGAAGCAA,21687657,3.077090,6.503944,2
4,splVo3,AAGTGTCT,GCAGGTTC,20426489,2.898152,6.125718,2
1,PhiX,GCACATAG,GACTACTA,1331071,0.188855,0.399179,2
1,PhiX,TGTGTCGA,TGTCTGAC,1303590,0.184956,0.390937,2
1,PhiX,ATGTCGCT,CTAGCTCG,1201061,0.170409,0.360189,2
1,PhiX,CACAGATC,ACGAGAGT,1198013,0.169977,0.359274,2""")

        self.tmp_undet_file = os.path.join(tmp_folder, unique_id + "_UnassignedSequences.csv")
        with open(self.tmp_undet_file, "w") as writer:
            writer.write("""I1,I2,% Polonies,Count,Lane
TAGACCAA,TCTTTCCC,0.022963,161841,1+2
CAAGGTGA,TCTTTCCC,0.019295,135993,1+2
AACTCTGA,AAGCCATC,0.007484,52751,1
TAGACCAA,TCTTTCCC,0.005980,42151,1
CAAGGTGA,TCTTTCCC,0.019295,135993,2
TAGACCAA,TCTTTCCC,0.016982,119690,2""")

    def tearDown(self):
        for path in [self.tmp_stat_file, self.tmp_undet_file]:
            if os.path.exists(path):
                os.remove(path)

    def testExpectedBarcodesCounts(self):
        expected = {
            "ATGTCGCT+CTAGCTCG": 1865140,
            "CACAGATC+ACGAGAGT": 1882306,
            "GCACATAG+GACTACTA": 2030740,
            "TGTGTCGA+TGTCTGAC": 2059054,
            "AACTCTTG+AAGCCATC": 80175404,
            "GTAGTCAT+AACTCTTG": 79627771,
            "AAGCCATC+GAGATTCA": 45058050,
            "CAAGGTGA+ATGGTTAG": 22634514,
            "GCAGGTTC+AGAAGCAA": 21687657,
            "AAGTGTCT+GCAGGTTC": 20426489
        }
        stats = DemultStat(self.tmp_stat_file)
        self.assertEqual(stats.expectedBarcodesCounts(), expected)

    def testSamplesCounts(self):
        expected = {
            "PhiX": 1865140 + 1882306 + 2030740 + 2059054,
            "splh91": 80175404,
            "splh92": 79627771,
            "splh90": 45058050,
            "splVo1": 22634514,
            "splVo2": 21687657,
            "splVo3": 20426489
        }
        stats = DemultStat(self.tmp_stat_file)
        self.assertEqual(stats.samplesCounts(), expected)

    def testUndeterminedCounts(self):
        expected = {
            "AACTCTGA+AAGCCATC": 52751,
            "TAGACCAA+TCTTTCCC": 42151 + 119690,
            "CAAGGTGA+TCTTTCCC": 135993
        }
        stats = DemultStat(self.tmp_stat_file)
        with self.assertRaises(Exception):
            self.assertEqual(stats.undeterminedCounts(), expected)
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.undeterminedCounts(), expected)

    def testUnexpectedBarcodes(self):
        with open(self.tmp_stat_file, "w") as writer:
            writer.write("""SampleNumber,SampleName,I1,I2,NumPoloniesAssigned,PercentPoloniesAssigned,Yield(Gb),Lane
1,splA,CTACTA,GTCGTC,182200,0.22775,,1+2
2,splB,ACTCTC,CTGTAC,148500,0.173750,,1+2
1,splA,CTACTA,GTCGTC,101700,0.254250,,1
2,splB,ACTCTC,CTGTAC,90000,0.225000,,1
1,splA,CTACTA,GTCGTC,80500,0.201250,,2
2,splB,ACTCTC,CTGTAC,58500,0.146250,,2""")
        dataset = [
            {"tag": "hidden", "line": "GTTATC,GTTGGG,0.0175,7000,1"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.0125,5000,1"},
            {"tag": "hidden", "line": "AAATCT,GTTGGG,0.006,2400,1"},
            {"tag": "hidden", "line": "TTTACT,GTTGGG,0.0055,2200,1"},
            {"tag": "hidden", "line": "GGGTTC,GGGAAG,0.005,2000,1"},
            {"tag": "spl", "line": "CTACTA,GTCGNN,0.004625,1850,1"},
            {"tag": "spl", "line": "CTACNN,GTCGTC,0.00425,1700,1"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.00375,1500,1"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.0035,1400,1"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.00325,1300,1"},
            {"tag": "noise", "line": "ATGCAT,AAATTT,0.003,1200,1"},
            {"tag": "noise", "line": "ATGCAT,GGCTAA,0.00275,1100,1"},
            {"tag": "noise", "line": "GGCTAA,ATGCAT,0.0025,1000,1"},
            {"tag": "hidden", "line": "GTTATC,GTTGGG,0.0175,7000,2"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.0125,5000,2"},
            {"tag": "hidden", "line": "AAATCT,GTTGGG,0.006,2400,2"},
            {"tag": "hidden", "line": "TTTACT,GTTGGG,0.0055,2200,2"},
            {"tag": "hidden", "line": "GGGTTC,GGGAAG,0.005,2000,2"},
            {"tag": "spl", "line": "CTACTA,GTCGNN,0.004625,1850,2"},
            {"tag": "spl", "line": "CTACNN,GTCGTC,0.00425,1700,2"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.00375,1500,2"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.0035,1400,2"},
            {"tag": "phiX", "line": "GGGGGG,AGTCAA,0.00325,1300,2"},
            {"tag": "noise", "line": "ATGCAT,AAATTT,0.003,1200,2"},
            {"tag": "noise", "line": "ATGCAT,GGCTAA,0.00275,1100,2"},
            {"tag": "noise", "line": "GGCTAA,ATGCAT,0.0025,1000,2"}
        ]
        # No hidden samples in undetermined
        expected = []
        with open(self.tmp_undet_file, "w") as writer:
            writer.write("I1,I2,% Polonies,Count,Lane\n")
            for rec in dataset:
                if rec["tag"] != "hidden":
                    writer.write(rec["line"] + "\n")
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # Hidden samples in undetermined
        expected = [
            {"spl": "GTTATC+GTTGGG", "ct": 14000, "parent": {"bc": "CTACTA+GTCGTC", "dist": 7}, "rate": 0.34},
            {"spl": "AAATCT+GTTGGG", "ct": 4800, "parent": {"bc": "CTACTA+GTCGTC", "dist": 8}, "rate": 0.92},
            {"spl": "TTTACT+GTTGGG", "ct": 4400, "parent": {"bc": "CTACTA+GTCGTC", "dist": 8}, "rate": 0.91},
            {"spl": "GGGTTC+GGGAAG", "ct": 4000, "parent": {"bc": "ACTCTC+CTGTAC", "dist": 8}, "rate": 0.6}
        ]
        with open(self.tmp_undet_file, "w") as writer:
            writer.write("I1,I2,% Polonies,Count,Lane\n")
            for rec in dataset:
                writer.write(rec["line"] + "\n")
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # Hidden samples in undetermined and error with 2nt of barcode for splC in samplesheet
        with open(self.tmp_stat_file, "w") as writer:
            writer.write("""SampleNumber,SampleName,I1,I2,NumPoloniesAssigned,PercentPoloniesAssigned,Yield(Gb),Lane
1,splA,CTACTA,GTCGTC,182200,0.22775,,1+2
2,splB,ACTCTC,CTGTAC,148500,0.185625,,1+2
3,splC,AAATCT,GAAGGG,100,0.000125,,1+2
1,splA,CTACTA,GTCGTC,101700,0.25425,,1
2,splB,ACTCTC,CTGTAC,90000,0.225,,1
3,splC,AAATCT,GAAGGG,100,0.00025,,1
1,splA,CTACTA,GTCGTC,80500,0.20125,,2
2,splB,ACTCTC,CTGTAC,58500,0.14625,,2
3,splC,AAATCT,GAAGGG,0,0,,2""")
        expected = [
            {"spl": "GTTATC+GTTGGG", "ct": 14000, "parent": {"bc": "CTACTA+GTCGTC", "dist": 7}, "rate": 0.34},
            {"spl": "AAATCT+GTTGGG", "ct": 4800, "parent": {"bc": "AAATCT+GAAGGG", "dist": 2}, "rate": 0.92},
            {"spl": "TTTACT+GTTGGG", "ct": 4400, "parent": {"bc": "AAATCT+GAAGGG", "dist": 6}, "rate": 0.91},
            {"spl": "GGGTTC+GGGAAG", "ct": 4000, "parent": {"bc": "ACTCTC+CTGTAC", "dist": 8}, "rate": 0.6}
        ]
        with open(self.tmp_undet_file, "w") as writer:
            writer.write("I1,I2,% Polonies,Count,Lane\n")
            for rec in dataset:
                writer.write(rec["line"] + "\n")
        stats = DemultStat(self.tmp_stat_file, self.tmp_undet_file)
        self.assertEqual(stats.unexpectedBarcodes(), expected)


if __name__ == "__main__":
    unittest.main()
