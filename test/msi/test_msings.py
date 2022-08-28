#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU-Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.msi.base import Status
from anacore.msi.locus import Locus, LocusDataDistrib
from anacore.msi.msings import MSINGSAnalysisIO, MSINGSEval, MSINGSReport
from anacore.msi.sample import MSISample


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestMSINGSAnalysisIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")
        with open(self.tmp_in, "w") as writer:
            writer.write("""Position	Name	Average_Depth	Number_of_Peaks	Standard_Deviation	IndelLength:AlleleFraction:SupportingCalls
2:95183613-95183636	NR24	767	3	1.025526	-7:0.0016051:1 -6:0.0096308:6 -5:0.0080257:5 -4:0.0064205:4 -3:0.0304976:19 -2:0.0786517:49 -1:0.0770465:48 0:1.0000000:623 1:0.0112360:7 2:0.0064205:4 3:0.0016051:1
4:54732045-54732070	BAT25	2191	2	0.915497	-8:0.0005402:1 -7:0:0 -6:0.0027012:5 -5:0.0070232:13 -4:0.0156672:29 -3:0.0183684:34 -2:0.0459211:85 -1:0.0605078:112 0:1.0000000:1851 1:0.0199892:37 2:0.0102647:19 3:0.0021610:4 4:0.0005402:1
""")
        self.records = [
            Locus.fromDict({
                "position": "2:95183613-95183636",
                "name": "NR24",
                "results": {
                    "mSINGS": {
                        "status": Status.undetermined,
                        "data": {
                            "avg_depth": 767,
                            "nb_peaks": 3,
                            "std_dev": 1.025526,
                            "lengths": {
                                "ct_by_len": {
                                    23 - 7: 1,
                                    23 - 6: 6,
                                    23 - 5: 5,
                                    23 - 4: 4,
                                    23 - 3: 19,
                                    23 - 2: 49,
                                    23 - 1: 48,
                                    23 - 0: 623,
                                    23 + 1: 7,
                                    23 + 2: 4,
                                    23 + 3: 1
                                }
                            }
                        }
                    }
                }
            }),
            Locus.fromDict({
                "position": "4:54732045-54732070",
                "name": "BAT25",
                "results": {
                    "mSINGS": {
                        "status": Status.undetermined,
                        "data": {
                            "avg_depth": 2191,
                            "nb_peaks": 2,
                            "std_dev": 0.915497,
                            "lengths": {
                                "ct_by_len": {
                                    25 - 8: 1,
                                    25 - 6: 5,
                                    25 - 5: 13,
                                    25 - 4: 29,
                                    25 - 3: 34,
                                    25 - 2: 85,
                                    25 - 1: 112,
                                    25 - 0: 1851,
                                    25 + 1: 37,
                                    25 + 2: 19,
                                    25 + 3: 4,
                                    25 + 4: 1
                                }
                            }
                        }
                    }
                }
            })
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testRead(self):
        with MSINGSAnalysisIO(self.tmp_in) as reader:
            self.assertEqual(
                [toDict(rec) for rec in reader],
                [toDict(rec) for rec in self.records]
            )

    def testWrite(self):
        with MSINGSAnalysisIO(self.tmp_out, "w") as writer:
            for rec in self.records:
                writer.write(rec)
        with open(self.tmp_in) as reader:
            expected = reader.readlines()
        expected[2] = expected[2].replace("-7:0:0", "-7:0.0000000:0")
        with open(self.tmp_out) as reader:
            observed = reader.readlines()
        self.assertEqual(observed, expected)


class TestMSINGSEval(unittest.TestCase):
    def testGetNbPeaks(self):
        records = [
            {
                "data": LocusDataDistrib({
                    23 - 7: 1,
                    23 - 6: 6,
                    23 - 5: 5,
                    23 - 4: 4,
                    23 - 3: 19,
                    23 - 2: 49,
                    23 - 1: 48,
                    23 - 0: 623,
                    23 + 1: 7,
                    23 + 2: 4,
                    23 + 3: 1
                }),
                "expected": 3
            },
            {
                "data": LocusDataDistrib({
                    25 - 8: 1,
                    25 - 6: 5,
                    25 - 5: 13,
                    25 - 4: 29,
                    25 - 3: 34,
                    25 - 2: 85,
                    25 - 1: 112,
                    25 - 0: 1851,
                    25 + 1: 37,
                    25 + 2: 19,
                    25 + 3: 4,
                    25 + 4: 1
                }),
                "expected": 2
            }
        ]
        expected = list()
        observed = list()
        for curr in records:
            expected.append(curr["expected"])
            observed.append(
                MSINGSEval.getNbPeaks(curr["data"])
            )
        self.assertEqual(observed, expected)

    def testGetThresholdFromNbPeaks(self):
        records = [
            {"data": [3, 3, 4], "expected": 0.47140452079103168 * 2.0 +	3.3333333333333335},
            {"data": [2, 2, 2], "expected": 0.0 * 2.0 + 2.0},
            {"data": [3, 2, 4], "expected": 0.81649658092772603 * 2.0 + 3.0},
            {"data": [5, 4, 5], "expected": 0.47140452079103168 * 2.0 + 4.666666666666667},
            {"data": [1, 3, 2], "expected": 0.81649658092772603 * 2.0 + 2.0}
        ]
        expected = list()
        observed = list()
        for curr in records:
            expected.append(round(curr["expected"], 7))
            observed.append(
                round(
                    MSINGSEval.getThresholdFromNbPeaks(curr["data"]),
                    7
                )
            )
        self.assertEqual(observed, expected)


class TestMSINGSReport(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")
        with open(self.tmp_in, "w") as writer:
            writer.write("""Position	splA	splB
unstable_loci	3	6
covered_loci	7	7
msings_score	0.4286	0.8571
msi_status	NEG	POS
2:47414420-47414447	Unstable	Unstable
2:95183613-95183636	Stable	Unstable
4:54732045-54732070	Unstable	Unstable
11:102322777-102322803	Unstable	Unstable
11:125620870-125620891	Stable	Unstable
13:31148483-31148500	Stable	Unstable
14:23183137-23183158	Stable	Stable
""")
        self.records = [
            MSISample.fromDict({
                "name": "splA",
                "results": {
                    "mSINGS": {
                        "method": "mSINGS",
                        "status": Status.stable,
                        "score": 1 - 0.4286
                    }
                },
                "loci": {
                    "2:47414420-47414447": {
                        "position": "2:47414420-47414447",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "2:95183613-95183636": {
                        "position": "2:95183613-95183636",
                        "results": {"mSINGS": {"status": Status.stable}}
                    },
                    "4:54732045-54732070": {
                        "position": "4:54732045-54732070",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "11:102322777-102322803": {
                        "position": "11:102322777-102322803",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "11:125620870-125620891": {
                        "position": "11:125620870-125620891",
                        "results": {"mSINGS": {"status": Status.stable}}
                    },
                    "13:31148483-31148500": {
                        "position": "13:31148483-31148500",
                        "results": {"mSINGS": {"status": Status.stable}}
                    },
                    "14:23183137-23183158": {
                        "position": "14:23183137-23183158",
                        "results": {"mSINGS": {"status": Status.stable}}
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splB",
                "results": {
                    "mSINGS": {
                        "method": "mSINGS",
                        "status": Status.unstable,
                        "score": 0.8571
                    }
                },
                "loci": {
                    "2:47414420-47414447": {
                        "position": "2:47414420-47414447",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "2:95183613-95183636": {
                        "position": "2:95183613-95183636",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "4:54732045-54732070": {
                        "position": "4:54732045-54732070",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "11:102322777-102322803": {
                        "position": "11:102322777-102322803",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "11:125620870-125620891": {
                        "position": "11:125620870-125620891",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "13:31148483-31148500": {
                        "position": "13:31148483-31148500",
                        "results": {"mSINGS": {"status": Status.unstable}}
                    },
                    "14:23183137-23183158": {
                        "position": "14:23183137-23183158",
                        "results": {"mSINGS": {"status": Status.stable}}
                    }
                }
            })
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testRead(self):
        report = MSINGSReport(self.tmp_in)
        self.assertEqual(
            [toDict(rec) for rec in report.samples.values()],
            [toDict(rec) for rec in self.records]
        )


def toDict(obj):
    return json.loads(
        json.dumps(obj, default=lambda o: getattr(o, '__dict__', str(o)))
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
