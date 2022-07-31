#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU-Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import json
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.msi.base import Status, toDict
from anacore.msi.reportIO import ReportIO


########################################################################
#
# FUNCTIONS
#
########################################################################
REPORT_CONTENT = """[
    {
        "loci": {
            "1:208454-208466": {
                "name": null,
                "position": "1:208454-208466",
                "results": {
                    "mSINGSUp": {
                        "data": {
                            "lengths": {
                                "ct_by_len": {
                                    "10": 140,
                                    "11": 100,
                                    "12": 180
                                },
                                "mode": "reads"
                            }
                        },
                        "score": 0.84,
                        "status": "MSS"
                    },
                    "SVC": {
                        "data": {
                            "lengths": {
                                "ct_by_len": {
                                    "10": 140,
                                    "11": 100,
                                    "12": 180
                                },
                                "mode": "reads"
                            }
                        },
                        "score": 0.99,
                        "status": "MSI"
                    }
                }
            },
            "8:578711-578723": {
                "name": "B12",
                "position": "8:578711-578723",
                "results": {
                    "mSINGSUp": {
                        "data": {
                            "lengths": {
                                "ct_by_len": {
                                    "10": 20,
                                    "11": 15,
                                    "12": 30
                                },
                                "mode": "reads"
                            }
                        },
                        "score": null,
                        "status": "Undetermined"
                    },
                    "SVC": {
                        "data": {
                            "lengths": {
                                "ct_by_len": {
                                    "10": 20,
                                    "11": 15,
                                    "12": 30
                                },
                                "mode": "reads"
                            }
                        },
                        "score": 0.99,
                        "status": "MSI"
                    }
                }
            }
        },
        "name": "splA",
        "results": {
            "SVC": {
                "method": "SVC",
                "param": {"aggregation_method": "instability count", "instability_threshold": 0.33, "min_voting_loci": 1},
                "score": 0.99,
                "status": "MSI",
                "version": "1.0.0"
            },
            "mSINGSUp": {
                "method": "mSINGSUp",
                "param": {"aggregation_method": "instability count", "instability_threshold": 0.4, "min_voting_loci": 1},
                "score": 0.84,
                "status": "MSS",
                "version": "1.0.0"
            }
        }
    }
]"""


class TestReportIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.json")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.json")

        with open(self.tmp_in, "w") as writer:
            writer.write(REPORT_CONTENT)

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGetIncompleteModels(self):
        with open(self.tmp_out, "w") as writer:
            writer.write(REPORT_CONTENT.replace('"SVC"', '"model"'))
        self.assertEqual(
            sorted(
                ReportIO.getIncompleteModels(self.tmp_out, 2),
                key=lambda elt: (elt["locus_id"], elt["status"])
            ),
            [
                {"locus_id": "1:208454-208466", "locus_name": None, "status": "MSI", "support": 1},
                {"locus_id": "1:208454-208466", "locus_name": None, "status": "MSS", "support": 0},
                {"locus_id": "8:578711-578723", "locus_name": "B12", "status": "MSI", "support": 1},
                {"locus_id": "8:578711-578723", "locus_name": "B12", "status": "MSS", "support": 0}
            ]
        )
        self.assertEqual(
            sorted(
                ReportIO.getIncompleteModels(self.tmp_out, 1),
                key=lambda elt: (elt["locus_id"], elt["status"])
            ),
            [
                {"locus_id": "1:208454-208466", "locus_name": None, "status": "MSS", "support": 0},
                {"locus_id": "8:578711-578723", "locus_name": "B12", "status": "MSS", "support": 0}
            ]
        )

    def testParse(self):
        report = ReportIO.parse(self.tmp_in)
        self.assertEqual(len(report), 1)
        self.assertEqual(report[0].results["SVC"].status, Status.unstable)
        self.assertEqual(report[0].results["mSINGSUp"].status, Status.stable)
        self.assertEqual(len(report[0].loci), 2)
        self.assertEqual(report[0].loci["1:208454-208466"].name, None)
        self.assertEqual(report[0].loci["8:578711-578723"].name, "B12")
        self.assertEqual(report[0].loci["1:208454-208466"].results["SVC"].data["lengths"].getCount(), 420)
        self.assertEqual(report[0].loci["8:578711-578723"].results["mSINGSUp"].data["lengths"].getCount(), 65)

    def testWrite(self):
        report = ReportIO.parse(self.tmp_in)
        ReportIO.write(report, self.tmp_out)
        observed = None
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        self.assertEqual(observed, json.loads(REPORT_CONTENT))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
