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
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.msi.base import Status, toDict
from anacore.msi.locus import Locus, LocusDataDistrib, LocusRes


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestLocus(unittest.TestCase):
    def testDelResult(self):
        locus = Locus.fromDict({
            "position": "4:54732045",
            "results": {"SVC": {"status": "MSS", "score": 0.99}}
        })
        self.assertEqual(locus.results["SVC"].__class__, LocusRes)
        locus.delResult("SVC")
        self.assertEqual(locus.results, {})

    def testFromDict(self):
        locus = Locus.fromDict({
            "position": "4:54732045",
            "results": {"SVC": {"status": "MSS", "score": 0.99}}
        })
        self.assertEqual(locus.name, None)
        self.assertEqual(locus.position, "4:54732045")
        self.assertEqual(locus.results["SVC"].__class__, LocusRes)
        self.assertEqual(locus.results["SVC"].status, Status.stable)

    def testToDict(self):
        locus = Locus.fromDict({
            "position": "4:54732045",
            "results": {"SVC": {"status": "MSS", "score": 0.99}}
        })
        locus_hash = json.dumps(locus, sort_keys=True, default=toDict)
        self.assertEqual(
            locus_hash,
            '{"name": null, "position": "4:54732045", "results": {"SVC": {"data": {}, "score": 0.99, "status": "MSS"}}}'
        )


class TestLocusDataDistrib(unittest.TestCase):
    def testGetCount(self):
        self.assertEqual(LocusDataDistrib().getCount(), 0)
        self.assertEqual(LocusDataDistrib({1: 5, 8: 0, 11: 8}).getCount(), 13)

    def testGetDenseCount(self):
        self.assertEqual(
            LocusDataDistrib().getDenseCount(0, 4),
            [0, 0, 0, 0, 0]
        )
        self.assertEqual(
            LocusDataDistrib({1: 5, 8: 0, 11: 8}).getDenseCount(0, 4),
            [0, 5, 0, 0, 0]
        )
        self.assertEqual(
            LocusDataDistrib({1: 5, 8: 0, 11: 8}).getDenseCount(),
            [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8]
        )

    def testGetDensePrct(self):
        self.assertEqual(
            LocusDataDistrib().getDensePrct(0, 4),
            [None, None, None, None, None]
        )
        self.assertEqual(
            LocusDataDistrib({1: 5, 8: 0, 11: 8}).getDensePrct(0, 4),
            [0, (5 * 100) / 13, 0, 0, 0]
        )
        self.assertEqual(
            LocusDataDistrib({1: 5, 8: 0, 11: 8}).getDensePrct(),
            [(5 * 100) / 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, (8 * 100) / 13]
        )

    def testGetMaxLength(self):
        with self.assertRaises(Exception):
            LocusDataDistrib().getMaxLength()

        self.assertEqual(
            LocusDataDistrib({1: 5, 8: 0, 11: 8, 20: 0}).getMaxLength(),
            11
        )

    def testGetMinLength(self):
        with self.assertRaises(Exception):
            LocusDataDistrib().getMinLength()

        self.assertEqual(
            LocusDataDistrib({0: 0, 1: 5, 8: 0, 11: 8}).getMinLength(),
            1
        )

    def testItems(self):
        self.assertEqual(
            list(LocusDataDistrib({0: 0, 1: 5, 8: 0, 11: 8}).items()),
            [(0, 0), (1, 5), (8, 0), (11, 8)],
        )

    def testToDict(self):
        distrib = LocusDataDistrib({1: 5, 8: 0, 11: 8}, "fragments")
        distrib_hash = json.dumps(distrib, sort_keys=True, default=toDict)
        self.assertEqual(
            distrib_hash,
            '{"ct_by_len": {"1": 5, "8": 0, "11": 8}, "mode": "fragments"}'
        )


class TestLocusRes(unittest.TestCase):
    def testFromDict(self):
        # Need to convert LocusDataDistrib
        locus_res = LocusRes.fromDict({
            "status": Status.unstable,
            "data": {"lengths": {"ct_by_len": {1: 5, 8: 0, 11: 8}}}
        })
        self.assertEqual(locus_res.status, Status.unstable)
        self.assertEqual(locus_res.data["lengths"].__class__, LocusDataDistrib)
        self.assertEqual(locus_res.data["lengths"].ct_by_len, {1: 5, 8: 0, 11: 8})
        # No need to convert LocusDataDistrib
        locus_res = LocusRes.fromDict({
            "status": Status.stable,
            "data": {"lengths": LocusDataDistrib({1: 5, 8: 0, "11": 8})}
        })
        self.assertEqual(locus_res.status, Status.stable)
        self.assertEqual(locus_res.data["lengths"].__class__, LocusDataDistrib)
        self.assertEqual(locus_res.data["lengths"].ct_by_len, {1: 5, 8: 0, 11: 8})

    def testToDict(self):
        locus_res = LocusRes.fromDict({
            "status": Status.unstable,
            "data": {"lengths": {"ct_by_len": {1: 5, 8: 0, 11: 8}}}
        })
        locus_hash = json.dumps(locus_res, sort_keys=True, default=toDict)
        self.assertEqual(
            locus_hash,
            '{"data": {"lengths": {"ct_by_len": {"1": 5, "8": 0, "11": 8}, "mode": "reads"}}, "score": null, "status": "MSI"}'
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
