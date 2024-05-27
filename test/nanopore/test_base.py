#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.nanopore.base import getInfFromSeqDesc


class TestGetInfFromSeqDesc(unittest.TestCase):
    def test(self):
        # Standard
        expected = {
            "run_id": "2cf9e7c4dd1888d9a1090ebc394dd5b7cb2fe3f0",
            "read_number": 126,
            "channel_id": 181,
            "start_time": "2017-10-17T06:43:19Z",
            "barcode_id": "barcode01",
            "parent_read_id": None
        }
        observed = getInfFromSeqDesc("runid=2cf9e7c4dd1888d9a1090ebc394dd5b7cb2fe3f0 read=126 ch=181 start_time=2017-10-17T06:43:19Z barcode=barcode01")
        self.assertEqual(expected, observed)
        # Without barcode and with additionnal info
        expected = {
            "run_id": "0a82096db062bcf8041358b407ab5e0aa5de7138",
            "read_number": 23606,
            "channel_id": 186,
            "start_time": "2019-03-11T16:16:33Z",
            "barcode_id": None,
            "parent_read_id": None,
            "flow_cell_id": None,
            "protocol_group_id": "190311_ang_1",
            "sample_id": "190311_ang_1"
        }
        observed = getInfFromSeqDesc("runid=0a82096db062bcf8041358b407ab5e0aa5de7138 read=23606 ch=186 start_time=2019-03-11T16:16:33Z flow_cell_id= protocol_group_id=190311_ang_1 sample_id=190311_ang_1")
        self.assertEqual(expected, observed)
        # Without barcode
        expected = {
            "run_id": "17838b1d08f30a031bf60afabb146a8b0fba7486",
            "read_number": 12217,
            "channel_id": 492,
            "start_time": "2017-07-04T06:42:43Z",
            "barcode_id": None,
            "parent_read_id": None
        }
        observed = getInfFromSeqDesc("runid=17838b1d08f30a031bf60afabb146a8b0fba7486 read=12217 ch=492 start_time=2017-07-04T06:42:43Z")
        self.assertEqual(expected, observed)


if __name__ == "__main__":
    unittest.main()
