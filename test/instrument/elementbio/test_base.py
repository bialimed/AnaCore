#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.elementbio.base import getInfFromSeqDesc, getInfFromSeqID


class TestGetInfFromSeqDesc(unittest.TestCase):
    def test(self):
        # Standard
        expected = {
            "reads_phases": 1,
            "is_kept": False,
            "control_bits": 18,
            "barcode": "ATCACG+ATTA"
        }
        observed = getInfFromSeqDesc("1:Y:18:ATCACG+ATTA")
        self.assertEqual(expected, observed)
        # Without barcode
        expected = {
            "reads_phases": 2,
            "is_kept": True,
            "control_bits": None,
            "barcode": None
        }
        observed = getInfFromSeqDesc("2:N:0:1")


class TestGetInfFromSeqID(unittest.TestCase):
    def test(self):
        # Whithout UMI
        expected = {
            "sequencer_name": "AVITI1",
            "run_name": "SIDEB_ADA_2510541191",
            "flowcell_id": "2510541191",
            "lane_id": 1,
            "tile_id": 10102,
            "x_pos": 144,
            "y_pos": 37,
            "umi": None
        }
        observed = getInfFromSeqID("AVITI1:SIDEB_ADA_2510541191:2510541191:1:10102:0144:0037")
        self.assertEqual(expected, observed)
        # With UMI
        expected = {
            "sequencer_name": "AVITI3",
            "run_name": "run_tera_3",
            "flowcell_id": "2510541192",
            "lane_id": 2,
            "tile_id": 10102,
            "x_pos": 144,
            "y_pos": 37,
            "umi": "ATTGC+GGTA"
        }
        observed = getInfFromSeqID("AVITI3:run_tera_3:2510541192:2:10102:0144:0037:ATTGC+GGTA")
        self.assertEqual(expected, observed)


if __name__ == "__main__":
    unittest.main()
