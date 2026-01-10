#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.illumina.base import getInfFromSeqDesc, getInfFromSeqID, getPlatformFromSerialNumber


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
        observed = getInfFromSeqDesc("2:N:0:")


class TestGetInfFromSeqID(unittest.TestCase):
    def test(self):
        # Whithout UMI
        expected = {
            "sequencer_id": "NDX550421_RUO",
            "run_id": "15",
            "flowcell_id": "HY53NBGXC",
            "lane_id": 1,
            "tile_id": 11101,
            "x_pos": 9518,
            "y_pos": 1037,
            "umi": None
        }
        observed = getInfFromSeqID("NDX550421_RUO:15:HY53NBGXC:1:11101:9518:1037")
        self.assertEqual(expected, observed)
        # With UMI
        expected = {
            "sequencer_id": "NDX550421_RUO",
            "run_id": "15",
            "flowcell_id": "HY53NBGXC",
            "lane_id": 1,
            "tile_id": 11101,
            "x_pos": 9518,
            "y_pos": 1037,
            "umi": "TTGCCNG+CCTCATA"
        }
        observed = getInfFromSeqID("NDX550421_RUO:15:HY53NBGXC:1:11101:9518:1037:TTGCCNG+CCTCATA")
        self.assertEqual(expected, observed)


class TestGetPlatformFromSerialNumber(unittest.TestCase):
    def test(self):
        self.assertEqual(getPlatformFromSerialNumber("ML-P2-01"), "MiniSeq")
        self.assertEqual(getPlatformFromSerialNumber("M70265"), "MiSeq")
        self.assertEqual(getPlatformFromSerialNumber("NDX550421"), "NextSeq")
        self.assertEqual(getPlatformFromSerialNumber("NB551452"), "NextSeq")
        self.assertEqual(getPlatformFromSerialNumber("NDX550500"), "NextSeq")
        self.assertEqual(getPlatformFromSerialNumber("D00154"), "HiSeq")
        self.assertEqual(getPlatformFromSerialNumber("A01789"), "NovaSeq")
        self.assertEqual(getPlatformFromSerialNumber("SH00093"), "MiSeq_i100")
        self.assertEqual(getPlatformFromSerialNumber("VH00159"), "NextSeq_2000")
        self.assertEqual(getPlatformFromSerialNumber("LH00101"), "NovaSeq_X")


if __name__ == "__main__":
    unittest.main()
