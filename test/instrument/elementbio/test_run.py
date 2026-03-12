#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.elementbio.run import Run


class TestRun(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in_folder = os.path.join(tmp_folder, unique_id + "_run")
        os.makedirs(self.tmp_in_folder)

    def tearDown(self):
        if os.path.exists(self.tmp_in_folder):
            for filename in os.listdir(self.tmp_in_folder):
                os.remove(os.path.join(self.tmp_in_folder, filename))
            os.rmdir(self.tmp_in_folder)

    def testIsCopied(self):
        run = Run(self.tmp_in_folder)
        with open(os.path.join(self.tmp_in_folder, "AvitiRunStats.json"), "w") as writer:
            writer.write("")
        self.assertFalse(run.isCopied())
        with open(os.path.join(self.tmp_in_folder, "RunUploaded.json"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isCopied())

    def testIsEnded(self):
        run = Run(self.tmp_in_folder)
        with open(os.path.join(self.tmp_in_folder, "AvitiRunStats.json"), "w") as writer:
            writer.write("")
        self.assertFalse(run.isEnded())
        with open(os.path.join(self.tmp_in_folder, "RunUploaded.json"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isEnded())

    def testIsSequenced(self):
        run = Run(self.tmp_in_folder)
        self.assertFalse(run.isSequenced())
        with open(os.path.join(self.tmp_in_folder, "AvitiRunStats.json"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isSequenced())


if __name__ == "__main__":
    unittest.main()
