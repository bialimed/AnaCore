#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFByNames(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_ids = os.path.join(tmp_folder, unique_id + "_in_ids.txt")
        self.tmp_in_variants = os.path.join(tmp_folder, unique_id + "_in.vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")
        self.tmp_list = [self.tmp_in_ids, self.tmp_in_variants, self.tmp_output]

        # Command
        self.cmd = [
            "filterVCFByNames.py",
            "--input-ids", self.tmp_in_ids,
            "--input-variants", self.tmp_in_variants,
            "--output-variants", self.tmp_output
        ]

        # Create input files
        with open(self.tmp_in_ids, "w") as FH_out:
            FH_out.write("""21:5030088=C/T
21:5030240=AC/A
21:5030241=C/AGT
21:5030319=C/G/T
22:18457=G/T""")
        with open(self.tmp_in_variants, "w") as FH_out:
            FH_out.write("""##fileformat=VCFv4.1
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
21	5030088	rs1455320509	C	T	.	.	AA=C
21	5030218	rs1416749598	TA	T	.	.	AA=A
21	5030240	rs1185110752	AC	A	.	.	AA=C
21	5030241	.	C	AGT	.	.	.
21	5030319	rs1231507777	C	G,T	.	.	AA=C
21	5030625	rs1380012774	GCAT	G	.	.	AA=CAT""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in self.tmp_list:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testSelect(self):
        # Stringency strict
        with self.assertRaises(Exception) as context:
            subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [elt for elt in self.cmd] + [
            "--stringency", "lenient"
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            "##fileformat=VCFv4.1",
            '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
            '##FILTER=<ID=selID,Description="Select variants by list of names">',
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
            "21	5030088	rs1455320509	C	T	.	.	AA=C",
            "21	5030240	rs1185110752	AC	A	.	.	AA=C",
            "21	5030241	.	C	AGT	.	.	.",
            "21	5030319	rs1231507777	C	G,T	.	.	AA=C",
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )

    def testRemove(self):
        # Stringency strict
        cmd = [elt for elt in self.cmd] + [
            "--mode", "remove"
        ]
        with self.assertRaises(Exception) as context:
            subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [elt for elt in self.cmd] + [
            "--stringency", "lenient",
            "--mode", "remove"
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            "##fileformat=VCFv4.1",
            '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
            '##FILTER=<ID=remvID,Description="Remove variants by list of names">',
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
            "21	5030218	rs1416749598	TA	T	.	.	AA=A",
            "21	5030625	rs1380012774	GCAT	G	.	.	AA=CAT"
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
