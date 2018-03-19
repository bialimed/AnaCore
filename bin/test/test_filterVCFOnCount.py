#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT-O
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
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.join(APP_DIR, "lib")
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, VCFRecord

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFPrimers(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_variants = os.path.join( tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join( tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "filterVCFOnCount.py",
            "--AF-threshold", str(0.1),
            "--DP-threshold", str(500),
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Create VCF
        vcf_content = """##fileformat=VCFv4.0
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	10	1	A	G	.	.	AF=0.1;DP=1000	DP:AF	500:0.1	500:0.1
1	10	2	A	G	.	.	AF=0.1;DP=1000	DP:AF	500:0.01	500:0.01
1	10	3	A	G	.	.	AF=0.1;DP=2	DP:AF	1:0.1	1:0.1
1	10	4	A	G	.	.	AF=0.01;DP=2	DP:AF	1:0.01	1:0.01
1	10	5	A	G	.	.	AF=0.1;DP=501	DP:AF	500:0.1	1:0.01
1	10	6	A	G	.	.	AF=0.01;DP=501	DP:AF	500:0.01	1:0.1
1	10	7	A	G	.	.	AF=0.1;DP=501	DP:AF	500:0.1	1:0.1
1	10	8	A	G	.	.	AF=0.01;DP=501	DP:AF	500:0.01	1:0.01
1	10	9	A	G	.	.	AF=0.01;DP=501	DP:AF	1:0.1	500:0.01
1	10	10	A	G	.	.	AF=0.1;DP=501	DP:AF	1:0.01	500:0.1
1	10	11	A	G	.	.	AF=0.1;DP=501	DP:AF	1:0.1	500:0.1
1	10	12	A	G	.	.	AF=0.01;DP=501	DP:AF	1:0.01	500:0.01
1	10	13	A	G	.	.	AF=0.055;DP=2	DP:AF	1:0.01	1:0.1
1	10	14	A	G	.	.	AF=0.55;DP=1000	DP:AF	500:0.01	500:0.1
1	10	15	A	G	.	.	AF=0.055;DP=2	DP:AF	1:0.1	1:0.01
1	10	16	A	G	.	.	AF=0.55;DP=1000	DP:AF	500:0.1	500:0.01"""
        with open(self.tmp_variants, "w") as FH_variants:
            FH_variants.write( vcf_content )

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = """##fileformat=VCFv4.0
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=incomplete,Description="The variant has a sufficient frequency in one sample (AF>=10.0% and DP>=500) but the depth in second sample it has an insufficiant depth (the AF cannot be take into account).">
##FILTER=<ID=invalid,Description="The variant has an insufficient frequency in one sample (AF<10.0% and DP>=500) but the depth in second sample it has an insufficiant depth (the AF cannot be take into account).">
##FILTER=<ID=libSpe,Description="The variant has a frequency lower than 10.0% in one of the sample and depth is superior than 500 in all the samples.">
##FILTER=<ID=lowAF,Description="The variant has a frequency lower than 10.0% in all the samples and depth is superior than 500 in all the samples.">
##FILTER=<ID=lowDP,Description="The variant has a depth lower than 500 in all the samples.">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	10	1	A	G	.	PASS	AF=0.1;DP=1000	DP:AF	500:0.1	500:0.1
1	10	2	A	G	.	lowAF	AF=0.1;DP=1000	DP:AF	500:0.01	500:0.01
1	10	3	A	G	.	lowDP	AF=0.1;DP=2	DP:AF	1:0.1	1:0.1
1	10	4	A	G	.	lowDP	AF=0.01;DP=2	DP:AF	1:0.01	1:0.01
1	10	5	A	G	.	incomplete	AF=0.1;DP=501	DP:AF	500:0.1	1:0.01
1	10	6	A	G	.	invalid	AF=0.01;DP=501	DP:AF	500:0.01	1:0.1
1	10	7	A	G	.	incomplete	AF=0.1;DP=501	DP:AF	500:0.1	1:0.1
1	10	8	A	G	.	invalid	AF=0.01;DP=501	DP:AF	500:0.01	1:0.01
1	10	9	A	G	.	invalid	AF=0.01;DP=501	DP:AF	1:0.1	500:0.01
1	10	10	A	G	.	incomplete	AF=0.1;DP=501	DP:AF	1:0.01	500:0.1
1	10	11	A	G	.	incomplete	AF=0.1;DP=501	DP:AF	1:0.1	500:0.1
1	10	12	A	G	.	invalid	AF=0.01;DP=501	DP:AF	1:0.01	500:0.01
1	10	13	A	G	.	lowDP	AF=0.055;DP=2	DP:AF	1:0.01	1:0.1
1	10	14	A	G	.	libSpe	AF=0.55;DP=1000	DP:AF	500:0.01	500:0.1
1	10	15	A	G	.	lowDP	AF=0.055;DP=2	DP:AF	1:0.1	1:0.01
1	10	16	A	G	.	libSpe	AF=0.55;DP=1000	DP:AF	500:0.1	500:0.01""".split("\n")
        observed = list()
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip().split("\n")
        self.assertEqual(
            expected[1:], # skip file format version
            observed[1:] # skip file format version
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
