#!/usr/bin/env python3
#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import pysam
import tempfile
import unittest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.join(APP_DIR, "lib")
sys.path.append(BIN_DIR)
sys.path.append(LIB_DIR)

from anacore.vcf import VCFRecord
from mergeCoOccurVar import mergedRecord, setSupportingReads

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class LoggerSilencer:
    def debug(self, args):
        pass

    def info(self, args):
        pass

class SetSupportingReads(unittest.TestCase):
    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sam_path, self.tmp_bam_path, self.tmp_bam_path + ".bai"]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_sam_path = os.path.join(tmp_folder, unique_id + ".sam")
        self.tmp_bam_path = os.path.join(tmp_folder, unique_id + ".bam")
        self.ref_seq = "ggaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgcattggggtg"
        #               | | | | | |  |  |  |  |
        #               1 3 5 7 9 11 14 17 20 23
        self.reads_content = """>subtit_AAA/CAC_1_alt
ggaagccctgatcACGCCACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_2_alt
aagccctgatcACGCCACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtit_AAA/CAC_4_mixUp
ggaagccctgatcACGCCAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_5_mixDown
ggaagccctgatcACGCAACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_1_alt
ggaagccctgatcACGCCCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_2_alt
aagccctgatcACGCCCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtitClose_AA/CC_4_mixUp
ggaagccctgatcACGCCAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_5_mixDown
ggaagccctgatcACGCACATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_1_alt
ggaagccctgatcACGCCCTTCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_2_alt
aagccctgatcACGCCCTTCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtit_AAATCTC/CCTTCGG_4_mixUp
ggaagccctgatcACGCCCTTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_5_mixDown
ggaagccctgatcACGCAAATCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_1_alt
ggaagccctgatcACGCTGGAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_2_alt
aagccctgatcACGCTGGAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insertion_A/TGGAGG_4_mixUp
ggaagccctgatcACGCTGGAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_5_mixDown
ggaagccctgatcACGCAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_1_alt
ggaagccctgatcACGCTGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_2_alt
aagccctgatcACGCTGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>deletion_AAATCTC/T_4_mixUp
ggaagccctgatcACGCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_5_mixDown
ggaagccctgatcACGCAAATGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_1_alt
ggaagccctgatcACGCTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_2_alt
aagccctgatcACGCTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>delIns_AAAT/TGA_4_mixUp
ggaagccctgatcACGCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_5_mixDown
ggaagccctgatcACGCAAATGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_1_alt
ggaagccctgatcACGCGGGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_2_alt
aagccctgatcACGCGGGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insDel_AAA/GGGA_4_mixUp
ggaagccctgatcACGCGGGAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_5_mixDown
ggaagccctgatcACGCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_1_alt
ggaagccctgatcACGCCTGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_2_alt
aagccctgatcACGCCTGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>delInsNoStd_AAATCTC/CTGGG_4_mixUp
ggaagccctgatcACGCCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_5_mixDown
ggaagccctgatcACGCAAATGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca"""
        self.test_cases = [
            [
                VCFRecord("chr1", 18, "subtit_AAA/CAC", "A", ["C"]),
                VCFRecord("chr1", 20, "subtit_AAA/CAC", "A", ["C"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtit_AAA/CAC_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:17A1A103	AS:i:113	XS:i:0
subtit_AAA/CAC_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:17A105	AS:i:118	XS:i:0
subtit_AAA/CAC_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCAACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:19A103	AS:i:118	XS:i:0
subtit_AAA/CAC_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtit_AAA/CAC_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:15A1A103	AS:i:111	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "subtitClose_AA/CC", "A", ["C"]),
                VCFRecord("chr1", 19, "subtitClose_AA/CC", "A", ["C"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtitClose_AA/CC_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:17A0A104	AS:i:113	XS:i:0
subtitClose_AA/CC_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:17A105	AS:i:118	XS:i:0
subtitClose_AA/CC_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCACATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:18A104	AS:i:118	XS:i:0
subtitClose_AA/CC_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtitClose_AA/CC_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:15A0A104	AS:i:111	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "subtit_AAATCTC/CCTTCGG", "AAA", ["CCT"]),
                VCFRecord("chr1", 23, "subtit_AAATCTC/CCTTCGG", "TC", ["GG"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtit_AAATCTC/CCTTCGG_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCTTCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17A0A0A2T0C99	AS:i:99	XS:i:0
subtit_AAATCTC/CCTTCGG_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCTTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17A0A0A103	AS:i:108	XS:i:0
subtit_AAATCTC/CCTTCGG_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCAAATCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:22T0C99	AS:i:113	XS:i:0
subtit_AAATCTC/CCTTCGG_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtit_AAATCTC/CCTTCGG_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCCTTCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15A0A0A2T0C99	AS:i:99	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "insertion_A/TGGAGG", "-", ["TGG"]),
                VCFRecord("chr1", 19, "insertion_A/TGGAGG", "-", ["GG"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insertion_A/TGGAGG_1_alt	0	chr1	1	60	17M3I1M2I105M	*	0	0	GGAAGCCCTGATCACGCTGGAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:123	AS:i:107	XS:i:0
insertion_A/TGGAGG_4_mixUp	0	chr1	1	60	17M3I106M	*	0	0	GGAAGCCCTGATCACGCTGGAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
insertion_A/TGGAGG_5_mixDown	0	chr1	1	60	18M2I105M	*	0	0	GGAAGCCCTGATCACGCAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:123	AS:i:115	XS:i:0
insertion_A/TGGAGG_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insertion_A/TGGAGG_2_alt	0	chr1	3	60	15M3I1M2I105M	*	0	0	AAGCCCTGATCACGCTGGAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:121	AS:i:105	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "deletion_AAATCTC/T", "AAA", ["-"]),
                VCFRecord("chr1", 22, "deletion_AAATCTC/T", "CTC", ["-"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
deletion_AAATCTC/T_1_alt	0	chr1	1	60	17M3D1M3D99M	*	0	0	GGAAGCCCTGATCACGCTGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17^AAA1^CTC99	AS:i:100	XS:i:0
deletion_AAATCTC/T_4_mixUp	0	chr1	1	60	17M3D103M	*	0	0	GGAAGCCCTGATCACGCTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^AAA103	AS:i:111	XS:i:0
deletion_AAATCTC/T_5_mixDown	0	chr1	1	60	21M3D99M	*	0	0	GGAAGCCCTGATCACGCAAATGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:21^CTC99	AS:i:111	XS:i:0
deletion_AAATCTC/T_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
deletion_AAATCTC/T_2_alt	0	chr1	3	60	15M3D1M3D99M	*	0	0	AAGCCCTGATCACGCTGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15^AAA1CTC99	AS:i:99	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "delIns_AAAT/TGA", "AAA", ["-"]),
                VCFRecord("chr1", 22, "delIns_AAAT/TGA", "-", ["GA"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
delIns_AAAT/TGA_1_alt	0	chr1	1	60	17M3D1M2I102M	*	0	0	GGAAGCCCTGATCACGCTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17^AAA103	AS:i:105	XS:i:0
delIns_AAAT/TGA_4_mixUp	0	chr1	1	60	17M3D103M	*	0	0	GGAAGCCCTGATCACGCTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^AAA103	AS:i:111	XS:i:0
delIns_AAAT/TGA_5_mixDown	0	chr1	1	60	21M2I102M	*	0	0	GGAAGCCCTGATCACGCAAATGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:123	AS:i:115	XS:i:0
delIns_AAAT/TGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
delIns_AAAT/TGA_2_alt	0	chr1	3	60	15M3D1M2I102M	*	0	0	AAGCCCTGATCACGCTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15^AAA103	AS:i:103	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "insDel_AAA/GGGA", "-", ["GGG"]),
                VCFRecord("chr1", 19, "insDel_AAA/GGGA", "AA", ["-"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insDel_AAA/GGGA_1_alt	0	chr1	1	60	17M3I1M2D103M	*	0	0	GGAAGCCCTGATCACGCGGGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:18^AA103	AS:i:106	XS:i:0
insDel_AAA/GGGA_4_mixUp	0	chr1	1	60	17M3I106M	*	0	0	GGAAGCCCTGATCACGCGGGAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
insDel_AAA/GGGA_5_mixDown	0	chr1	1	60	18M2D103M	*	0	0	GGAAGCCCTGATCACGCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:18^AA103	AS:i:113	XS:i:0
insDel_AAA/GGGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insDel_AAA/GGGA_2_alt	0	chr1	3	60	15M3I1M2D103M	*	0	0	AAGCCCTGATCACGCGGGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:16^AA103	AS:i:104	XS:i:0"""
            ],
            # [
            #     VCFRecord("chr1", 18, "delInsNoStd_AAATCTC/CTGGG", "AAA", ["C"]),
            #     VCFRecord("chr1", 22, "delInsNoStd_AAATCTC/CTGGG", "-", ["GGG"]),
            #     """"""
            # ]
            # VCFRecord("chr1", 18, "insDelNoStd_AAATCTC/T", "AAA", "-")
        ]

    def testSetSupportingReads(self):
        for first, second, aln_content in self.test_cases:
            # Write BAM
            with open(self.tmp_sam_path, "w") as FH_sam:
                FH_sam.write(aln_content)
            with pysam.AlignmentFile(self.tmp_sam_path) as FH_sam:
                with pysam.AlignmentFile(self.tmp_bam_path, "wb", template=FH_sam) as FH_bam:
                    for rec in FH_sam:
                        FH_bam.write(rec)
            pysam.index(self.tmp_bam_path)
            # Eval
            with pysam.AlignmentFile(self.tmp_bam_path) as FH_aln:
                first.end = int(first.refEnd())
                second.start = int(second.refStart())
                first.isIns = first.isInsertion()
                second.isIns = second.isInsertion()
                print("{}/{}\t{}/{}".format(first.ref, first.alt[0], second.ref, second.alt[0]))
                setSupportingReads(first, second, FH_aln, LoggerSilencer())
                # Check supporting first
                expected = sorted([
                    "{}_{}".format(first.id, curr_suffix) for curr_suffix in ["1_alt", "2_alt", "4_mixUp"]
                ])
                self.assertEqual(
                    sorted(first.supporting_reads),
                    expected
                )
                # Check supporting first
                expected = sorted([
                    "{}_{}".format(second.id, curr_suffix) for curr_suffix in ["1_alt", "2_alt", "5_mixDown"]
                ])
                self.assertEqual(
                    sorted(second.supporting_reads),
                    expected
                )


class MergeCoOccurVar(unittest.TestCase):
    def setUp(self):
        self.ref_seq = "ACGCAAATCTCGGCATGCCGATT"
        #               | | | | | |  |  |  |  |
        #               1 3 5 7 9 11 14 17 20 23
        self.variant_1 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            "artificial_1",  # id
            None,  # ref
            None,  # alt
            10,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.05]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [10], "DP": 100},
                "splB": {"AD": [40], "DP": 4900},
            }
        )
        self.variant_2 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            30,  # qual
            ["PASS"],  # filter
            {"AF": [0.06]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )
        self.expected_merge = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            20,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=A/T", "chr1:20=G/C"]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )

    def testMergedRecord_1_substit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "A"
        self.variant_1.alt = ["T"]
        # Variant 2
        self.variant_2.pos = 20
        self.variant_2.ref = "G"
        self.variant_2.alt = ["C"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCTCGGCATGCCG"
        self.expected_merge.alt = ["TAATCTCGGCATGCCC"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=A/T", "chr1:20=G/C"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_2_largeSubstit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["TGCA"]
        # Variant 2
        self.variant_2.pos = 10
        self.variant_2.ref = "TC"
        self.variant_2.alt = ["GG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCTC"
        self.expected_merge.alt = ["TGCACGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/TGCA", "chr1:10=TC/GG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_3_largeCloseSubstit(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["TGCA"]
        # Variant 2
        self.variant_2.pos = 9
        self.variant_2.ref = "CT"
        self.variant_2.alt = ["GG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATCT"
        self.expected_merge.alt = ["TGCAGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/TGCA", "chr1:9=CT/GG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_4_delIns(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["-"]
        # Variant 2
        self.variant_2.pos = 10
        self.variant_2.ref = "-"
        self.variant_2.alt = ["GGCATCT"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATC"
        self.expected_merge.alt = ["CGGCATCT"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/-", "chr1:10=-/GGCATCT"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_5_coDelIns(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "AAAT"
        self.variant_1.alt = ["-"]
        # Variant 2
        self.variant_2.pos = 9
        self.variant_2.ref = "-"
        self.variant_2.alt = ["AGG"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAAT"
        self.expected_merge.alt = ["AGG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=AAAT/-", "chr1:9=-/AGG"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_6_insDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 7
        self.variant_2.ref = "ATC"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAATC"
        self.expected_merge.alt = ["GTGTGAA"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:7=ATC/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_7_closeInsDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 6
        self.variant_2.ref = "AA"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AAA"
        self.expected_merge.alt = ["GTGTGA"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:6=AA/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    def testMergedRecord_8_coInsDel(self):
        # Variant 1
        self.variant_1.pos = 5
        self.variant_1.ref = "-"
        self.variant_1.alt = ["GTGTG"]
        # Variant 2
        self.variant_2.pos = 5
        self.variant_2.ref = "AA"
        self.variant_2.alt = ["-"]
        # Expected merge
        self.expected_merge.pos = 5
        self.expected_merge.ref = "AA"
        self.expected_merge.alt = ["GTGTG"]
        self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=-/GTGTG", "chr1:5=AA/-"]}
        # Eval
        observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
        self.assertEqual(
            strVariant(observed_merge),
            strVariant(self.expected_merge)
        )

    # def testMergedRecord_9_delInsOtherSyntax(self):
    #     # Variant 1
    #     self.variant_1.pos = 4
    #     self.variant_1.ref = "CAAAT"
    #     self.variant_1.alt = ["C"]
    #     # Variant 2
    #     self.variant_2.pos = 9
    #     self.variant_2.ref = "C"
    #     self.variant_2.alt = ["CGCATCT"]
    #     # Expected merge
    #     self.expected_merge.pos = 5
    #     self.expected_merge.ref = "AAATC"
    #     self.expected_merge.alt = ["CGGCATCT"]
    #     self.expected_merge.info = {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:4=CAAAT/-", "chr1:9=C/CGGCATCT"]}
    #     # Eval
    #     observed_merge = mergedRecord(self.variant_1, self.variant_2, self.ref_seq)
    #     self.assertEqual(
    #         strVariant(observed_merge),
    #         strVariant(self.expected_merge)
    #     )


def strVariant(var):
    info = []
    for info_key, info_val in sorted(var.info.items()):
        info.append(
            "{}:{}".format(info_key, info_val)
        )
    samples = []
    for spl_name, spl_val in sorted(var.samples.items()):
        for key, val in sorted(spl_val.items()):
            samples.append(
                "{}:{}:{}".format(spl_name, key, val)
            )
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        var.getName(),
        var.id,
        var.qual,
        sorted(var.filter),
        info,
        sorted(var.format),
        samples
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
