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
import sys
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.join(APP_DIR, "lib")
sys.path.append(LIB_DIR)

from anacore.bed import BEDIO, BEDRecord
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr
from anacore.sequenceIO import FastaIO, Sequence

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class StandardizeVCF(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_sequences = os.path.join(tmp_folder, unique_id + ".fasta")
        self.tmp_faidx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "standardizeVCF.py",
            "--trace-unstandard",
            "--input-reference", self.tmp_sequences,
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Create fasta
        with FastaIO(self.tmp_sequences, "w") as FH_seq:
            # Repeats:                                       ****....            ...***
            # Region:                                 |----|        |------------|         |------|
            FH_seq.write(Sequence("artificial_chr1", "CTCAGTCATGTATGTATGTGCTCACAAAGTAGTAGATCATGGCAC"))
            #                                         123456789| | | | | | | | | | | | | | | | | |
            #                                                  10| 14| 18| 22| 26| 30| 34| 38| 42|
            #                                                    12  16  20  24  28  32  36  40  44
            FH_seq.write(Sequence("artificial_chr2", "CGATNNNCGAT"))
            #                                         123456789|
            #                                                  10

        # Create faidx
        with open(self.tmp_faidx, "w") as FH_fai:
            FH_fai.write("""artificial_chr1	45	17	45	46
artificial_chr2	11	80	11	12""")

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {
                "expected": HeaderInfoAttr("expected", "Standardized version of {chrom}:{pos}={ref}/{alt}.", type="String", number="."),
                "ANN": HeaderInfoAttr("ANN", "Annotation of variants Format: Allele|Annotation_id|Alt_allele_idx", type="String", number="."),
                "expectedANN": HeaderInfoAttr("expectedANN", "Standardized version of annotations Format: Allele|Annotation_id|Alt_allele_idx", type="String", number=".")
            }
            FH_var.writeHeader()
            self.variants = [
                # Substit single nt
                VCFRecord("artificial_chr1", 14, "sub_01", "G", ["T"], None, None, {
                    "expected": ["artificial_chr1:14=G/T"],
                    "ANN": ["T|ann_1|0", "T|ann_2|0", "A|ann_3|"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 19, "sub_02", "T", ["A", "C"], None, None, {
                    "expected": ["artificial_chr1:19=T/A", "artificial_chr1:19=T/C"],
                    "ANN": ["A|ann_1|0", "A|ann_2|0", "T|ann_3|"],
                    "expectedANN": ["A|ann_1|0", "A|ann_2|0"]
                }),
                # Substit multi nt
                VCFRecord("artificial_chr1", 7, "sub_03", "CATGTATG", ["GTACCCGC"], None, None, {
                    "expected": ["artificial_chr1:7=CATGTATG/GTACCCGC"],
                    "ANN": ["GTACCCGC|ann_1|0", "GTACCCGC|ann_2|0", "GTGT|ann_3|"],
                    "expectedANN": ["GTACCCGC|ann_1|0", "GTACCCGC|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 11, "sub_04", "TATGTATG", ["GTACCCGC", "GTACCCAA"], None, None, {
                    "expected": ["artificial_chr1:11=TATGTATG/GTACCCGC", "artificial_chr1:11=TATGTATG/GTACCCAA"],
                    "ANN": ["GTACCCGC|ann_1|0", "GTACCCGC|ann_2|0", "GTACCCAA|ann_3|1"],
                    "expectedANN": ["GTACCCGC|ann_1|0", "GTACCCGC|ann_2|0", "GTACCCAA|ann_3|1"]
                }),
                # Insertion single nt
                VCFRecord("artificial_chr1", 14, "ins_01", "G", ["GA"], None, None, {
                    "expected": ["artificial_chr1:14=G/GA"],
                    "ANN": ["GA|ann_1|0", "GA|ann_2|0", "GT|ann_3|"],
                    "expectedANN": ["GA|ann_1|0", "GA|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 20, "ins_02", "-", ["A"], None, None, {
                    "expected": ["artificial_chr1:19=T/TA"],
                    "ANN": ["A|ann_1|0", "A|ann_2|0", "T|ann_3|"],
                    "expectedANN": ["TA|ann_1|0", "TA|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 14, "ins_03", "G", ["GA", "GC"], None, None, {
                    "expected": ["artificial_chr1:14=G/GA", "artificial_chr1:14=G/GC"],
                    "ANN": ["GA|ann_1|0", "GA|ann_2|0", "GC|ann_3|1", "GT|ann_4|"],
                    "expectedANN": ["GA|ann_1|0", "GA|ann_2|0", "GC|ann_3|1"]
                }),
                VCFRecord("artificial_chr1", 20, "ins_04", "-", ["A", "C"], None, None, {
                    "expected": ["artificial_chr1:19=T/TA", "artificial_chr1:19=T/TC"],
                    "ANN": ["A|ann_1|0", "A|ann_2|0", "C|ann_3|1", "T|ann_4|"],
                    "expectedANN": ["TA|ann_1|0", "TA|ann_2|0", "TC|ann_3|1"]
                }),
                # Insertion multi nt
                VCFRecord("artificial_chr1", 14, "ins_05", "G", ["GATGC"], None, None, {
                    "expected": ["artificial_chr1:14=G/GATGC"],
                    "ANN": ["GATGC|ann_1|0", "GATGC|ann_2|0", "GAAAC|ann_3|"],
                    "expectedANN": ["GATGC|ann_1|0", "GATGC|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 20, "ins_06", "-", ["AAATC"], None, None, {
                    "expected": ["artificial_chr1:19=T/TAAATC"],
                    "ANN": ["AAATC|ann_1|0", "AAATC|ann_2|0", "GAAAC|ann_3|"],
                    "expectedANN": ["TAAATC|ann_1|0", "TAAATC|ann_2|0"]
                }),
                # Movable insertion multi nt
                VCFRecord("artificial_chr1", 14, "ins_07", "G", ["GTG"], None, None, {
                    "expected": ["artificial_chr1:12=A/ATG"],
                    "ANN": ["GTG|ann_1|0", "GTG|ann_2|0", "GAAAC|ann_3|"],
                    "expectedANN": ["ATG|ann_1|0", "ATG|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 27, "ins_08", "A", ["AAAA"], None, None, {
                    "expected": ["artificial_chr1:25=C/CAAA"],
                    "ANN": ["AAAA|ann_1|0", "AAAA|ann_2|0", "CAAA|ann_3|"],
                    "expectedANN": ["CAAA|ann_1|0", "CAAA|ann_2|0"]
                }),
                # Deletion single nt
                VCFRecord("artificial_chr1", 14, "del_01", "G", [""], None, None, {
                    "expected": ["artificial_chr1:13=TG/T"],
                    "ANN": ["-|ann_1|0", "-|ann_2|0", "T|ann_3|"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 14, "del_02", "G", ["-"], None, None, {
                    "expected": ["artificial_chr1:13=TG/T"],
                    "ANN": ["-|ann_1|0", "-|ann_2|0", "T|ann_3|"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 13, "del_03", "TG", ["T"], None, None, {
                    "expected": ["artificial_chr1:13=TG/T"],
                    "ANN": ["T|ann_1|0", "T|ann_2|0", "A|ann_3|"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0"]
                }),
                VCFRecord("artificial_chr1", 13, "del_04", "TG", ["T", "-"], None, None, {
                    "expected": ["artificial_chr1:13=TG/T", "artificial_chr1:12=ATG/A"],
                    "ANN": ["T|ann_1|0", "T|ann_2|0", "-|ann_3|1"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0", "A|ann_3|1"]
                }),
                # Movable deletion multi nt
                VCFRecord("artificial_chr1", 11, "del_05", "TATG", ["T", "TA", "-"], None, None, {
                    "expected": ["artificial_chr1:11=TATG/T", "artificial_chr1:12=ATG/A", "artificial_chr1:7=CATGT/C"],
                    "ANN": ["T|ann_1|0", "T|ann_2|0", "TA|ann_3|1", "-|ann_4|2"],
                    "expectedANN": ["T|ann_1|0", "T|ann_2|0", "A|ann_3|1", "C|ann_4|2"]
                }),
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sequences, self.tmp_faidx, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testVCFIO(self):
        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = {}
        for record in self.variants:
            for idx, alt in enumerate(record.alt):
                id = "{} {}:{}={}/{}".format(record.id, record.chrom, record.pos, record.ref, alt)
                expected[id] = record.info["expected"][idx]
        observed = {}
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed[record.id + " " + record.info["UNSTD"]] = record.getName()
        self.assertEqual(
            expected,
            observed
        )

    def testAnnotVCFIO(self):
        # Execute command
        subprocess.check_call(self.cmd + ["--annotations-field", "ANN"], stderr=subprocess.DEVNULL)

        # Validate results
        expected = {}
        for record in self.variants:
            for idx, alt in enumerate(record.alt):
                id = "{} {}:{}={}/{}".format(record.id, record.chrom, record.pos, record.ref, alt)
                expected[id] = record.info["expected"][idx]
        observed = {}
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed[record.id + " " + record.info["UNSTD"]] = record.getName()
        self.assertEqual(
            expected,
            observed
        )

        # Validate annotations
        expected = {}
        for record in self.variants:
            for idx, alt in enumerate(record.alt):
                id = "{} {}:{}={}/{}".format(record.id, record.chrom, record.pos, record.ref, alt)
                expected[id] = sorted([ann for ann in record.info["expectedANN"] if ann.split("|")[2] == str(idx)])
        observed = {}
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                id = record.id + " " + record.info["UNSTD"]
                observed[id] = []
                if "ANN" in record.info:
                    observed[id] = sorted([ann for ann in record.info["ANN"]])
        self.assertEqual(
            expected,
            observed
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
