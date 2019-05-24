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
__version__ = '1.1.0'
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

from anacore.vcf import VCFRecord
from anacore.annotVcf import AnnotVCFIO

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFByAnnot(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_selected_rna = os.path.join(tmp_folder, unique_id + "_rna.tsv")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Create RNA ref
        with open(self.tmp_selected_rna, "w") as FH_rna:
            FH_rna.write("#Gene\tTranscript\n")
            FH_rna.write("Gene_1\tENST_selected1\n")
            FH_rna.write("Gene_1\tENST_selected2\n")

        # Create VCF
        with AnnotVCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.ANN_titles = ["Allele", "Consequence", "Feature", "EUR_AF", "gnomAD_AF", "expected_filter"]
            FH_var.info = {
                "ANN": {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|gnomAD_AF|expected_filter."},
                "expected_filter": {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "The expected filters."}
            }
            FH_var._writeHeader()
            self.variants = [
                VCFRecord(
                    "artificial_chr1", 14, "alt_00", "G", ["T"], None, None,
                    {
                        "ANN": [{"Allele": "T", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "PASS"}],
                        "expected_filter": ["PASS"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_01", "G", ["T"], None, None,
                    {
                        "expected_filter": ["CSQ"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_02", "G", ["T"], None, None,
                    {
                        "ANN": [{"Allele": "T", "Consequence": "synonymous_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.CSQ"}],
                        "expected_filter": ["CSQ"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_03", "G", ["T"], None, None,
                    {
                        "ANN": [{"Allele": "T", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.popAF"}],
                        "expected_filter": ["popAF"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_04", "G", ["T"], None, None,
                    {
                        "ANN": [{"Allele": "T", "Consequence": "missense_variant", "Feature": "other", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.RNA"}],
                        "expected_filter": ["CSQ"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_05", "G", ["T"], None, None,
                    {
                        "ANN": [{"Allele": "G", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}],
                        "expected_filter": ["CSQ"]
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_06", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "PASS"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["PASS"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_07", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.popAF"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["popAF"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_08", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "synonymous_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.CSQ"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["CSQ"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_09", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "missense_variant", "Feature": "other", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.RNA"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["CSQ"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_10", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "synonymous_variant", "Feature": "other", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.CSQ&ANN.RNA&ANN.popAF"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["CSQ", "popAF"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_11", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "synonymous_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.CSQ&ANN.popAF"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["CSQ", "popAF"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_12", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "synonymous_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.CSQ&ANN.popAF"},
                            {"Allele": "T", "Consequence": "missense_variant", "Feature": "other", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.01, "expected_filter": "ANN.RNA&ANN.popAF"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.001&0.001", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC"}
                        ],
                        "expected_filter": ["CSQ", "popAF"],
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 14, "alt_13", "G", ["T"], None, None,
                    {
                        "ANN": [
                            {"Allele": "T", "Consequence": "synonymous_variant", "Feature": "ENST_selected1", "EUR_AF": "0.01&0.01", "gnomAD_AF": 0.001, "expected_filter": "ANN.CSQ&ANN.popAF"},
                            {"Allele": "C", "Consequence": "missense_variant", "Feature": "ENST_selected1", "EUR_AF": "0.05&0.05", "gnomAD_AF": 0.001, "expected_filter": "ANN.COLLOC&ANN.popAF"}
                        ],
                        "expected_filter": ["CSQ", "popAF"],
                    }
                )
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_selected_rna, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testResultsRecordRemove(self):
        cmd = [
            "filterVCFByAnnot.py",
            "--mode", "remove",
            "--input-selected-RNA", self.tmp_selected_rna,
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            if record.info["expected_filter"] == ["PASS"]:
                expected.append(record.id)
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                if record.info["expected_filter"] == ["PASS"]:
                    observed.append(record.id)
        self.assertEqual(expected, observed)

    def testResultsRecordTag(self):
        cmd = [
            "filterVCFByAnnot.py",
            "--mode", "tag",
            "--input-selected-RNA", self.tmp_selected_rna,
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            for curr_filter in sorted(record.info["expected_filter"]):
                expected.append("{}:{}".format(record.id, curr_filter))
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                for curr_filter in sorted(record.filter):
                    observed.append("{}:{}".format(record.id, curr_filter))
        self.assertEqual(expected, observed)

    def testResultsAnnotTag(self):
        cmd = [
            "filterVCFByAnnot.py",
            "--mode", "tag",
            "--input-selected-RNA", self.tmp_selected_rna,
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            if "ANN" not in record.info:
                expected.append("{}:None".format(record.id))
            else:
                for annot_idx, annot in enumerate(record.info["ANN"]):
                    for curr_filter in sorted(annot["expected_filter"].split("&")):
                        expected.append("{}:{}:{}".format(record.id, annot_idx, curr_filter))
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                if "ANN" not in record.info or len(record.info["ANN"]) == 0:
                    observed.append("{}:None".format(record.id))
                else:
                    for annot_idx, annot in enumerate(record.info["ANN"]):
                        for curr_filter in sorted(annot["FILTER"].split("&")):
                            observed.append("{}:{}:{}".format(record.id, annot_idx, curr_filter))
        self.assertEqual(expected, observed)

    def testResultsAnnotRemove(self):
        cmd = [
            "filterVCFByAnnot.py",
            "--mode", "remove",
            "--input-selected-RNA", self.tmp_selected_rna,
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            if record.info["expected_filter"] == ["PASS"]:
                annot_idx = 0
                for annot in record.info["ANN"]:
                    if annot["expected_filter"] == "PASS":
                        expected.append("{}:{}:PASS".format(record.id, annot_idx))
                    annot_idx += 1
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                for annot_idx, annot in enumerate(record.info["ANN"]):
                    for curr_filter in sorted(annot["FILTER"].split("&")):
                        observed.append("{}:{}:{}".format(record.id, annot_idx, curr_filter))
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
