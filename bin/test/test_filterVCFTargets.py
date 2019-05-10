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

from anacore.bed import BEDIO, BEDRecord
from anacore.vcf import VCFIO, VCFRecord
from anacore.sequenceIO import FastaIO, Sequence

BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFTargets(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_sequences = os.path.join(tmp_folder, unique_id + ".fasta")
        self.tmp_faidx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        self.tmp_regions = os.path.join(tmp_folder, unique_id + ".bed")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "filterVCFTargets.py",
            "--mode", "remove",
            "--input-variants", self.tmp_variants,
            "--input-targets", self.tmp_regions,
            "--input-reference", self.tmp_sequences,
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

        # Create targets
        with BEDIO(self.tmp_regions, "w", write_nb_col=4) as FH_bed:
            FH_bed.write(BEDRecord("artificial_chr1", 1, 6, "target_1"))
            FH_bed.write(BEDRecord("artificial_chr1", 15, 28, "target_2"))
            FH_bed.write(BEDRecord("artificial_chr1", 38, 45, "target_3"))

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {"target": {"type": str, "type_tag": "String", "number": 1, "number_tag": 1, "description": "The ID of the overlapped target."}}
            FH_var._writeHeader()
            self.variants = [
                # Substit single nt
                VCFRecord("artificial_chr1", 14, "alt_00", "G", ["T"], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 15, "alt_01", "G", ["T"], None, None, {"target": "target_2"}),  # On target ; first nt of target
                VCFRecord("artificial_chr1", 21, "alt_02", "C", ["G"], None, None, {"target": "target_2"}),  # On target
                VCFRecord("artificial_chr1", 28, "alt_03", "A", ["G"], None, None, {"target": "target_2"}),  # On target ; last nt
                VCFRecord("artificial_chr1", 29, "alt_04", "G", ["C"], None, None, {"target": None}),  # After target ; first nt after target
                # Substit multi nt
                VCFRecord("artificial_chr1", 7, "alt_05", "CATGTATG", ["GTACCCGC"], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 11, "alt_06", "TATGTATG", ["GTACCCGC"], None, None, {"target": "target_2"}),  # Overlap target start
                VCFRecord("artificial_chr1", 13, "alt_07", "TGTATGTGCTCACAAAGTA", ["CCCGCCCCTACATTGCAGT"], None, None, {"target": "target_2"}),  # Include target
                VCFRecord("artificial_chr1", 15, "alt_08", "TATGTGCTCACAAA", ["CGCCCCTACATTGC"], None, None, {"target": "target_2"}),  # Exact target
                VCFRecord("artificial_chr1", 21, "alt_09", "CTCACAA", ["GTACCCG"], None, None, {"target": "target_2"}),  # Included by target
                VCFRecord("artificial_chr1", 24, "alt_10", "ACAAAGTA", ["GTACCCG"], None, None, {"target": "target_2"}),  # Overlap target end
                VCFRecord("artificial_chr1", 29, "alt_11", "GTAGTAGAT", ["GTACCCGA"], None, None, {"target": None}),  # After target ; first nt after target
                # Ins single nt
                VCFRecord("artificial_chr1", 14, "alt_12", "G", ["GA"], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 15, "alt_12.2", "-", ["A"], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 15, "alt_13", "A", ["TG"], None, None, {"target": "target_2"}),  # On target ; first nt of target
                VCFRecord("artificial_chr1", 21, "alt_14", "C", ["CG"], None, None, {"target": "target_2"}),  # On target
                VCFRecord("artificial_chr1", 27, "alt_15", "A", ["AT"], None, None, {"target": "target_2"}),  # On target ; last nt
                VCFRecord("artificial_chr1", 28, "alt_15.2", "-", ["T"], None, None, {"target": "target_2"}),  # On target ; last nt
                VCFRecord("artificial_chr1", 28, "alt_16", "A", ["AT"], None, None, {"target": None}),  # After target ; first nt afetr target
                # Movable del multi nt
                VCFRecord("artificial_chr1", 14, "alt_17", "G", ["GT"], None, None, {"target": "target_2"}),  # Movable to first nt of target
                VCFRecord("artificial_chr1", 28, "alt_18", "A", ["AA"], None, None, {"target": "target_2"}),  # Movable to last nt of target
                # Del single nt
                VCFRecord("artificial_chr1", 14, "alt_19", "G", [""], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 15, "alt_20", "T", [""], None, None, {"target": "target_2"}),  # On target ; first nt of target
                VCFRecord("artificial_chr1", 21, "alt_21", "C", [""], None, None, {"target": "target_2"}),  # On target
                VCFRecord("artificial_chr1", 28, "alt_22", "A", [""], None, None, {"target": "target_2"}),  # On target ; last nt
                VCFRecord("artificial_chr1", 29, "alt_23", "G", [""], None, None, {"target": None}),  # After target ; first nt afetr target
                # Del multi nt
                VCFRecord("artificial_chr1", 11, "alt_24", "TATG", ["T"], None, None, {"target": None}),  # Before target ; first nt before target
                VCFRecord("artificial_chr1", 13, "alt_25", "TGTA", ["T"], None, None, {"target": "target_2"}),  # On target ; first nt of target
                VCFRecord("artificial_chr1", 20, "alt_26", "GCTC", ["G"], None, None, {"target": "target_2"}),  # On target
                VCFRecord("artificial_chr1", 27, "alt_27", "AAGT", ["A"], None, None, {"target": "target_2"}),  # On target ; last nt
                VCFRecord("artificial_chr1", 28, "alt_28", "AGT", ["A"], None, None, {"target": None}),  # After target ; first nt afetr target
                # Movable del multi nt
                VCFRecord("artificial_chr1", 7, "alt_29", "CATGT", ["C"], None, None, {"target": "target_2"}),  # On repeat and movable to first nt of target
                VCFRecord("artificial_chr1", 12, "alt_30", "ATG", ["A"], None, None, {"target": "target_2"}),  # Movable to first nt of target
                VCFRecord("artificial_chr1", 28, "alt_31", "AGTA", ["A"], None, None, {"target": "target_2"}),  # Movable to last nt of target
                VCFRecord("artificial_chr1", 30, "alt_32", "TAGT", ["T"], None, None, {"target": "target_2"}),  # On repeat and movable to last nt of target
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sequences, self.tmp_faidx, self.tmp_regions, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = [curr_var.id for curr_var in self.variants if curr_var.info["target"] == "target_2"]
        observed = list()
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
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
