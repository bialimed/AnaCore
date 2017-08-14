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
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(BIN_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from bed import BEDIO, BEDRecord
from vcf import VCFIO, VCFRecord
from sequenceIO import FastaIO, Sequence


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
        self.tmp_sequences = os.path.join( tmp_folder, unique_id + ".fasta")
        self.tmp_regions = os.path.join( tmp_folder, unique_id + ".bed")
        self.tmp_variants = os.path.join( tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join( tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "filterVCFPrimers.py",
            "--input-variants", self.tmp_variants,
            "--input-regions", self.tmp_regions,
            "--input-sequences", self.tmp_sequences,
            "--output-variants", self.tmp_output
        ]

        # Create fasta
        with FastaIO(self.tmp_sequences, "w") as FH_seq:
            FH_seq.write( Sequence("artificial_chr1", "NNNAAAATTTGGGGGGGGGGTTTAAANNN") )
            #                                          123456789| | | | | | | | | |
            #                                                   10| 14| 18| 22| 26|
            #                                                     12  16  20  24  28
            FH_seq.write( Sequence("artificial_chr2", "CGATNNNCGAT") )
            #                                          123456789|
            #                                                   10

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {"ZOI": {"type": str, "type_tag": "String", "number": 1, "number_tag": 1, "description": "If the variant can be in interest area."} }
            FH_var._writeHeader()
            self.variants = [
                VCFRecord( "artificial_chr1", 6, None, "A", ["AA"], None, None, {"ZOI": "no"} ),
                VCFRecord( "artificial_chr1", 8, None, "TT", ["T"], None, None, {"ZOI": "no"} ),
                VCFRecord( "artificial_chr1", 8, None, "T", ["TT"], None, None, {"ZOI": "yes"} ),
                VCFRecord( "artificial_chr1", 9, None, "TTGG", ["TT"], None, None, {"ZOI": "yes"} ),
                VCFRecord( "artificial_chr1", 14, None, "G", ["GG"], None, None, {"ZOI": "yes"} ),
                VCFRecord( "artificial_chr1", 18, None, "GGG", ["G"], None, None, {"ZOI": "yes"} ), # ZOI downstream limit deletion
                VCFRecord( "artificial_chr1", 22, None, "T", ["TT"], None, None, {"ZOI": "yes"} ),

				VCFRecord( "artificial_chr1", 9, None, "TT", ["TC"], None, None, {"ZOI": "no"} ), # Substitution before end of upstream primer
                VCFRecord( "artificial_chr1", 10, None, "TG", ["TC"], None, None, {"ZOI": "yes"} ), # Substitution in upstream limit of ZOI
				VCFRecord( "artificial_chr1", 15, None, "GG", ["GC"], None, None, {"ZOI": "yes"} ), # Substitution in dosnstream limit of ZOI
				VCFRecord( "artificial_chr1", 20, None, "GT", ["GC"], None, None, {"ZOI": "no"} ), # Substitution after start of downstream primer
				VCFRecord( "artificial_chr1", 21, None, "TT", ["TC"], None, None, {"ZOI": "no"} ), # Substitution in downstream primer

                VCFRecord( "artificial_chr2", 1, None, "C", ["CTT"], None, None, {"ZOI": "no"} ), # Insertion before end of upstream primer
                VCFRecord( "artificial_chr2", 2, None, "G", ["GCC"], None, None, {"ZOI": "yes"} ), # Insertion in upstream limit of ZOI
                VCFRecord( "artificial_chr2", 3, None, "AT", ["CCGC"], None, None, {"ZOI": "yes"} ), # Insertion in upstream limit of ZOI and without standardization
                VCFRecord( "artificial_chr2", 9, None, "G", ["GCC"], None, None, {"ZOI": "yes"} ), # Insertion in downstream limit of ZOI
                VCFRecord( "artificial_chr2", 9, None, "G", ["NNN"], None, None, {"ZOI": "yes"} ), # Insertion in downstream limit of ZOI and without standardization
                VCFRecord( "artificial_chr2", 10, None, ".", ["CC"], None, None, {"ZOI": "yes"} ), # Insertion in downstream limit of ZOI
                VCFRecord( "artificial_chr2", 10, None, "A", ["ATT"], None, None, {"ZOI": "no"} ), # Insertion after start of downstream primer

                VCFRecord( "artificial_chr2", 1, None, "CG", ["C"], None, None, {"ZOI": "no"} ), # Deletion before end of upstream primer
                VCFRecord( "artificial_chr2", 2, None, "GA", ["G"], None, None, {"ZOI": "yes"} ), # Deletion in upstream limit of ZOI
                VCFRecord( "artificial_chr2", 3, None, "AT", ["C"], None, None, {"ZOI": "yes"} ), # Deletion in upstream limit of ZOI and without standardization
                VCFRecord( "artificial_chr2", 6, None, "NNCG", ["N"], None, None, {"ZOI": "yes"} ), # Deletion in downstream limit of ZOI
                VCFRecord( "artificial_chr2", 8, None, "CG", ["C"], None, None, {"ZOI": "yes"} ), # Deletion in downstream limit of ZOI
                VCFRecord( "artificial_chr2", 8, None, "CG", ["T"], None, None, {"ZOI": "yes"} ), # Deletion in downstream limit of ZOI and without standardization
                VCFRecord( "artificial_chr2", 9, None, "GA", ["G"], None, None, {"ZOI": "no"} ), # Insertion after start of downstream primer
                VCFRecord( "artificial_chr2", 10, None, "A", ["."], None, None, {"ZOI": "no"} ), # Insertion after start of downstream primer
                VCFRecord( "artificial_chr2", 10, None, "AT", ["A"], None, None, {"ZOI": "no"} ), # Insertion after start of downstream primer
            ]
            for idx, curr_var in enumerate(self.variants):
                curr_var.id = "alt_" + str(idx)
                FH_var.write( curr_var )

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sequences, self.tmp_regions, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testOverlapException(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord( "artificial_chr1", 5, 25, "ampl1", None, "+", 11, 20 )
            FH_reg.write( ampl1 )
            ampl2 = BEDRecord( "artificial_chr2", 1, 11, "ampl2", None, "+", 3, 9 )
            FH_reg.write( ampl2 )
            ampl3 = BEDRecord( "artificial_chr1", 23, 28, "ampl2", None, "+", 25, 26 )
            FH_reg.write( ampl3 )

        # Execute command
        with self.assertRaises(subprocess.CalledProcessError) as context:
            subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

    def testResults(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord( "artificial_chr1", 5, 25, "ampl1", None, "+", 11, 20 )
            FH_reg.write( ampl1 )
            ampl2 = BEDRecord( "artificial_chr2", 1, 11, "ampl2", None, "+", 3, 9 )
            FH_reg.write( ampl2 )

        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = [curr_var.id for curr_var in self.variants if curr_var.info["ZOI"] == "yes"]
        observed = list()
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
        self.assertEqual(
            "_".join(sorted(expected)),
            "_".join(sorted(observed))
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
