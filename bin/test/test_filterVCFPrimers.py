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
            chr1 = Sequence( "artificial_chr1", "NNNAAAATTTGGGGGGGGGGTTTAAANNN" )
            FH_seq.write( chr1 )
            chr2 = Sequence( "artificial_chr2", "NNNAAAATTTGGGGGGGGGGTTTAAANNN" )
            FH_seq.write( chr2 )

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {"ZOI": {"type": str, "type_tag": "String", "number": 1, "number_tag": 1, "description": "If the variant can be in interest area."} }
            FH_var._writeHeader()
            FH_var.write(
                VCFRecord( "artificial_chr1", 6, "arti1", "A", ["AA"], None, None, {"ZOI": "no"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr1", 8, "arti2", "T", ["TT"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr1", 14, "arti3", "G", ["GG"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr1", 22, "arti4", "T", ["TT"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr1", 24, "arti5", "A", ["AA"], None, None, {"ZOI": "no"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr2", 6, "arti6", "A", ["AA"], None, None, {"ZOI": "no"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr2", 8, "arti7", "T", ["TT"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr2", 14, "arti8", "G", ["GG"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr2", 22, "arti9", "T", ["TT"], None, None, {"ZOI": "yes"} )
            )
            FH_var.write(
                VCFRecord( "artificial_chr2", 24, "arti10", "A", ["AA"], None, None, {"ZOI": "no"} )
            )

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sequences, self.tmp_regions, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testOverlapException(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord( "artificial_chr1", 5, 25, "ampl1", None, "+", 10, 20 )
            FH_reg.write( ampl1 )
            ampl2 = BEDRecord( "artificial_chr2", 5, 25, "ampl2", None, "+", 10, 20 )
            FH_reg.write( ampl2 )
            ampl3 = BEDRecord( "artificial_chr1", 23, 28, "ampl2", None, "+", 25, 26 )
            FH_reg.write( ampl3 )

        # Execute command
        with self.assertRaises(subprocess.CalledProcessError) as context:
            subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

    def testResults(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord( "artificial_chr1", 5, 25, "ampl1", None, "+", 10, 20 )
            FH_reg.write( ampl1 )
            ampl2 = BEDRecord( "artificial_chr2", 5, 25, "ampl2", None, "+", 10, 20 )
            FH_reg.write( ampl2 )

        # Execute command
        subprocess.check_call(self.cmd)

        # Validate results
        expected = ["arti2", "arti3", "arti4", "arti7", "arti8", "arti9"]
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
