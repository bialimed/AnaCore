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
__version__ = '1.1.0'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.dirname(CURRENT_DIR))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFRecord, VCFIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestVCFRecord(unittest.TestCase):
    def testStandardizeSingleAllele(self):
        # Test substitution one nt
        substitution = VCFRecord( "artificial_1", 18, None, "TC", ["TA"], 230 )
        substitution.standardizeSingleAllele()
        self.assertTrue( substitution.ref == "C" and substitution.alt[0] == "A" and substitution.pos == 19 )
        # Test substitution multi nt
        substitution = VCFRecord( "artificial_1", 18, None, "TCtgA", ["TaGCc"], 230 )
        substitution.standardizeSingleAllele()
        self.assertTrue( substitution.ref == "CTGA" and substitution.alt[0] == "AGCC" and substitution.pos == 19 )
        # Test substitution multi nt with possible remove
        substitution = VCFRecord( "artificial_1", 18, None, "TCtgATT", ["TaGGctT"], 230 )
        substitution.standardizeSingleAllele()
        self.assertTrue( substitution.ref == "CTGA" and substitution.alt[0] == "AGGC" and substitution.pos == 19 )

        # Test insertion
        insertion = VCFRecord( "artificial_1", 18, None, "T", ["TA"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "." and insertion.alt[0] == "A" and insertion.pos == 19 )
        # Test insertion multi nt and remove upstream and downstream
        insertion = VCFRecord( "artificial_1", 18, None, "TGAT", ["TCGAGAT"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "." and insertion.alt[0] == "CGA" and insertion.pos == 19 )
        # Test insertion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord( "artificial_1", 18, None, "TCtgATTAGC", ["TaGGctTATGCGC"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "CTGATTA" and insertion.alt[0] == "AGGCTTATGC" and insertion.pos == 19 )

        # Test deletion
        insertion = VCFRecord( "artificial_1", 18, None, "TA", ["T"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "A" and insertion.alt[0] == "." and insertion.pos == 19 )
        # Test insertion multi nt and remove upstream and downstream
        insertion = VCFRecord( "artificial_1", 18, None, "TCGAGAT", ["TGAT"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "CGA" and insertion.alt[0] == "." and insertion.pos == 19 )
        # Test insertion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord( "artificial_1", 18, None, "TaGGctTATGCGC", ["TCtgATTAGC"], 230 )
        insertion.standardizeSingleAllele()
        self.assertTrue( insertion.ref == "AGGCTTATGC" and insertion.alt[0] == "CTGATTA" and insertion.pos == 19 )

    def testGetMostUpstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT" # length: 33
        # Test fix deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TtTaAGC", ["T"], 230 )
        upstream = deletion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 19 and upstream.ref == "TTAAGC" and upstream.alt[0] == "." )
        # Test homopolymer deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TTT", ["T"], 230 )
        upstream = deletion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 17 and upstream.ref == "TT" and upstream.alt[0] == "." )
        # Test complex deletion
        deletion = VCFRecord( "artificial_1", 25, None, "CGAgCC", ["C"], 230 )
        upstream = deletion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 22 and upstream.ref == "AGCCG" and upstream.alt[0] == "." )
        # Test deletion on start
        deletion = VCFRecord( "artificial_1", 1, None, ref[0], ["."], 230 )
        upstream = deletion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 1 and upstream.ref == "N" and upstream.alt[0] == "." )
        # Test deletion on end
        deletion = VCFRecord( "artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230 )
        upstream = deletion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "." )
        # Test fix insertion
        insertion = VCFRecord( "artificial_1", 4, None, "N", ["NGGTT"], 230 )
        upstream = insertion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 5 and upstream.ref == "." and upstream.alt[0] == "GGTT" )
        # Test homopolymer insertion
        insertion = VCFRecord( "artificial_1", 18, None, "T", ["TT"], 230 )
        upstream = insertion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 17 and upstream.ref == "." and upstream.alt[0] == "T" )
        # Test complex insertion
        insertion = VCFRecord( "artificial_1", 25, None, "C", ["CGAgCC"], 230 )
        upstream = insertion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 22 and upstream.ref == "." and upstream.alt[0] == "AGCCG" )
        # Test insertion on start
        insertion = VCFRecord( "artificial_1", 1, None, ".", ["A"], 230 )
        upstream = insertion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 1 and upstream.ref == "." and upstream.alt[0] == "A" )
        # Test insertion on end
        insertion = VCFRecord( "artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230 )
        upstream = insertion.getMostUpstream( ref )
        self.assertTrue( upstream.pos == 34 and upstream.ref == "." and upstream.alt[0] == "G" )

    def testGetMostDownstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT" # length: 33
        # Test fix deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TtTaAGC", ["T"], 230 )
        downstream = deletion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 19 and downstream.ref == "TTAAGC" and downstream.alt[0] == "." )
        # Test homopolymer deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TTT", ["T"], 230 )
        downstream = deletion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 19 and downstream.ref == "TT" and downstream.alt[0] == "." )
        # Test complex deletion
        deletion = VCFRecord( "artificial_1", 25, None, "CGAgCC", ["C"], 230 )
        downstream = deletion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 28 and downstream.ref == "GCCGA" and downstream.alt[0] == "." )
        # Test deletion on start
        deletion = VCFRecord( "artificial_1", 1, None, ref[0], ["."], 230 )
        downstream = deletion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "." )
        # Test deletion on end
        deletion = VCFRecord( "artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230 )
        downstream = deletion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "." )
        # Test fix insertion
        insertion = VCFRecord( "artificial_1", 4, None, "N", ["NGGTT"], 230 )
        downstream = insertion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 5 and downstream.ref == "." and downstream.alt[0] == "GGTT" )
        # Test homopolymer insertion
        insertion = VCFRecord( "artificial_1", 18, None, "T", ["TT"], 230 )
        downstream = insertion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 21 and downstream.ref == "." and downstream.alt[0] == "T" )
        # Test complex insertion
        insertion = VCFRecord( "artificial_1", 25, None, "C", ["CGAgCC"], 230 )
        downstream = insertion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 33 and downstream.ref == "." and downstream.alt[0] == "GCCGA" )
        # Test insertion on start
        insertion = VCFRecord( "artificial_1", 1, None, ".", ["A"], 230 )
        downstream = insertion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 1 and downstream.ref == "." and downstream.alt[0] == "A" )
        # Test insertion on end
        insertion = VCFRecord( "artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230 )
        downstream = insertion.getMostDownstream( ref )
        self.assertTrue( downstream.pos == 34 and downstream.ref == "." and downstream.alt[0] == "G" )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
