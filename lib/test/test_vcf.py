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
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.dirname(CURRENT_DIR))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFRecord, VCFIO



class TestVCFRecord:
    def getMostUpstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT" # length: 33
        error_msg = "Error in method getMostUpstream() for class VCFRecord."
        # Test fix deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TtTaAGC", ["T"], 230 )
        upstream = deletion.getMostUpstream( ref )
        assert (upstream.pos == 19 and upstream.ref == "TTAAGC" and upstream.alt[0] == "."), error_msg
        # Test homopolymer deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TTT", ["T"], 230 )
        upstream = deletion.getMostUpstream( ref )
        assert (upstream.pos == 17 and upstream.ref == "TT" and upstream.alt[0] == "."), error_msg
        # Test complex deletion
        deletion = VCFRecord( "artificial_1", 25, None, "CGAgCC", ["C"], 230 )
        upstream = deletion.getMostUpstream( ref )
        assert (upstream.pos == 22 and upstream.ref == "AGCCG" and upstream.alt[0] == "."), error_msg
        # Test deletion on start
        deletion = VCFRecord( "artificial_1", 1, None, ref[0], ["."], 230 )
        upstream = deletion.getMostUpstream( ref )
        assert (upstream.pos == 1 and upstream.ref == "N" and upstream.alt[0] == "."), error_msg
        # Test deletion on end
        deletion = VCFRecord( "artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230 )
        upstream = deletion.getMostUpstream( ref )
        assert (upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "."), error_msg
        # Test fix insertion
        insertion = VCFRecord( "artificial_1", 4, None, "N", ["NGGTT"], 230 )
        upstream = insertion.getMostUpstream( ref )
        assert (upstream.pos == 5 and upstream.ref == "." and upstream.alt[0] == "GGTT"), error_msg
        # Test homopolymer insertion
        insertion = VCFRecord( "artificial_1", 18, None, "T", ["TT"], 230 )
        upstream = insertion.getMostUpstream( ref )
        assert (upstream.pos == 17 and upstream.ref == "." and upstream.alt[0] == "T"), error_msg
        # Test complex insertion
        insertion = VCFRecord( "artificial_1", 25, None, "C", ["CGAgCC"], 230 )
        upstream = insertion.getMostUpstream( ref )
        assert (upstream.pos == 22 and upstream.ref == "." and upstream.alt[0] == "AGCCG"), error_msg
        # Test insertion on start
        insertion = VCFRecord( "artificial_1", 1, None, ".", ["A"], 230 )
        upstream = insertion.getMostUpstream( ref )
        assert (upstream.pos == 1 and upstream.ref == "." and upstream.alt[0] == "A"), error_msg
        # Test insertion on end
        insertion = VCFRecord( "artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230 )
        upstream = insertion.getMostUpstream( ref )
        assert (upstream.pos == 34 and upstream.ref == "." and upstream.alt[0] == "G"), error_msg

    def getMostDownstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT" # length: 33
        error_msg = "Error in method getMostDownstream() for class VCFRecord."
        # Test fix deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TtTaAGC", ["T"], 230 )
        downstream = deletion.getMostDownstream( ref )
        assert (downstream.pos == 19 and downstream.ref == "TTAAGC" and downstream.alt[0] == "."), error_msg
        # Test homopolymer deletion
        deletion = VCFRecord( "artificial_1", 18, None, "TTT", ["T"], 230 )
        downstream = deletion.getMostDownstream( ref )
        assert (downstream.pos == 19 and downstream.ref == "TT" and downstream.alt[0] == "."), error_msg
        # Test complex deletion
        deletion = VCFRecord( "artificial_1", 25, None, "CGAgCC", ["C"], 230 )
        downstream = deletion.getMostDownstream( ref )
        assert (downstream.pos == 28 and downstream.ref == "GCCGA" and downstream.alt[0] == "."), error_msg
        # Test deletion on start
        deletion = VCFRecord( "artificial_1", 1, None, ref[0], ["."], 230 )
        downstream = deletion.getMostDownstream( ref )
        assert (downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "."), error_msg
        # Test deletion on end
        deletion = VCFRecord( "artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230 )
        downstream = deletion.getMostDownstream( ref )
        assert (downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "."), error_msg
        # Test fix insertion
        insertion = VCFRecord( "artificial_1", 4, None, "N", ["NGGTT"], 230 )
        downstream = insertion.getMostDownstream( ref )
        assert (downstream.pos == 5 and downstream.ref == "." and downstream.alt[0] == "GGTT"), error_msg
        # Test homopolymer insertion
        insertion = VCFRecord( "artificial_1", 18, None, "T", ["TT"], 230 )
        downstream = insertion.getMostDownstream( ref )
        assert (downstream.pos == 21 and downstream.ref == "." and downstream.alt[0] == "T"), error_msg
        # Test complex insertion
        insertion = VCFRecord( "artificial_1", 25, None, "C", ["CGAgCC"], 230 )
        downstream = insertion.getMostDownstream( ref )
        assert (downstream.pos == 33 and downstream.ref == "." and downstream.alt[0] == "GCCGA"), error_msg
        # Test insertion on start
        insertion = VCFRecord( "artificial_1", 1, None, ".", ["A"], 230 )
        downstream = insertion.getMostDownstream( ref )
        assert (downstream.pos == 1 and downstream.ref == "." and downstream.alt[0] == "A"), error_msg
        # Test insertion on end
        insertion = VCFRecord( "artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230 )
        downstream = insertion.getMostDownstream( ref )
        assert (downstream.pos == 34 and downstream.ref == "." and downstream.alt[0] == "G"), error_msg

    def processTest(self):
        self.getMostUpstream()
        self.getMostDownstream()


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    vcf_record = TestVCFRecord()
    vcf_record.processTest()
