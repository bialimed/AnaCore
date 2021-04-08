#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.6.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.genomicRegion import Gene, Transcript, Exon, Intron, Protein, CDS
from anacore.region import Region
from anacore.sequenceIO import IdxFastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestTranscript(unittest.TestCase):
    def setUp(self):
        # Forward
        self.fwd = {
            "exon_1": Exon(10, 30, "+", "chr1"),
            "exon_2": Exon(40, 70, "+", "chr1"),
            "exon_3": Exon(80, 100, "+", "chr1"),
            "intron_1": Intron(31, 39, "+", "chr1"),
            "intron_2": Intron(71, 79, "+", "chr1")
        }
        self.fwd["transcript"] = Transcript(
            10, 100, "+", "chr1", children=[
                self.fwd["exon_1"], self.fwd["exon_3"], self.fwd["exon_2"]
            ]
        )
        # Reverse
        self.rvs = {
            "exon_1": Exon(10, 30, "-", "chr1"),
            "exon_2": Exon(40, 70, "-", "chr1"),
            "exon_3": Exon(80, 100, "-", "chr1"),
            "intron_1": Intron(31, 39, "-", "chr1"),
            "intron_2": Intron(71, 79, "-", "chr1")
        }
        self.rvs["transcript"] = Transcript(
            reference="chr1",
            children=[
                self.rvs["exon_1"], self.rvs["exon_3"], self.rvs["exon_2"]
            ]
        )

    def addProtein(self):
        protein_1 = Protein(10, 30, "+", "chr1", "p1")
        protein_2 = Protein(32, 50, "+", "chr1", "p2")
        transcript_1 = Transcript(name="tr1")
        # Empty
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            []
        )
        # Add protein_1
        transcript_1.addProtein(protein_1)
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            ["p1"]
        )  # Check from tr
        self.assertEqual(
            protein_1.transcript,
            transcript_1
        )  # Check from prot
        # Add protein_2
        transcript_1.addProtein(protein_2)
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            ["p1", "p2"]
        )  # Check from tr
        self.assertEqual(
            protein_1.transcript,
            transcript_1
        )  # Check from prot
        self.assertEqual(
            protein_2.transcript,
            transcript_1
        )  # Check from prot

    def delProtein(self):
        protein_1 = Protein(10, 30, "+", "chr1", "p1")
        protein_2 = Protein(32, 50, "+", "chr1", "p2")
        transcript_1 = Transcript(name="tr1", proteins=[protein_1, protein_2])
        # Init
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            ["p1", "p2"]
        )  # Check from tr
        self.assertEqual(
            protein_1.transcript,
            transcript_1
        )  # Check from prot
        self.assertEqual(
            protein_2.transcript,
            transcript_1
        )  # Check from prot
        # Delete protein 2
        transcript_1.delProtein(protein_2)
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            ["p1"]
        )  # Check from tr
        self.assertEqual(
            protein_1.transcript,
            transcript_1
        )  # Check from prot
        self.assertEqual(
            protein_2.transcript,
            None
        )  # Check from prot
        # Delete protein 1
        transcript_1.addProtein(protein_1)
        self.assertEqual(
            [elt.name for elt in transcript_1.proteins],
            []
        )  # Check from tr
        self.assertEqual(
            protein_1.transcript,
            None
        )  # Check from prot
        self.assertEqual(
            protein_2.transcript,
            None
        )  # Check from prot

    def testSetProteins(self):
        # By init
        protein_1 = Protein(10, 30, "+", "chr1", "p1")
        protein_2 = Protein(32, 50, "+", "chr1", "p2")
        transcript_1 = Transcript(name="tr1", proteins=[protein_1, protein_2])
        self.assertEqual(
            [prot.name for prot in transcript_1.proteins],
            [protein_1.name, protein_2.name]
        )
        self.assertEqual(
            [protein_1.transcript.name, protein_2.transcript.name],
            [transcript_1.name, transcript_1.name]
        )
        # By method
        protein_1 = Protein(10, 30, "+", "chr1", "p1")
        protein_2 = Protein(32, 50, "+", "chr1", "p2")
        transcript_1 = Transcript(name="tr1")
        transcript_1.proteins = [protein_1, protein_2]
        self.assertEqual(
            [prot.name for prot in transcript_1.proteins],
            [protein_1.name, protein_2.name]
        )
        self.assertEqual(
            [protein_1.transcript.name, protein_2.transcript.name],
            [transcript_1.name, transcript_1.name]
        )
        # Replace proteins
        protein_1 = Protein(10, 30, "+", "chr1", "p1")
        protein_2 = Protein(32, 50, "+", "chr1", "p2")
        protein_3 = Protein(54, 70, "+", "chr1", "p3")
        transcript_1 = Transcript(name="tr1", proteins=[protein_1, protein_3])
        transcript_1.proteins = [protein_2, protein_3]
        self.assertIsNone(protein_1.transcript)
        self.assertNotIn(protein_1, transcript_1.proteins)
        self.assertEqual(protein_2.transcript, transcript_1)
        self.assertIn(protein_2, transcript_1.proteins)
        self.assertEqual(protein_3.transcript, transcript_1)
        self.assertIn(protein_3, transcript_1.proteins)

    def testLength(self):
        self.assertEqual(self.fwd["transcript"].length(), 73)
        self.assertEqual(self.rvs["transcript"].length(), 73)

    def testGetSubFromRegionPos(self):
        # Forward
        self.assertEqual(
            self.fwd["transcript"].getSubFromRegionPos(54),
            (self.fwd["exon_3"], 2)
        )
        # Reverse
        self.assertEqual(
            self.rvs["transcript"].getSubFromRegionPos(73),
            (self.rvs["exon_1"], 21)
        )

    def testGetPosOnRegion(self):
        # Forward
        self.assertEqual(self.fwd["transcript"].getPosOnRegion(10), 1)
        self.assertEqual(self.fwd["transcript"].getPosOnRegion(20), 11)
        self.assertEqual(self.fwd["transcript"].getPosOnRegion(30), 21)
        self.assertEqual(self.fwd["transcript"].getPosOnRegion(40), 22)
        self.assertEqual(self.fwd["transcript"].getPosOnRegion(80), 53)
        # Reverse
        self.assertEqual(self.rvs["transcript"].getPosOnRegion(100), 1)
        self.assertEqual(self.rvs["transcript"].getPosOnRegion(90), 11)
        self.assertEqual(self.rvs["transcript"].getPosOnRegion(80), 21)
        self.assertEqual(self.rvs["transcript"].getPosOnRegion(70), 22)
        self.assertEqual(self.rvs["transcript"].getPosOnRegion(30), 53)
        # Exception out of transcript
        with self.assertRaises(Exception):
            self.fwd["transcript"].getPosOnRegion(8)
        # Exception in intron
        with self.assertRaises(Exception):
            self.fwd["transcript"].getPosOnRegion(31)
        # Exception strand
        unstranded_tr = Transcript(10, 30, None, "chr1")
        with self.assertRaises(Exception):
            unstranded_tr.getPosOnRegion(20)

    def testGetPosOnRef(self):
        # Forward
        self.assertEqual(self.fwd["transcript"].getPosOnRef(1), 10)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(10), 19)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(21), 30)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(22), 40)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(35), 53)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(52), 70)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(53), 80)
        self.assertEqual(self.fwd["transcript"].getPosOnRef(73), 100)
        # Reverse
        self.assertEqual(self.rvs["transcript"].getPosOnRef(1), 100)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(10), 91)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(21), 80)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(22), 70)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(35), 57)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(52), 40)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(53), 30)
        self.assertEqual(self.rvs["transcript"].getPosOnRef(73), 10)

    def testGetSubFromRefPos(self):
        # Forward
        fw_tests = [
            [10, (self.fwd["exon_1"], 1)],
            [20, (self.fwd["exon_1"], 1)],
            [30, (self.fwd["exon_1"], 1)],
            [31, (self.fwd["intron_1"], 1)],
            [35, (self.fwd["intron_1"], 1)],
            [39, (self.fwd["intron_1"], 1)],
            [40, (self.fwd["exon_2"], 2)],
            [50, (self.fwd["exon_2"], 2)],
            [70, (self.fwd["exon_2"], 2)],
            [71, (self.fwd["intron_2"], 2)],
            [75, (self.fwd["intron_2"], 2)],
            [79, (self.fwd["intron_2"], 2)],
            [80, (self.fwd["exon_3"], 3)],
            [90, (self.fwd["exon_3"], 3)],
            [100, (self.fwd["exon_3"], 3)]
        ]
        for param, expected in fw_tests:
            observed = self.fwd["transcript"].getSubFromRefPos(param)
            if isinstance(expected[0], Exon):
                self.assertEqual(observed, expected)
            else:
                self.assertEqual(observed[1], expected[1])
                self.assertEqual(observed[0].getCoordinatesStr(), expected[0].getCoordinatesStr())
            self.assertEqual(observed[0].annot["siblings_idx"], expected[1])  # Test siblings_idx
        # Reverse
        rvs_tests = [
            [10, (self.rvs["exon_1"], 3)],
            [20, (self.rvs["exon_1"], 3)],
            [30, (self.rvs["exon_1"], 3)],
            [31, (self.rvs["intron_1"], 2)],
            [35, (self.rvs["intron_1"], 2)],
            [39, (self.rvs["intron_1"], 2)],
            [40, (self.rvs["exon_2"], 2)],
            [50, (self.rvs["exon_2"], 2)],
            [70, (self.rvs["exon_2"], 2)],
            [71, (self.rvs["intron_2"], 1)],
            [75, (self.rvs["intron_2"], 1)],
            [79, (self.rvs["intron_2"], 1)],
            [80, (self.rvs["exon_3"], 1)],
            [90, (self.rvs["exon_3"], 1)],
            [100, (self.rvs["exon_3"], 1)]
        ]
        for param, expected in rvs_tests:
            observed = self.rvs["transcript"].getSubFromRefPos(param)
            if isinstance(expected[0], Exon):
                self.assertEqual(observed, expected)
            else:
                self.assertEqual(observed[1], expected[1])
                self.assertEqual(observed[0].getCoordinatesStr(), expected[0].getCoordinatesStr())
            self.assertEqual(observed[0].annot["siblings_idx"], expected[1])  # Test siblings_idx


class TestProtein(unittest.TestCase):
    def setUp(self):
        # Forward
        self.fwd = {
            "cds_1": CDS(10, 30, "+", "chr1"),
            "cds_2": CDS(40, 70, "+", "chr1"),
            "cds_3": CDS(80, 99, "+", "chr1")
        }
        self.fwd["protein"] = Protein(
            children=[
                self.fwd["cds_1"], self.fwd["cds_3"], self.fwd["cds_2"]
            ]
        )
        # Reverse
        self.rvs = {
            "cds_1": CDS(10, 30, "-", "chr1"),
            "cds_2": CDS(40, 70, "-", "chr1"),
            "cds_3": CDS(80, 99, "-", "chr1")
        }
        self.rvs["protein"] = Protein(
            children=[
                self.rvs["cds_1"], self.rvs["cds_3"], self.rvs["cds_2"]
            ]
        )

    def testLength(self):
        self.assertEqual(self.fwd["protein"].length(), 72)

    def testAaLength(self):
        self.assertEqual(self.fwd["protein"].aaLength(), 24)

    def testGetPosOnRef(self):
        # Forward
        self.assertEqual(self.fwd["protein"].getPosOnRef(1, 1), 10)
        self.assertEqual(self.fwd["protein"].getPosOnRef(1, 3), 12)
        self.assertEqual(self.fwd["protein"].getPosOnRef(2, 1), 13)
        self.assertEqual(self.fwd["protein"].getPosOnRef(7, 3), 30)
        self.assertEqual(self.fwd["protein"].getPosOnRef(8, 1), 40)
        self.assertEqual(self.fwd["protein"].getPosOnRef(18, 1), 70)
        self.assertEqual(self.fwd["protein"].getPosOnRef(18, 2), 80)
        self.assertEqual(self.fwd["protein"].getPosOnRef(24, 1), 97)
        self.assertEqual(self.fwd["protein"].getPosOnRef(24, 3), 99)
        # Reverse
        self.assertEqual(self.rvs["protein"].getPosOnRef(1, 1), 99)
        self.assertEqual(self.rvs["protein"].getPosOnRef(1, 3), 97)
        self.assertEqual(self.rvs["protein"].getPosOnRef(2, 1), 96)
        self.assertEqual(self.rvs["protein"].getPosOnRef(7, 2), 80)
        self.assertEqual(self.rvs["protein"].getPosOnRef(7, 3), 70)
        self.assertEqual(self.rvs["protein"].getPosOnRef(17, 3), 40)
        self.assertEqual(self.rvs["protein"].getPosOnRef(18, 1), 30)
        self.assertEqual(self.rvs["protein"].getPosOnRef(24, 1), 12)
        self.assertEqual(self.rvs["protein"].getPosOnRef(24, 3), 10)

    def testGetPosOnRegion(self):
        # Forward
        self.assertEqual(self.fwd["protein"].getPosOnRegion(10), (1, 1))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(12), (1, 3))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(13), (2, 1))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(30), (7, 3))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(40), (8, 1))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(70), (18, 1))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(80), (18, 2))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(97), (24, 1))
        self.assertEqual(self.fwd["protein"].getPosOnRegion(99), (24, 3))
        # Reverse
        self.assertEqual(self.rvs["protein"].getPosOnRegion(99), (1, 1))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(97), (1, 3))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(96), (2, 1))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(80), (7, 2))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(70), (7, 3))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(40), (17, 3))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(30), (18, 1))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(12), (24, 1))
        self.assertEqual(self.rvs["protein"].getPosOnRegion(10), (24, 3))

    def testGetNtPosFromRegionPos(self):
        # Forward
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(1, 1), 1)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(1, 3), 3)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(2, 1), 4)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(7, 3), 21)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(8, 1), 22)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(18, 1), 52)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(18, 2), 53)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(24, 1), 70)
        self.assertEqual(self.fwd["protein"].getNtPosFromRegionPos(24, 3), 72)
        # Reverse
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(1, 1), 1)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(1, 3), 3)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(2, 1), 4)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(7, 2), 20)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(7, 3), 21)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(17, 3), 51)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(18, 1), 52)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(24, 1), 70)
        self.assertEqual(self.rvs["protein"].getNtPosFromRegionPos(24, 3), 72)

    def testGetNtPosFromRefPos(self):
        # Forward
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(10), 1)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(12), 3)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(13), 4)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(30), 21)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(40), 22)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(70), 52)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(80), 53)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(97), 70)
        self.assertEqual(self.fwd["protein"].getNtPosFromRefPos(99), 72)
        # Reverse
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(99), 1)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(97), 3)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(96), 4)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(80), 20)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(70), 21)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(40), 51)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(30), 52)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(12), 70)
        self.assertEqual(self.rvs["protein"].getNtPosFromRefPos(10), 72)

    def testGetSubFromRegionPos(self):
        # Forward
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(1, 1),
            (self.fwd["cds_1"], 1)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(1, 3),
            (self.fwd["cds_1"], 3)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(2, 1),
            (self.fwd["cds_1"], 4)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(7, 3),
            (self.fwd["cds_1"], 21)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(8, 1),
            (self.fwd["cds_2"], 1)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(18, 1),
            (self.fwd["cds_2"], 31)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(18, 2),
            (self.fwd["cds_3"], 1)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(24, 1),
            (self.fwd["cds_3"], 18)
        )
        self.assertEqual(
            self.fwd["protein"].getSubFromRegionPos(24, 3),
            (self.fwd["cds_3"], 20)
        )
        # Reverse
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(1, 1),
            (self.rvs["cds_3"], 1)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(1, 3),
            (self.rvs["cds_3"], 3)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(2, 1),
            (self.rvs["cds_3"], 4)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(7, 2),
            (self.rvs["cds_3"], 20)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(7, 3),
            (self.rvs["cds_2"], 1)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(17, 3),
            (self.rvs["cds_2"], 31)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(18, 1),
            (self.rvs["cds_1"], 1)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(24, 1),
            (self.rvs["cds_1"], 19)
        )
        self.assertEqual(
            self.rvs["protein"].getSubFromRegionPos(24, 3),
            (self.rvs["cds_1"], 21)
        )

    def testHasOverlap(self):
        # Forward
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(8, 9, "+", "chr1")), False)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(9, 9, "+", "chr1")), False)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(9, 10, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(10, 10, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(12, 12, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(12, 29, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(30, 30, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(30, 31, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(31, 31, "+", "chr1")), False)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(37, 39, "+", "chr1")), False)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(39, 39, "+", "chr1")), False)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(39, 40, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(40, 40, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(98, 99, "+", "chr1")), True)
        self.assertEqual(self.fwd["protein"].hasOverlap(Region(100, 102, "+", "chr1")), False)
        # Reverse
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(8, 9, "-", "chr1")), False)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(9, 9, "-", "chr1")), False)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(9, 10, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(10, 10, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(12, 12, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(12, 29, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(30, 30, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(30, 31, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(31, 31, "-", "chr1")), False)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(37, 39, "-", "chr1")), False)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(39, 39, "-", "chr1")), False)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(39, 40, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(40, 40, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(98, 99, "-", "chr1")), True)
        self.assertEqual(self.rvs["protein"].hasOverlap(Region(100, 102, "-", "chr1")), False)

    def testGetCDSFromTranscript(self):
        res = []
        # One exon forward
        transcript_1 = Transcript(None, None, "+", "chr1", children=[
            Exon(100, 150, "+", "chr1")
        ])
        res.append({
            "expected": [CDS(100, 150, "+", "chr1")],
            "observed": Protein(100, 150, "+", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 140, "+", "chr1")],
            "observed": Protein(110, 140, "+", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "+", "chr1")],
            "observed": Protein(110, 150, "+", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 140, "+", "chr1")],
            "observed": Protein(100, 140, "+", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        # Three exons forward
        transcript_2 = Transcript(None, None, "+", "chr1", children=[
            Exon(30, 80, "+", "chr1"),
            Exon(100, 150, "+", "chr1"),
            Exon(170, 200, "+", "chr1")
        ])
        res.append({
            "expected": [CDS(100, 150, "+", "chr1")],
            "observed": Protein(100, 150, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 140, "+", "chr1")],
            "observed": Protein(110, 140, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "+", "chr1")],
            "observed": Protein(110, 150, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 140, "+", "chr1")],
            "observed": Protein(100, 140, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(30, 80, "+", "chr1"), CDS(100, 150, "+", "chr1")],
            "observed": Protein(30, 150, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(60, 80, "+", "chr1"), CDS(100, 150, "+", "chr1")],
            "observed": Protein(60, 150, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(80, 80, "+", "chr1"), CDS(100, 150, "+", "chr1")],
            "observed": Protein(80, 150, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(80, 80, "+", "chr1"), CDS(100, 150, "+", "chr1"), CDS(170, 170, "+", "chr1")],
            "observed": Protein(80, 170, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 150, "+", "chr1"), CDS(170, 200, "+", "chr1")],
            "observed": Protein(100, 200, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 150, "+", "chr1"), CDS(170, 190, "+", "chr1")],
            "observed": Protein(100, 190, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "+", "chr1"), CDS(170, 200, "+", "chr1")],
            "observed": Protein(110, 200, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "+", "chr1"), CDS(170, 190, "+", "chr1")],
            "observed": Protein(110, 190, "+", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        # One exon reverse
        transcript_1 = Transcript(None, None, "-", "chr1", children=[
            Exon(100, 150, "-", "chr1")
        ])
        res.append({
            "expected": [CDS(100, 150, "-", "chr1")],
            "observed": Protein(100, 150, "-", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 140, "-", "chr1")],
            "observed": Protein(110, 140, "-", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "-", "chr1")],
            "observed": Protein(110, 150, "-", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 140, "-", "chr1")],
            "observed": Protein(100, 140, "-", "chr1", transcript=transcript_1).getCDSFromTranscript()
        })
        # Three exons reverse
        transcript_2 = Transcript(None, None, "-", "chr1", children=[
            Exon(170, 200, "-", "chr1"),
            Exon(100, 150, "-", "chr1"),
            Exon(30, 80, "-", "chr1")
        ])
        res.append({
            "expected": [CDS(100, 150, "-", "chr1")],
            "observed": Protein(100, 150, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 140, "-", "chr1")],
            "observed": Protein(110, 140, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(110, 150, "-", "chr1")],
            "observed": Protein(110, 150, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 140, "-", "chr1")],
            "observed": Protein(100, 140, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 150, "-", "chr1"), CDS(30, 80, "-", "chr1")],
            "observed": Protein(30, 150, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 150, "-", "chr1"), CDS(60, 80, "-", "chr1")],
            "observed": Protein(60, 150, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(100, 150, "-", "chr1"), CDS(80, 80, "-", "chr1")],
            "observed": Protein(80, 150, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(170, 170, "-", "chr1"), CDS(100, 150, "-", "chr1"), CDS(80, 80, "-", "chr1")],
            "observed": Protein(80, 170, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(170, 200, "-", "chr1"), CDS(100, 150, "-", "chr1")],
            "observed": Protein(100, 200, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(170, 190, "-", "chr1"), CDS(100, 150, "-", "chr1")],
            "observed": Protein(100, 190, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(170, 200, "-", "chr1"), CDS(110, 150, "-", "chr1")],
            "observed": Protein(110, 200, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        res.append({
            "expected": [CDS(170, 190, "-", "chr1"), CDS(110, 150, "-", "chr1")],
            "observed": Protein(110, 190, "-", "chr1", transcript=transcript_2).getCDSFromTranscript()
        })
        # Launch evaluation
        for eval_pair in res:
            self.assertEqual(
                ", ".join([curr_cds.getCoordinatesStr() for curr_cds in eval_pair["expected"]]),
                ", ".join([curr_cds.getCoordinatesStr() for curr_cds in eval_pair["observed"]]),
            )

    def testGetCodonRefPos(self):
        tr_1 = Transcript(None, None, "+", "chr1", children=[
            Exon(10, 20, "+", "chr1"),
            Exon(30, 40, "+", "chr1"),
            Exon(45, 50, "+", "chr1")
        ])
        prot_1 = Protein(16, 36, "+", "chr1", transcript=tr_1)
        tr_2 = Transcript(None, None, "-", "chr1", children=[
            Exon(10, 20, "-", "chr1"),
            Exon(30, 40, "-", "chr1"),
            Exon(45, 50, "-", "chr1")
        ])
        prot_2 = Protein(16, 36, "-", "chr1", transcript=tr_2)
        data = [
            {"prot": prot_1, "aa_pos": 1, "expected": [16, 17, 18]},
            {"prot": prot_1, "aa_pos": 2, "expected": [19, 20, 30]},
            {"prot": prot_1, "aa_pos": 3, "expected": [31, 32, 33]},
            {"prot": prot_1, "aa_pos": 4, "expected": [34, 35, 36]},
            {"prot": prot_2, "aa_pos": 1, "expected": [36, 35, 34]},
            {"prot": prot_2, "aa_pos": 2, "expected": [33, 32, 31]},
            {"prot": prot_2, "aa_pos": 3, "expected": [30, 20, 19]},
            {"prot": prot_2, "aa_pos": 4, "expected": [18, 17, 16]},
        ]
        for curr in data:
            self.assertEqual(
                curr["prot"].getCodonRefPos(curr["aa_pos"]),
                curr["expected"]
            )
            with self.assertRaises(Exception):
                self.prot_1.getCodonRefPos(8)  # Not in protein


class TestProteinSeq(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fasta_idx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        self.tmp_fasta = os.path.join(tmp_folder, unique_id + ".fasta")

        # Create sequence file
        content_fasta = """>one
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
>two another chromosome
ATGCATGCATGCAT
GCATGCATGCATGC"""
        with open(self.tmp_fasta, "w") as FH_out:
            FH_out.write(content_fasta)

        # Proteins
        # 1 3 5 7 9  11 14 16 19 21        30 33 36 39 41  44 47 50
        # ATGCATGCAT GCATG CATGC ATGCATGCA TGCATGCATGC ATGCATGCATGCATGCATGC
        #            ..... *****           *******....      ......
        # prot_1           12345           6789 11
        tr_1 = Transcript(None, None, "+", "one", children=[
            Exon(10, 20, "+", "one"),
            Exon(30, 40, "+", "one"),
            Exon(45, 50, "+", "one")
        ])
        self.prot_1 = Protein(16, 36, "+", "one", transcript=tr_1)
        tr_2 = Transcript(None, None, "-", "one", children=[
            Exon(10, 20, "-", "one"),
            Exon(30, 40, "-", "one"),
            Exon(45, 50, "-", "one")
        ])
        self.prot_2 = Protein(16, 36, "-", "one", transcript=tr_2)

        # Create index
        content_fasta_idx = """one	66	5	30	31
two	28	98	14	15"""
        with open(self.tmp_fasta_idx, "w") as FH_out:
            FH_out.write(content_fasta_idx)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta, self.tmp_fasta_idx]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGetCodonSeqFromProtPos(self):
        data = [
            {"prot": self.prot_1, "aa_pos": 1, "expected": "CAT"},
            {"prot": self.prot_1, "aa_pos": 2, "expected": "GCT"},
            {"prot": self.prot_1, "aa_pos": 4, "expected": "TGC"},
            {"prot": self.prot_2, "aa_pos": 1, "expected": "GCA"},
            {"prot": self.prot_2, "aa_pos": 3, "expected": "AGC"},
            {"prot": self.prot_2, "aa_pos": 4, "expected": "ATG"},
        ]
        with IdxFastaIO(self.tmp_fasta) as reader:
            for curr in data:
                self.assertEqual(
                    curr["prot"].getCodonSeqFromProtPos(curr["aa_pos"], reader),
                    curr["expected"]
                )
            with self.assertRaises(Exception):
                self.prot_1.getCodonSeqFromProtPos(8, reader)  # Not in protein

    def testGetCodonInfo(self):
        data = [
            {"prot": self.prot_1, "ref_pos": 16, "expected": (1, 1, "CAT")},
            {"prot": self.prot_1, "ref_pos": 18, "expected": (1, 3, "CAT")},
            {"prot": self.prot_1, "ref_pos": 19, "expected": (2, 1, "GCT")},
            {"prot": self.prot_1, "ref_pos": 20, "expected": (2, 2, "GCT")},
            {"prot": self.prot_1, "ref_pos": 34, "expected": (4, 1, "TGC")},
            {"prot": self.prot_1, "ref_pos": 35, "expected": (4, 2, "TGC")},
            {"prot": self.prot_1, "ref_pos": 36, "expected": (4, 3, "TGC")},
            {"prot": self.prot_2, "ref_pos": 36, "expected": (1, 1, "GCA")},
            {"prot": self.prot_2, "ref_pos": 34, "expected": (1, 3, "GCA")},
            {"prot": self.prot_2, "ref_pos": 30, "expected": (3, 1, "AGC")},
            {"prot": self.prot_2, "ref_pos": 20, "expected": (3, 2, "AGC")},
            {"prot": self.prot_2, "ref_pos": 19, "expected": (3, 3, "AGC")},
            {"prot": self.prot_2, "ref_pos": 18, "expected": (4, 1, "ATG")},
            {"prot": self.prot_2, "ref_pos": 16, "expected": (4, 3, "ATG")}
        ]
        with IdxFastaIO(self.tmp_fasta) as reader:
            for curr in data:
                self.assertEqual(
                    curr["prot"].getCodonInfo(curr["ref_pos"], reader),
                    curr["expected"]
                )
            with self.assertRaises(Exception):
                self.prot_1.getCodonInfo(14, reader)  # In exon but not in CDS


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
