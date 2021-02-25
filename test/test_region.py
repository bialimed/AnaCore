#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.region import Region, RegionTree, RegionList, splittedByRef, iterOverlapped, consolidated


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestRegion(unittest.TestCase):
    def testLength(self):
        self.assertEqual(Region(9, 15, None).length(), 7)
        self.assertEqual(Region(9, 15, "+").length(), 7)
        self.assertEqual(Region(9, 15, "-").length(), 7)
        self.assertEqual(Region(9, None, "-").length(), 1)
        self.assertEqual(Region(9, 9, "-").length(), 1)

    def testGetPosOnRef(self):
        # Forward strand
        f_stranded_region = Region(9, 15, "+")
        self.assertEqual(f_stranded_region.getPosOnRef(1), 9)  # First nt
        self.assertEqual(f_stranded_region.getPosOnRef(3), 11)
        self.assertEqual(f_stranded_region.getPosOnRef(7), 15)  # Last nt
        # self.assertRaises(ValueError, f_stranded_region.getPosOnRef(8))  # Out of region
        # Reverse strand
        r_stranded_region = Region(9, 15, "-")
        self.assertEqual(r_stranded_region.getPosOnRef(1), 15)  # First nt
        self.assertEqual(r_stranded_region.getPosOnRef(3), 13)
        self.assertEqual(r_stranded_region.getPosOnRef(7), 9)  # Last nt
        # self.assertRaises(ValueError, r_stranded_region.getPosOnRef(8))  # Out of region

    def testGetPosOnRegion(self):
        # Forward strand
        f_stranded_region = Region(9, 15, "+", "chr1")
        self.assertEqual(f_stranded_region.getPosOnRegion(9), 1)  # First nt
        self.assertEqual(f_stranded_region.getPosOnRegion(11), 3)
        self.assertEqual(f_stranded_region.getPosOnRegion(15), 7)  # Last nt
        with self.assertRaises(ValueError):
            f_stranded_region.getPosOnRegion(8)  # Out of region
        with self.assertRaises(ValueError):
            f_stranded_region.getPosOnRegion(16)  # Out of region
        # Reverse strand
        r_stranded_region = Region(9, 15, "-", "chr1")
        self.assertEqual(r_stranded_region.getPosOnRegion(15), 1)  # First nt
        self.assertEqual(r_stranded_region.getPosOnRegion(13), 3)
        self.assertEqual(r_stranded_region.getPosOnRegion(9), 7)  # Last nt
        with self.assertRaises(ValueError):
            r_stranded_region.getPosOnRegion(8)  # Out of region
        with self.assertRaises(ValueError):
            r_stranded_region.getPosOnRegion(16)  # Out of region

    def testContains(self):
        container_region = Region(9, 15, "+", "chr1")
        self.assertEqual(
            container_region.contains(Region(9, 9, "+", "chr1")),
            True
        )
        self.assertEqual(
            container_region.contains(Region(15, 15, "+", "chr1")),
            True
        )
        self.assertEqual(
            container_region.contains(Region(12, 13, "+", "chr1")),
            True
        )
        self.assertEqual(
            container_region.contains(Region(9, 15, "+", "chr1")),
            True
        )
        self.assertEqual(
            container_region.contains(Region(8, 14, "+", "chr1")),
            False
        )
        self.assertEqual(
            container_region.contains(Region(10, 16, "+", "chr1")),
            False
        )
        self.assertEqual(
            container_region.contains(Region(8, 16, "+", "chr1")),
            False
        )
        self.assertEqual(
            container_region.contains(Region(8, 8, "+", "chr1")),
            False
        )
        self.assertEqual(
            container_region.contains(Region(16, 16, "+", "chr1")),
            False
        )
        self.assertEqual(
            container_region.contains(Region(12, 13, "+", "chr2")),
            False
        )

    def testStrandedContains(self):
        # Forward strand
        f_region = Region(9, 15, "+", "chr1")
        self.assertEqual(
            f_region.strandedContains(Region(9, 12, "+", "chr1")),
            True
        )
        self.assertEqual(
            f_region.strandedContains(Region(9, 12, "-", "chr1")),
            False
        )
        # Reverse strand
        r_region = Region(9, 15, "-", "chr1")
        self.assertEqual(
            r_region.strandedContains(Region(9, 12, "+", "chr1")),
            False
        )
        self.assertEqual(
            r_region.strandedContains(Region(9, 12, "-", "chr1")),
            True
        )

    def testHasOverlap(self):
        region = Region(9, 15, "+", "chr1")
        self.assertEqual(
            region.hasOverlap(Region(9, 9, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(15, 15, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(12, 13, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(9, 15, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(8, 14, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(10, 16, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(8, 16, "+", "chr1")),
            True
        )
        self.assertEqual(
            region.hasOverlap(Region(8, 8, "+", "chr1")),
            False
        )
        self.assertEqual(
            region.hasOverlap(Region(16, 16, "+", "chr1")),
            False
        )
        self.assertEqual(
            region.hasOverlap(Region(12, 13, "+", "chr2")),
            False
        )

    def testHasStrandedOverlap(self):
        # Forward strand
        f_region = Region(9, 15, "+", "chr1")
        self.assertEqual(
            f_region.hasStrandedOverlap(Region(5, 12, "+", "chr1")),
            True
        )
        self.assertEqual(
            f_region.hasStrandedOverlap(Region(5, 12, "-", "chr1")),
            False
        )
        # Reverse strand
        r_region = Region(9, 15, "-", "chr1")
        self.assertEqual(
            r_region.hasStrandedOverlap(Region(5, 12, "+", "chr1")),
            False
        )
        self.assertEqual(
            r_region.hasStrandedOverlap(Region(5, 12, "-", "chr1")),
            True
        )

    def testGetMinDist(self):
        region = Region(9, 15, "+", "chr1")
        self.assertEqual(
            region.getMinDist(Region(14, 18, "+", "chr1")),
            0
        )
        self.assertEqual(
            region.getMinDist(Region(16, 18, "+", "chr1")),
            1
        )
        self.assertEqual(
            region.getMinDist(Region(1, 5, "+", "chr1")),
            4
        )
        with self.assertRaises(Exception):
            region.getMinDist(Region(1, 5, "+", "chr2"))

    def testIterOverlapped_1(self):
        """Case where it does not exist any overlap between subjects."""
        # Init test data
        sbjct_1 = Region(7, 10, "+", "chr1", "sbjct_1")
        sbjct_2 = Region(14, 17, "+", "chr1", "sbjct_2")
        sbjct_3 = Region(24, 29, "+", "chr1", "sbjct_3")
        subjects = RegionList([sbjct_1, sbjct_2, sbjct_3])
        queries_info = [
            {"query": Region(11, 11, "+", "chr1", "query_l1_01"), "overlapped": []},
            {"query": Region(12, 12, "+", "chr1", "query_l1_02"), "overlapped": []},
            {"query": Region(13, 13, "+", "chr1", "query_l1_03"), "overlapped": []},
            {"query": Region(14, 14, "+", "chr1", "query_l1_04"), "overlapped": [sbjct_2]},
            {"query": Region(15, 15, "+", "chr1", "query_l1_05"), "overlapped": [sbjct_2]},
            {"query": Region(16, 16, "+", "chr1", "query_l1_06"), "overlapped": [sbjct_2]},
            {"query": Region(17, 17, "+", "chr1", "query_l1_07"), "overlapped": [sbjct_2]},
            {"query": Region(18, 18, "+", "chr1", "query_l1_08"), "overlapped": []},
            {"query": Region(19, 19, "+", "chr1", "query_l1_09"), "overlapped": []},
            {"query": Region(20, 20, "+", "chr1", "query_l1_10"), "overlapped": []},
            {"query": Region(21, 21, "+", "chr1", "query_l1_11"), "overlapped": []},
            {"query": Region(22, 22, "+", "chr1", "query_l1_12"), "overlapped": []},
            {"query": Region(11, 13, "+", "chr1", "query_l3_01"), "overlapped": []},
            {"query": Region(12, 14, "+", "chr1", "query_l3_02"), "overlapped": [sbjct_2]},
            {"query": Region(13, 15, "+", "chr1", "query_l3_03"), "overlapped": [sbjct_2]},
            {"query": Region(14, 16, "+", "chr1", "query_l3_04"), "overlapped": [sbjct_2]},
            {"query": Region(15, 17, "+", "chr1", "query_l3_05"), "overlapped": [sbjct_2]},
            {"query": Region(16, 18, "+", "chr1", "query_l3_06"), "overlapped": [sbjct_2]},
            {"query": Region(17, 19, "+", "chr1", "query_l3_07"), "overlapped": [sbjct_2]},
            {"query": Region(18, 20, "+", "chr1", "query_l3_08"), "overlapped": []},
            {"query": Region(19, 21, "+", "chr1", "query_l3_09"), "overlapped": []},
            {"query": Region(20, 22, "+", "chr1", "query_l3_10"), "overlapped": []},
            {"query": Region(21, 23, "+", "chr1", "query_l3_11"), "overlapped": []},
            {"query": Region(13, 17, "+", "chr1", "query_l5_01"), "overlapped": [sbjct_2]},
            {"query": Region(15, 19, "+", "chr1", "query_l5_02"), "overlapped": [sbjct_2]},
            {"query": Region(17, 21, "+", "chr1", "query_l5_03"), "overlapped": [sbjct_2]},
            {"query": Region(13, 19, "+", "chr1", "query_l5_04"), "overlapped": [sbjct_2]},
            {"query": Region(13, 18, "+", "chr1", "query_l6_01"), "overlapped": [sbjct_2]},
            {"query": Region(14, 19, "+", "chr1", "query_l6_02"), "overlapped": [sbjct_2]},
            {"query": Region(15, 20, "+", "chr1", "query_l6_03"), "overlapped": [sbjct_2]},
            {"query": Region(16, 21, "+", "chr1", "query_l6_04"), "overlapped": [sbjct_2]},
            {"query": Region(17, 22, "+", "chr1", "query_l6_05"), "overlapped": [sbjct_2]},
            {"query": Region(18, 23, "+", "chr1", "query_l6_06"), "overlapped": []},
            {"query": Region(13, 19, "+", "chr1", "query_l7_01"), "overlapped": [sbjct_2]},
            {"query": Region(14, 20, "+", "chr1", "query_l7_02"), "overlapped": [sbjct_2]},
            {"query": Region(15, 21, "+", "chr1", "query_l7_03"), "overlapped": [sbjct_2]},
            {"query": Region(16, 22, "+", "chr1", "query_l7_04"), "overlapped": [sbjct_2]},
            {"query": Region(17, 23, "+", "chr1", "query_l7_05"), "overlapped": [sbjct_2]},
            {"query": Region(13, 20, "+", "chr1", "query_l8_01"), "overlapped": [sbjct_2]},
            {"query": Region(14, 21, "+", "chr1", "query_l8_02"), "overlapped": [sbjct_2]},
            {"query": Region(15, 22, "+", "chr1", "query_l8_03"), "overlapped": [sbjct_2]},
            {"query": Region(13, 21, "+", "chr1", "query_l9_01"), "overlapped": [sbjct_2]},
            {"query": Region(14, 22, "+", "chr1", "query_l9_02"), "overlapped": [sbjct_2]},
            {"query": Region(13, 22, "+", "chr1", "query_l10_01"), "overlapped": [sbjct_2]}
        ]
        queries_info = sorted(queries_info, key=lambda x: (x["query"].start, x["query"].end))
        # Independant evaluation
        for curr_eval in queries_info:
            obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped([curr_eval["query"]], subjects)]
            self.assertEqual(obs_overlapped, [curr_eval["overlapped"]])
        # Grouped evaluation
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and overlapped subject starts the list of subject
        shift_subjects = subjects[1:]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, shift_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and overlapped subject ends the list of subject
        pop_subjects = subjects[:-1]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, pop_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)

    def testIterOverlapped_2(self):
        """Case where it exists a partial overlap between two subjects."""
        # Init test datas
        sbjct_1 = Region(7, 10, "+", "chr1", "sbjct_1")
        sbjct_2 = Region(14, 18, "+", "chr1", "sbjct_2")
        sbjct_3 = Region(16, 20, "+", "chr1", "sbjct_3")
        sbjct_4 = Region(24, 29, "+", "chr1", "sbjct_4")
        subjects = RegionList([sbjct_1, sbjct_2, sbjct_3, sbjct_4])
        queries_info = [
            {"query": Region(11, 11, "+", "chr1", "query_l1_01"), "overlapped": []},
            {"query": Region(12, 12, "+", "chr1", "query_l1_02"), "overlapped": []},
            {"query": Region(13, 13, "+", "chr1", "query_l1_03"), "overlapped": []},
            {"query": Region(14, 14, "+", "chr1", "query_l1_04"), "overlapped": [sbjct_2]},
            {"query": Region(15, 15, "+", "chr1", "query_l1_05"), "overlapped": [sbjct_2]},
            {"query": Region(16, 16, "+", "chr1", "query_l1_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 17, "+", "chr1", "query_l1_07"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 18, "+", "chr1", "query_l1_08"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 19, "+", "chr1", "query_l1_09"), "overlapped": [sbjct_3]},
            {"query": Region(20, 20, "+", "chr1", "query_l1_10"), "overlapped": [sbjct_3]},
            {"query": Region(21, 21, "+", "chr1", "query_l1_11"), "overlapped": []},
            {"query": Region(22, 22, "+", "chr1", "query_l1_12"), "overlapped": []},
            {"query": Region(23, 23, "+", "chr1", "query_l1_13"), "overlapped": []},
            {"query": Region(24, 24, "+", "chr1", "query_l1_14"), "overlapped": [sbjct_4]},
            {"query": Region(11, 13, "+", "chr1", "query_l3_01"), "overlapped": []},
            {"query": Region(12, 14, "+", "chr1", "query_l3_02"), "overlapped": [sbjct_2]},
            {"query": Region(13, 15, "+", "chr1", "query_l3_03"), "overlapped": [sbjct_2]},
            {"query": Region(14, 16, "+", "chr1", "query_l3_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 17, "+", "chr1", "query_l3_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 18, "+", "chr1", "query_l3_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 19, "+", "chr1", "query_l3_07"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 20, "+", "chr1", "query_l3_08"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 21, "+", "chr1", "query_l3_09"), "overlapped": [sbjct_3]},
            {"query": Region(20, 22, "+", "chr1", "query_l3_10"), "overlapped": [sbjct_3]},
            {"query": Region(21, 23, "+", "chr1", "query_l3_11"), "overlapped": []},
            {"query": Region(13, 17, "+", "chr1", "query_l5_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 19, "+", "chr1", "query_l5_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 21, "+", "chr1", "query_l5_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 22, "+", "chr1", "query_l5_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 23, "+", "chr1", "query_l5_05"), "overlapped": [sbjct_3]},
            {"query": Region(20, 24, "+", "chr1", "query_l5_06"), "overlapped": [sbjct_3, sbjct_4]},
            {"query": Region(13, 18, "+", "chr1", "query_l6_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 19, "+", "chr1", "query_l6_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 20, "+", "chr1", "query_l6_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 21, "+", "chr1", "query_l6_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 22, "+", "chr1", "query_l6_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 23, "+", "chr1", "query_l6_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 24, "+", "chr1", "query_l6_07"), "overlapped": [sbjct_3, sbjct_4]},
            {"query": Region(13, 19, "+", "chr1", "query_l7_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 20, "+", "chr1", "query_l7_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 21, "+", "chr1", "query_l7_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 22, "+", "chr1", "query_l7_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 23, "+", "chr1", "query_l7_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 24, "+", "chr1", "query_l7_06"), "overlapped": [sbjct_2, sbjct_3, sbjct_4]},
            {"query": Region(19, 25, "+", "chr1", "query_l7_07"), "overlapped": [sbjct_3, sbjct_4]},
            {"query": Region(13, 20, "+", "chr1", "query_l8_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 21, "+", "chr1", "query_l8_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 22, "+", "chr1", "query_l8_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(13, 21, "+", "chr1", "query_l9_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 22, "+", "chr1", "query_l9_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(13, 22, "+", "chr1", "query_l10_01"), "overlapped": [sbjct_2, sbjct_3]}
        ]
        queries_info = sorted(queries_info, key=lambda x: (x["query"].start, x["query"].end))
        # Independant evaluation
        for curr_eval in queries_info:
            obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped([curr_eval["query"]], subjects)]
            self.assertEqual(obs_overlapped, [curr_eval["overlapped"]])
        # Grouped evaluation
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and overlap between subjects starts the list of subjects
        shifted_subjects = subjects[1:]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, shifted_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and overlap between subjects ends the list of subjects
        poped_subjects = subjects[:-1]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = []
        for curr_info in queries_info:
            expec_overlapped.append([elt for elt in curr_info["overlapped"] if elt != sbjct_4])
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, poped_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)

    def testIterOverlapped_3(self):
        """Case where a subject is included in another."""
        # Init test data
        sbjct_1 = Region(7, 10, "+", "chr1", "sbjct_1")
        sbjct_2 = Region(14, 20, "+", "chr1", "sbjct_2")
        sbjct_3 = Region(16, 18, "+", "chr1", "sbjct_3")
        sbjct_4 = Region(24, 29, "+", "chr1", "sbjct_4")
        subjects = RegionList([sbjct_1, sbjct_2, sbjct_3, sbjct_4])
        queries_info = [
            {"query": Region(11, 11, "+", "chr1", "query_l1_01"), "overlapped": []},
            {"query": Region(12, 12, "+", "chr1", "query_l1_02"), "overlapped": []},
            {"query": Region(13, 13, "+", "chr1", "query_l1_03"), "overlapped": []},
            {"query": Region(14, 14, "+", "chr1", "query_l1_04"), "overlapped": [sbjct_2]},
            {"query": Region(15, 15, "+", "chr1", "query_l1_05"), "overlapped": [sbjct_2]},
            {"query": Region(16, 16, "+", "chr1", "query_l1_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 17, "+", "chr1", "query_l1_07"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 18, "+", "chr1", "query_l1_08"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 19, "+", "chr1", "query_l1_09"), "overlapped": [sbjct_2]},
            {"query": Region(20, 20, "+", "chr1", "query_l1_10"), "overlapped": [sbjct_2]},
            {"query": Region(21, 21, "+", "chr1", "query_l1_11"), "overlapped": []},
            {"query": Region(22, 22, "+", "chr1", "query_l1_12"), "overlapped": []},
            {"query": Region(11, 13, "+", "chr1", "query_l3_01"), "overlapped": []},
            {"query": Region(12, 14, "+", "chr1", "query_l3_02"), "overlapped": [sbjct_2]},
            {"query": Region(13, 15, "+", "chr1", "query_l3_03"), "overlapped": [sbjct_2]},
            {"query": Region(14, 16, "+", "chr1", "query_l3_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 17, "+", "chr1", "query_l3_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 18, "+", "chr1", "query_l3_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 19, "+", "chr1", "query_l3_07"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 20, "+", "chr1", "query_l3_08"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 21, "+", "chr1", "query_l3_09"), "overlapped": [sbjct_2]},
            {"query": Region(20, 22, "+", "chr1", "query_l3_10"), "overlapped": [sbjct_2]},
            {"query": Region(21, 23, "+", "chr1", "query_l3_11"), "overlapped": []},
            {"query": Region(13, 17, "+", "chr1", "query_l5_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 19, "+", "chr1", "query_l5_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 21, "+", "chr1", "query_l5_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 22, "+", "chr1", "query_l5_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 23, "+", "chr1", "query_l5_05"), "overlapped": [sbjct_2]},
            {"query": Region(20, 24, "+", "chr1", "query_l5_06"), "overlapped": [sbjct_2, sbjct_4]},
            {"query": Region(13, 18, "+", "chr1", "query_l6_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 19, "+", "chr1", "query_l6_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 20, "+", "chr1", "query_l6_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 21, "+", "chr1", "query_l6_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 22, "+", "chr1", "query_l6_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 23, "+", "chr1", "query_l6_06"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(19, 24, "+", "chr1", "query_l6_07"), "overlapped": [sbjct_2, sbjct_4]},
            {"query": Region(13, 19, "+", "chr1", "query_l7_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 20, "+", "chr1", "query_l7_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 21, "+", "chr1", "query_l7_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(16, 22, "+", "chr1", "query_l7_04"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(17, 23, "+", "chr1", "query_l7_05"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(18, 24, "+", "chr1", "query_l7_06"), "overlapped": [sbjct_2, sbjct_3, sbjct_4]},
            {"query": Region(19, 24, "+", "chr1", "query_l7_07"), "overlapped": [sbjct_2, sbjct_4]},
            {"query": Region(13, 20, "+", "chr1", "query_l8_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 21, "+", "chr1", "query_l8_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(15, 22, "+", "chr1", "query_l8_03"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(13, 21, "+", "chr1", "query_l9_01"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(14, 22, "+", "chr1", "query_l9_02"), "overlapped": [sbjct_2, sbjct_3]},
            {"query": Region(13, 22, "+", "chr1", "query_l10_01"), "overlapped": [sbjct_2, sbjct_3]}
        ]
        queries_info = sorted(queries_info, key=lambda x: (x["query"].start, x["query"].end))
        # Independant evaluation
        for curr_eval in queries_info:
            obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped([curr_eval["query"]], subjects)]
            self.assertEqual(obs_overlapped, [curr_eval["overlapped"]])
        # Grouped evaluation
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and inclusion between subjects starts the list of subjects
        shifted_subjects = subjects[1:]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = [curr_info["overlapped"] for curr_info in queries_info]
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, shifted_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)
        # Grouped evaluation and inclusion between subjects ends the list of subjects
        poped_subjects = subjects[:-1]
        queries = [curr_info["query"] for curr_info in queries_info]
        expec_overlapped = []
        for curr_info in queries_info:
            expec_overlapped.append([elt for elt in curr_info["overlapped"] if elt != sbjct_4])
        obs_overlapped = [overlapped_subjects for query, overlapped_subjects in iterOverlapped(queries, poped_subjects)]
        self.assertEqual(obs_overlapped, expec_overlapped)

    def testFromStr(self):
        observed = Region.fromStr("12:1534187-1534287")
        expected = Region(1534187, 1534287, None, "12")
        self.assertEqual(str(observed), str(expected))


class TestRegionTree(unittest.TestCase):
    def testInit(self):
        exon1 = RegionTree(10, 30, "+", "chr1")
        exon2 = RegionTree(40, 70, "+", "chr1")
        exon3 = RegionTree(80, 100, "+", "chr1")
        RegionTree(children=[exon2, exon3, exon1])
        RegionTree(10, 100, "+", "chr1", children=[exon2, exon3, exon1])

    def testStart(self):
        tr1_region = RegionTree(10, 100, "+", "chr1")
        tr2_region = RegionTree(None, None, None, "chr1")
        exon1_region = RegionTree(10, 30, "+", "chr1")
        exon2_region = RegionTree(40, 70, "+", "chr1")
        exon3_region = RegionTree(80, 100, "+", "chr1")

        tr1_region.addChild(exon1_region)
        tr1_region.addChild(exon2_region)
        tr1_region.addChild(exon3_region)
        self.assertEqual(tr1_region.start, 10)

        tr2_region.addChild(exon1_region)
        tr2_region.addChild(exon2_region)
        tr2_region.addChild(exon3_region)
        self.assertEqual(tr2_region.start, 10)

    def testEnd(self):
        tr1_region = RegionTree(10, 100, "+", "chr1")
        tr2_region = RegionTree(None, None, None, "chr1")
        exon1_region = RegionTree(10, 30, "+", "chr1")
        exon2_region = RegionTree(40, 70, "+", "chr1")
        exon3_region = RegionTree(80, 100, "+", "chr1")

        tr1_region.addChild(exon1_region)
        tr1_region.addChild(exon2_region)
        tr1_region.addChild(exon3_region)
        self.assertEqual(tr1_region.end, 100)

        tr2_region.addChild(exon1_region)
        tr2_region.addChild(exon2_region)
        tr2_region.addChild(exon3_region)
        self.assertEqual(tr2_region.end, 100)

    def testStrand(self):
        # Forward strand
        tr1_region = RegionTree(10, 100, "+", "chr1")
        tr2_region = RegionTree(None, None, None, "chr1")
        exon1_region = RegionTree(10, 30, "+", "chr1")
        exon2_region = RegionTree(40, 70, "+", "chr1")
        exon3_region = RegionTree(80, 100, "+", "chr1")

        tr1_region.addChild(exon1_region)
        tr1_region.addChild(exon2_region)
        tr1_region.addChild(exon3_region)
        self.assertEqual(tr1_region.strand, "+")

        tr2_region.addChild(exon1_region)
        tr2_region.addChild(exon2_region)
        tr2_region.addChild(exon3_region)
        self.assertEqual(tr2_region.strand, "+")

        # Reverse strand
        tr1_region = RegionTree(10, 100, "-", "chr1")
        tr2_region = RegionTree(None, None, None, "chr1")
        exon1_region = RegionTree(10, 30, "-", "chr1")
        exon2_region = RegionTree(40, 70, "-", "chr1")
        exon3_region = RegionTree(80, 100, "-", "chr1")

        tr1_region.addChild(exon1_region)
        tr1_region.addChild(exon2_region)
        tr1_region.addChild(exon3_region)
        self.assertEqual(tr1_region.strand, "-")

        tr2_region.addChild(exon1_region)
        tr2_region.addChild(exon2_region)
        tr2_region.addChild(exon3_region)
        self.assertEqual(tr2_region.strand, "-")

    def testAddChild(self):
        tr1_region = RegionTree(None, None, None, "chr1")
        exon1_region = RegionTree(10, 30, "-", "chr1", "exon1")
        exon2_region = RegionTree(40, 70, "-", "chr1", "exon2")
        exon3_region = RegionTree(80, 100, "-", "chr1", "exon3")
        tr1_region.addChild(exon2_region)
        tr1_region.addChild(exon1_region)
        tr1_region.addChild(exon3_region)
        self.assertEqual(
            [str(child) for child in tr1_region.children],
            [str(exon3_region), str(exon2_region), str(exon1_region)]
        )
        self.assertEqual(
            [str(child.parent) for child in tr1_region.children],
            [str(tr1_region), str(tr1_region), str(tr1_region)]
        )
        exon4_region = RegionTree(105, 120, "+", "chr1")
        with self.assertRaises(ValueError):
            tr1_region.addChild(exon4_region)
        exon5_region = RegionTree(105, 120, "-", "chr2")
        with self.assertRaises(ValueError):
            tr1_region.addChild(exon5_region)


class TestSplittedByRef(unittest.TestCase):
    def testSplittedByRef(self):
        reg_list = RegionList([
            Region(10, 30, "-", "chr1", "region1"),
            Region(40, 70, "-", "chr1", "region2"),
            Region(80, 100, "-", "chr2", "region3")
        ])
        reg_by_chr = splittedByRef(reg_list)
        expected = ["chr1:region1", "chr1:region2", "chr2:region3"]
        observed = []
        for chrom, regions in sorted(reg_by_chr.items()):
            named_regions = []
            for curr_region in regions:
                named_regions.append("{}:{}".format(chrom, curr_region.name))
            observed.extend(named_regions)
        self.assertEqual(expected, observed)


class TestConsolidate(unittest.TestCase):
    def testConsolidate(self):
        reg_list = RegionList([
            Region(5, 9, "-", "chr1", "region1"),
            Region(10, 30, "-", "chr1", "region2"),
            Region(30, 40, "-", "chr1", "region3"),
            Region(35, 39, "-", "chr1", "region4"),
            Region(40, 70, "-", "chr1", "region5"),
            Region(71, 90, "-", "chr1", "region6"),
            Region(92, 100, "-", "chr1", "region7"),
            Region(100, 100, "+", "chr1", "region8"),
            Region(80, 100, "-", "chr2", "region9")
        ])
        # Merge overlapping
        consolidated_reg = consolidated(reg_list, False)
        expected = ["chr1:5-9[-]", "chr1:10-70[-]", "chr1:71-90[-]", "chr1:92-100[None]", "chr2:80-100[-]"]
        observed = [curr.getCoordinatesStr() for curr in consolidated_reg]
        self.assertEqual(expected, observed)
        # Merge overlapping and contiguous
        consolidated_reg = consolidated(reg_list, True)
        expected = ["chr1:5-90[-]", "chr1:92-100[None]", "chr2:80-100[-]"]
        observed = [curr.getCoordinatesStr() for curr in consolidated_reg]
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
