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
import unittest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.dirname(os.path.dirname(CURRENT_DIR)))
sys.path.append(LIB_DIR)

from anacore.region import Region, RegionTree, RegionList, splittedByRef


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


class TestRegionTree(unittest.TestCase):
    def tesStart(self):
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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
