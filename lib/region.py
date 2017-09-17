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
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'


class Region:
    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None):
        """
        @param start: [int] The start position on the reference. This position is 1-based and ascending (start <= end).
        @param end: [int] The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        @param strand: [str] The strand of the instance ("+" or "-").
        @param reference: [Region|str] The region object or the region name of the reference.
        @param name: [str] The name of the region.
        @param annot: [dict] The annotations of the region.
        """
        self.start = None if start is None else int(start)
        self.end = self.start if end is None else int(end)
        self.strand = strand
        self.setReference(reference)
        self.name = name
        self.annot = dict() if annot is None else annot

    def length( self ):
        """
        @summary: Returns the length of the region.
        @return: [int] The length of the region.
        """
        return( self.end - self.start + 1 )

    def setReference(self, reference):
        """
        @summary: Changes the reference region of the instance.
        @param reference: [Region|str] The region object or the region name of the reference.
        """
        if reference is None:
            self.reference = None
        elif isinstance(reference, str):
            self.reference = Region(name=reference)
        elif isinstance(reference, Region):
            self.reference = reference
        else:
            raise Exception( "The type " + str(type(reference)) + " is not authorized for reference attribute in Region class." )

    def contains(self, eval_region):
        """
        @summary: Returns True if the region contains the eval_region.
        @param eval_region: The evaluated region.
        @return: [bool] True if the region contains the evaluated region.
        """
        contains = False
        if self.reference.name == eval_region.reference.name:
            if self.start <= eval_region.start and self.end >= eval_region.end:
                contains = True
        return( contains )

    def strandedContains(self, eval_region):
        """
        @summary: Returns True if the region contains the eval_region and their are on the same strand.
        @param eval_region: The evaluated region.
        @return: [bool] True if the region contains the evaluated region.
        """
        contains = False
        if self.reference.strand == eval_region.reference.strand:
            contains = self.contains( eval_region )
        return( contains )

    def hasOverlap(self, eval_region):
        """
        @summary: Returns True if the region has an overlap with eval_region.
        @param eval_region: The evaluated region.
        @return: [bool] True if the region has an overlap with evaluated region.
        """
        has_overlap = False
        if self.reference.name == eval_region.reference.name:
            if not self.start > eval_region.end and not self.end < eval_region.start:
                has_overlap = True
        return( has_overlap )

    def strandedHasoverlap(self, eval_region):
        """
        @summary: Returns True if the region has an overlap with eval_region and their are on the same strand.
        @param eval_region: The evaluated region.
        @return: [bool] True if the region has an overlap with evaluated region.
        """
        has_overlap = False
        if self.reference.strand == eval_region.reference.strand:
            has_overlap = self.hasOverlap( eval_region )
        return( has_overlap )

    def getMinDist(self, eval_region):
        """
        @summary: Returns the distance with the eval_region.
        @param eval_region: The evaluated region.
        @return: [int] The distance between the instance and the evaluated region.
        """
        if self.reference.name != eval_region.reference.name:
            raise Exception( 'The minimal distance between regions cannot be processed because their are located on diffrents reference ("' + self.reference.name+ '" vs "' + eval_region.reference.name + '").' )
        min_dist = None
        if self.hasOverlap(eval_region): # Eval region overlaps region
            min_dist = 0
        elif eval_region.end < self.start: # Eval region is located before region
            min_dist = self.start - eval_region.end
        else: # Eval region is located after region
            min_dist = eval_region.start - self.end
        return( min_dist )


class RegionList(list):
    def __init__(self, regions=None):
        """
        @param regions: [list] The list of regions.
        """
        if regions is not None:
            for curr_region in regions:
                self.append( curr_region )

    def getContainers(self, eval_region):
        """
        @summary: Returns all the regions that contains the eval_region.
        @param eval_region: The evaluated region.
        @return: [list] The regions that contains the evaluated region.
        """
        containers = list()
        for curr_region in self:
            if curr_region.contains( eval_region ):
                containers.append( curr_region )
        return containers

    def getOverlapped(self, eval_region):
        """
        @summary: Returns all the regions that have an overlap with eval_region.
        @param eval_region: The evaluated region.
        @return: [list] The regions that have an overlap with evaluated region.
        """
        overlapped = list()
        for curr_region in self:
            if curr_region.hasOverlap( eval_region ):
                overlapped.append( curr_region )
        return overlapped

    def getNearests(self, eval_region, select_fct=None):
        """
        @summary: Returns the nearest region to eval_region.
        @param eval_region: The evaluated region.
        @return: [list] The distance with the nearest region and the nearest region himself.
        """
        nearests = list()
        min_dist = None
        for curr_region in self:
            selected = True
            if select_fct is not None:
                selected = select_fct(curr_region)
            if selected:
                if curr_region.reference.name == eval_region.reference.name:
                    curr_dist = curr_region.getMinDist(eval_region)
                    if min_dist is None or curr_dist <= min_dist:
                        if min_dist is None or curr_dist < min_dist:
                            nearests = list()
                        min_dist = curr_dist
                        nearests.append( curr_region )
        return( min_dist, nearests )
