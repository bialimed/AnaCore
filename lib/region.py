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
__status__ = 'dev'


class Region:
    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None):
        self.start = None if start is None else int(start) # The start position with 1-based and ascending positions on reference
        self.end = self.start if end is None else int(end) # The end position with 1-based and ascending positions on reference
        self.strand = strand
        self.setReference(reference)
        self.name = name
        self.annot = annot

    def length( self ):
        return( self.end - self.start + 1 )

    def setReference(self, reference):
        if reference is None:
            self.reference = None
        elif isinstance(reference, str):
            self.reference = Region(name=reference)
        elif isinstance(reference, Region):
            self.reference = reference
        else:
            raise Exception( "The type " + str(type(reference)) + " is not authorized for reference attribute in Region class." )

    def contains(self, eval_region):
        contains = False
        if self.reference.name == eval_region.reference.name:
            if self.start <= eval_region.start and self.end >= eval_region.end:
                contains = True
        return( contains )

    def strandedContains(self, eval_region):
        contains = False
        if self.reference.strand == eval_region.reference.strand:
            contains = self.contains( eval_region )
        return( contains )        

    def hasOverlap(self, eval_region):
        has_overlap = False
        if self.reference.name == eval_region.reference.name:
            if not self.start > eval_region.end and not self.end < eval_region.start:
                has_overlap = True
        return( has_overlap )

    def strandedHasoverlap(self, eval_region):
        has_overlap = False
        if self.reference.strand == eval_region.reference.strand:
            has_overlap = self.hasOverlap( eval_region )
        return( has_overlap )

    def getMinDist(self, eval_region):
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

#~ def getSelectedArea( input_panel ):
    #~ """
    #~ @summary: Returns the list of selected areas from a BED file.
    #~ @param input_panel: [str] The path to the selected areas description (format: BED).
    #~ @return: [list] The list of BED's areas. Each area is represented by a dictionary with this format: {"region":"chr1", "start":501, "end":608, "id":"gene_98"}.
    #~ """
    #~ selected_areas = list()
    #~ with open(input_panel) as FH_panel:
        #~ for line in FH_panel:
            #~ if not line.startswith("browser ") and not line.startswith("track ") and not line.startswith("#"):
                #~ fields = [elt.strip() for elt in line.split("\t")]
                #~ selected_areas.append({
                    #~ "region": fields[0],
                    #~ "start": int(fields[1]) +1, # Start in BED is 0-based
                    #~ "end": int(fields[2]),
                    #~ "id": fields[3]
                #~ })
    #~ return( selected_areas )

#~ def hasNoOverlap( area_A, area_B ):
    #~ """
    #~ @summary: Returns True if the area_A does not overlap area_B.
    #~ @param area_A: The first evaluated area.
    #~ @param area_B: The second evaluated area.
    #~ @return: [bool] True if the area does not overlap the last area in area_thread.
    #~ """


class RegionList(list):
    def __init__(self, regions=None):
        if regions is not None:
            for curr_region in regions:
                self.append( curr_region )

    def getContainer(self, eval_region):
        containers = list()
        for curr_region in self:
            if curr_region.contains( eval_region ):
                containers.append( curr_region )
        return containers

    def getOverlapped(self, eval_region):
        overlapped = list()
        for curr_region in self:
            if curr_region.hasOverlap( eval_region ):
                overlapped.append( curr_region )
        return overlapped

    def getNearests(self, eval_region, select_fct=None):
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
