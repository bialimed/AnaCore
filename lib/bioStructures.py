# 
# Copyright (C) 2014 INRA
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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2014 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'


class Reference:
    def __init__(self, id, name, sequence, description=None):
        self.id          = id
        self.name        = name
        self.sequence    = sequence
        self.description = description if (description is not None) else ""

class Region:
    def __init__(self, start, end, strand, reference):
        self.start     = int(start)
        self.end       = int(end)
        self.strand    = strand
        self.reference = reference # object reference

    def length( self ):
        return( self.end - self.start + 1 )


class BiologicalRegion( Region ):
    def __init__(self, start, end, strand, reference, id, name):
        Region.__init__( self, start, end, strand, reference )
        self.id   = id
        self.name = name

class Exon( BiologicalRegion ):
    def __str__(self):
        string = ""
        if self.strand == "+":
            string = self.id + " : " + self.reference.name + "[" + str(self.start) + ":" + str(self.end) + "]"
        else:
            string = self.id + " : " + self.reference.name + "[" + str(self.end) + ":" + str(self.start) + "]"
        return string

    def getPosOnRef( self, region_pos ):
        reference_pos = None
        if self.strand == '+':
            reference_pos = region_pos + self.start - 1
        else:
            reference_pos = self.end - region_pos + 1
        return reference_pos

class CDS( BiologicalRegion ):
    def __str__(self):
        string = ""
        if self.strand == "+":
            string = self.id + " : " + self.reference.name + "[" + str(self.start) + ":" + str(self.end) + "]"
        else:
            string = self.id + " : " + self.reference.name + "[" + str(self.end) + ":" + str(self.start) + "]"
        return string

    def getPosOnRef( self, region_pos ):
        reference_pos = None
        if self.strand == '+':
            reference_pos = region_pos + self.start - 1
        else:
            reference_pos = self.end - region_pos + 1
        return reference_pos

class Transcript:
    def __init__(self, id, name, gene, strand, regions):
        self.id      = id
        self.name    = name
        self.gene    = gene # object gene
        self.strand  = strand
        self.regions = regions # objects Exon (order : start increase if strand + OR start decrease if strand -)

    def __str__(self):
        string = ""
        for current_region in self.regions:
            string += " {" + current_region.__class__.__name__ + " " + current_region.__str__() + "} "
        return self.id + " (" + self.strand + ") : [ " + string + " ]"

    def addRegion( self, exon ):
        if exon.strand == "+":
            if len(self.regions) == 0 or exon.start > self.regions[-1].start:
                self.regions.append( exon )
            else:
                self.regions.insert( 0, exon )
        else: # strand -
            if len(self.regions) == 0 or exon.end < self.regions[0].end:
                self.regions.append( exon )
            else:
                self.regions.insert( 0, exon )

    def length(self):
        length = 0
        for current_region in self.regions:
            length += current_region.length()
        return length

    def getEltByPos( self, transcript_pos ):
        walk_pos = 0
        region_idx = 0
        find = False
        pos_on_region = None
        while not find:
            if walk_pos <= transcript_pos and (walk_pos + self.regions[region_idx].length()) >= transcript_pos:
                find = True
                pos_on_region = transcript_pos - walk_pos
            else:
                walk_pos += self.regions[region_idx].length()
                region_idx += 1
        return self.regions[region_idx], pos_on_region

    def getPosOnRef( self, transcript_pos ):
        reference_pos = None
        region, region_pos = self.getEltByPos( transcript_pos )
        return region.getPosOnRef( region_pos )

class Protein:
    def __init__(self, id, name, transcript, strand, regions, function=None, evidence=None):
        self.id      = id
        self.name    = name
        self.transcript = transcript # object transcript
        self.strand  = strand
        self.regions = regions # objects CDS
        self.function = None
        self.evidence = None

    def __str__(self):
        string = ""
        for current_region in self.regions:
            string += " {" + current_region.__class__.__name__ + " " + current_region.__str__() + "} "
        return self.id + " (" + self.strand + ") : [ " + string + " ]"

    def addRegion( self, cds ):
        if cds.strand == "+":
            if len(self.regions) == 0 or cds.start > self.regions[-1].start:
                self.regions.append( cds )
            else:
                self.regions.insert( 0, cds )
        else: # strand -
            if len(self.regions) == 0 or cds.end < self.regions[0].end:
                self.regions.append( cds )
            else:
                self.regions.insert( 0, cds )

    def length(self, aa_unit=True):
        length = 0
        for current_region in self.regions:
            length += current_region.length()
        if aa_unit:
            return length/3
        return length

    def getEltByPos( self, protein_pos, codon_pos ):
        general_cds_pos = self.getPosOnGlobalCDS( protein_pos, codon_pos )
        walk_pos = 0
        region_idx = 0
        find = False
        pos_on_region = None
        while not find:
            if walk_pos <= general_cds_pos and (walk_pos + self.regions[region_idx].length()) >= general_cds_pos:
                find = True
                pos_on_region = general_cds_pos - walk_pos
            else:
                walk_pos += self.regions[region_idx].length()
                region_idx += 1
        return self.regions[region_idx], pos_on_region

    def getPosOnRef( self, protein_pos, codon_pos ):
        reference_pos = None
        region, region_pos = self.getEltByPos( proetin_pos, codon_pos )
        return region.getPosOnRef( region_pos )

    def getPosOnGlobalCDS( self, protein_pos, codon_pos ):
        return ((protein_pos-1)*3) + codon_pos

    def getFivePrimUTR( self ):
        first_CDS_start = self.regions[0].getPosOnRef(1)
        find = False
        five_prim_UTR_length = 0
        exon_idx = 0
        while not find:
            if first_CDS_start >= self.transcript.regions[exon_idx].start and first_CDS_start <= self.transcript.regions[exon_idx].end:
                find = True
                if self.transcript.regions[exon_idx].strand == "+":
                    five_prim_UTR_length += first_CDS_start - self.transcript.regions[exon_idx].start
                else:
                    five_prim_UTR_length += self.transcript.regions[exon_idx].end - first_CDS_start
            else:
                five_prim_UTR_length += self.transcript.regions[exon_idx].length()
                exon_idx += 1
        return five_prim_UTR_length

