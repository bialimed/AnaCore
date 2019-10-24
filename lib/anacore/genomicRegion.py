# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing genomic regions defines on references sequences (Gene, Transcripts, Proteins, Regulation sequences, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.region import RegionTree, RegionList


class Gene(RegionTree):
    pass


class Exon(RegionTree):
    pass


class Intron(RegionTree):
    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None, parent=None, children=None):
        """
        Build and return an instance of Intron.

        :param start: The start position on the reference. This position is 1-based and ascending (start <= end).
        :type start: int
        :param end: The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        :type end: int
        :param strand: The strand of the instance ("+" or "-").
        :type strand: str
        :param reference: The region object or the region name of the reference.
        :type reference: region.Region | str
        :param name: The name of the region.
        :type name: str
        :param annot: The annotations of the region.
        :type annot: dict
        :param parent: The parent region (example: gene for a transcript).
        :type parent: region.RegionTree
        :param children: The list of sub-regions (example: exons in transcript).
        :type children: region.RegionList
        :return: The new instance.
        :rtype: region.Intron
        """
        RegionTree.__init__(self, start, end, strand, reference, name, annot, None, children)
        self.parent = parent  # Prevent add of this intron as child of transcript


class CDS(RegionTree):
    pass


class Transcript(RegionTree):
    """Class to manage a transcript region."""

    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None, parent=None, children=None, proteins=None):
        """
        Build and return an instance of Transcript.

        :param start: The start position on the reference. This position is 1-based and ascending (start <= end).
        :type start: int
        :param end: The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        :type end: int
        :param strand: The strand of the instance ("+" or "-").
        :type strand: str
        :param reference: The region object or the region name of the reference.
        :type reference: region.Region | str
        :param name: The name of the region.
        :type name: str
        :param annot: The annotations of the region.
        :type annot: dict
        :param parent: The gene.
        :type parent: region.RegionTree
        :param children: The list of exons.
        :type children: region.RegionList
        :param proteins: The list of proteins produced on transcript. It can exists several proteins on transcript: operon and readthrough.
        :type proteins: region.RegionList
        :return: The new instance.
        :rtype: genomicRegion.Transcript
        """
        RegionTree.__init__(self, start, end, strand, reference, name, annot, parent, children)
        self.proteins = RegionList()
        if proteins is not None:
            self.setProteins(proteins)

    def setProteins(self, proteins):
        """
        Change proteins linked with transcript. Before the add of the new proteins the old are unlinked (prot.transcript is set to None).

        :param proteins: The new linked proteins.
        :type proteins: list
        """
        # Remove previous not in new proteins
        for idx_prot, prot in enumerate(self.proteins):
            prot.transcript = None
        self.proteins = RegionList()
        # Add new proteins
        for curr_prot in proteins:
            curr_prot.setTranscript(self)

    def addProtein(self, protein):
        """
        Add link with a protein.

        :param protein: The protein to add.
        :type protein: genomicRegion.Protein
        """
        protein.setTranscript(self)

    def length(self):
        """
        Return length of the transcript (sum of exons lengths).

        :return: The length of the transcript.
        :rtype: int
        """
        length = 0
        for curr_exon in self.children:
            length += curr_exon.length()
        return length

    def getPosOnRef(self, transcript_pos):
        """
        Return coordinate on reference sequence from the coordinate on transcript.

        :param transcript_pos: Coordinate on transcript (1-based).
        :type transcript_pos: int.
        :return: The coordinate on reference sequence (1-based).
        :rtype: int
        """
        exon, region_pos = self.getSubFromRegionPos(transcript_pos)
        return exon.getPosOnRef(region_pos)

    def getPosOnRegion(self, chr_pos):
        """
        Return coordinate on transcript from the coordinate on reference sequence.

        :param chr_pos: The coordinate on reference sequence (1-based).
        :type chr_pos: int
        :return: Coordinate on region (1-based).
        :rtype: int
        """
        raise NotImplementedError

    def getSubFromRegionPos(self, transcript_pos):
        """
        Return exon and coordinate on his from the coordinate on transcript.

        Coordinate on exon is stranded:

        - pos 2 in region test:10-20(+) correspond to ref position 12
        - pos 2 in region test:10-20(-) correspond to ref position 18

        :param transcript_pos: The coordinate on transcript (1-based).
        :type transcript_pos: int
        :return: Exon and coordinate on exon (1-based).
        :rtype: (genomicRegion.Exon, int)
        """
        if transcript_pos > self.length():
            raise ValueError("The position {} is out of transcript {}.".format(transcript_pos, self))
        walk_pos = 0
        region_idx = 0
        child = None
        while child is None:
            curr_child_length = self.children[region_idx].length()
            if walk_pos <= transcript_pos and (walk_pos + curr_child_length) >= transcript_pos:
                child = self.children[region_idx]
            else:
                walk_pos += curr_child_length
                region_idx += 1
        pos_on_child = transcript_pos - walk_pos
        return child, pos_on_child

    def getSubFromRefPos(self, chr_pos):
        """
        Return exon on intron where the chr_pos is located.

        :param chr_pos: The coordinate on reference sequence (1-based).
        :type chr_pos: int
        :return: The index of the sub-region containing the chr_pos (1-based) and the object of this sub-region.
        :rtype: int, (Intron | Exon)
        """
        if self.strand is None:
            raise Exception("Cannot return a region position from the reference position because the strand is None ({}).".format(
                self
            ))
        if chr_pos < self.start or chr_pos > self.end:
            raise ValueError("The position {} is out of transcript {}.".format(chr_pos, self))
        sub_region_idx = None
        sub_region = None
        exons = self.children
        exon_idx = 0
        if self.strand == "+":
            while chr_pos > exons[exon_idx].end:
                exon_idx += 1
            if chr_pos < exons[exon_idx].start:
                intron_start = exons[exon_idx - 1].end + 1
                intron_end = exons[exon_idx].start - 1
                sub_region_idx = exon_idx
                sub_region = Intron(intron_start, intron_end, self.strand, self.reference, "intron_{}".format(sub_region_idx), {"siblings_idx": sub_region_idx}, parent=self)
            else:
                sub_region_idx = exon_idx + 1
                sub_region = exons[exon_idx]
        else:
            while chr_pos < exons[exon_idx].start:
                exon_idx += 1
            if chr_pos > exons[exon_idx].end:
                intron_end = exons[exon_idx - 1].start - 1
                intron_start = exons[exon_idx].end + 1
                sub_region_idx = exon_idx
                sub_region = Intron(intron_start, intron_end, self.strand, self.reference, "intron_{}".format(sub_region_idx), {"siblings_idx": sub_region_idx}, parent=self)
            else:
                sub_region_idx = exon_idx + 1
                sub_region = exons[exon_idx]
        return sub_region, sub_region_idx


class Protein(RegionTree):
    """Class to manage a protein region."""

    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None, parent=None, children=None, transcript=None):
        """
        Build and return an instance of Protein.

        :param start: The start position on the reference. This position is 1-based and ascending (start <= end).
        :type start: int
        :param end: The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        :type end: int
        :param strand: The strand of the instance ("+" or "-").
        :type strand: str
        :param reference: The region object or the region name of the reference.
        :type reference: region.Region | str
        :param name: The name of the region.
        :type name: str
        :param annot: The annotations of the region.
        :type annot: dict
        :param parent: The gene.
        :type parent: region.RegionTree
        :param children: The list of CDS.
        :type children: region.RegionList
        :param transcript: The transcript corresponding to the protein.
        :type transcript: region.Transcript
        :return: The new instance.
        :rtype: genomicRegion.Protein
        """
        RegionTree.__init__(self, start, end, strand, reference, name, annot, parent, children)
        self.setTranscript(transcript)

    def setTranscript(self, transcript):
        """
        Add link between self and his transcript.

        :param transcript: The transcript to link.
        :type transcript: genomicRegion.Transcript
        """
        self.transcript = transcript
        if transcript is not None:
            if self not in transcript.proteins:
                transcript.proteins.append(self)
        else:
            prot_idx = None
            try:
                prot_idx = transcript.proteins.index(self)
            except Exception:  # Not found
                pass
            if prot_idx is not None:
                del transcript.proteins[prot_idx]

    def length(self):
        """
        Return length of the protein in nucleotids (sum of CDS lengths).

        :return: The length of the protein.
        :rtype: int
        """
        length = 0
        for curr_cds in self.children:
            length += curr_cds.length()
        return length

    def aaLength(self):
        """
        Return length of the protein in amino acids (sum of CDS lengths).

        :return: The length of the protein.
        :rtype: int
        """
        return int(self.length() / 3)

    def getPosOnRef(self, protein_pos, codon_pos):
        """
        Return coordinate on reference sequence from the coordinate on protein.

        :param protein_pos: Coordinate of the amino acid on protein (1-based).
        :type protein_pos: int.
        :param codon_pos: Coordinate on the codon (1-based).
        :type codon_pos: int.
        :return: The coordinate on reference sequence (1-based).
        :rtype: int
        """
        cds, cds_pos = self.getSubFromRegionPos(protein_pos, codon_pos)
        return cds.getPosOnRef(cds_pos)

    def getPosOnRegion(self, chr_pos):
        """
        Return coordinate on protein from the coordinate on reference sequence.

        :param chr_pos: The coordinate on reference sequence (1-based).
        :type chr_pos: int
        :return: Coordinate of the amino acid on protein (1-based) and the coordinate on the codon (1-based). Return None if coordinates are in UTRs or in introns.
        :rtype: (int, int) | None
        """
        if chr_pos < self.start or chr_pos > self.end:
            raise ValueError("The position {} is out of protein {}.".format(chr_pos, self))
        chained_cds_pos = self.getNtPosFromRefPos(chr_pos)
        if chained_cds_pos is None:
            return [None, None]
        return(
            int((chained_cds_pos - 1) / 3) + 1,  # aa position (1-based)
            ((chained_cds_pos - 1) % 3) + 1  # position in aa (1-based)
        )

    def getNtPosFromRegionPos(self, protein_pos, codon_pos):
        """
        Return nucleotids coordinate on protein from the amino acids coordinate.

        :param protein_pos: Coordinate of the amino acid on protein (1-based).
        :type protein_pos: int.
        :param codon_pos: Coordinate on the codon (1-based).
        :type codon_pos: int.
        :return: The nucleotids coordinate on protein (1-based).
        :rtype: int
        """
        if protein_pos > self.aaLength():
            raise ValueError("The position {} is out of protein {}.".format(protein_pos, self))
        return ((protein_pos - 1) * 3) + codon_pos

    def getNtPosFromRefPos(self, chr_pos):
        """
        Return nucleotids coordinate on protein from coordinate on reference sequence.

        :param chr_pos: The coordinate on reference sequence (1-based).
        :type chr_pos: int
        :return: The nucleotids coordinate on protein (1-based). Return None if coordinates are in UTRs or in introns.
        :rtype: int | None
        """
        if self.strand is None:
            raise Exception("Cannot return a region position from the reference position because the strand is None ({}).".format(
                self
            ))
        chained_cds_pos = None
        cds = self.children
        cds_idx = 0
        if cds[cds_idx].strand == "+":
            if chr_pos <= cds[-1].end and chr_pos >= cds[0].start:
                chained_cds_pos = 0
                while chr_pos > cds[cds_idx].end:
                    chained_cds_pos += cds[cds_idx].length()
                    cds_idx += 1
                chained_cds_pos += chr_pos - cds[cds_idx].start + 1
        else:
            if chr_pos <= cds[0].end and chr_pos >= cds[-1].start:
                chained_cds_pos = 0
                while chr_pos < cds[cds_idx].start:
                    chained_cds_pos += cds[cds_idx].length()
                    cds_idx += 1
                chained_cds_pos += cds[cds_idx].end - chr_pos + 1
        return chained_cds_pos

    def getSubFromRegionPos(self, protein_pos, codon_pos):
        """
        Return CDS and coordinate on his from the coordinate on protein.
        The coordinate on CDS is stranded:
            With region test:10-20(+) the pos 2 correspond to ref position 12
            With region test:10-20(-) the pos 2 correspond to ref position 18

        :param protein_pos: Coordinate of the amino acid on protein (1-based).
        :type protein_pos: int.
        :param codon_pos: Coordinate on the codon (1-based).
        :type codon_pos: int.
        :return: CDS and coordinate on CDS (1-based).
        :rtype: genomicRegion.CDS, int
        """
        chained_cds_pos = self.getNtPosFromRegionPos(protein_pos, codon_pos)
        walk_pos = 0
        region_idx = 0
        child = None
        while child is None:
            curr_child_length = self.children[region_idx].length()
            if walk_pos <= chained_cds_pos and (walk_pos + curr_child_length) >= chained_cds_pos:
                child = self.children[region_idx]
            else:
                walk_pos += curr_child_length
                region_idx += 1
        pos_on_child = chained_cds_pos - walk_pos
        return child, pos_on_child

    def getCDSFromTranscript(self):
        """
        Return CDS of the protein from the transcript and his exons. This function is used when CDS are not defined in the protein but exons and protein start and end are defined.

        :return: The list of CDS of the protein in protein strand order.
        :rtype: region.Regionlist
        """
        # Check information completion
        if self.transcript is None:
            raise Exception("A link with a transcript is required to return CDS for {}.".format(self))
        if self.start is None or self.end is None:
            raise Exception("Start and end for {} are required to return CDS from transcript {}.".format(self, self.transcript))
        # Exons to CDS
        exons = sorted(self.transcript.children, key=lambda x: (x.start, x.end))
        nb_exons = len(exons)
        idx_exon = 0
        curr_exon = exons[idx_exon]
        while self.start > curr_exon.end:
            idx_exon += 1
            curr_exon = exons[idx_exon]
        cds = RegionList()
        while curr_exon is not None and self.end >= curr_exon.start:
            cds_start = max(self.start, curr_exon.start)
            cds_end = min(self.end, curr_exon.end)
            cds.append(CDS(cds_start, cds_end, self.strand, self.reference))
            idx_exon += 1
            curr_exon = None
            if idx_exon < nb_exons:
                curr_exon = exons[idx_exon]
        # Sort by strand order
        if self.strand == "-":
            cds = RegionList(sorted(cds, key=lambda x: (x.end, x.start), reverse=True))
        # Return
        return cds

    def contains(self, eval_region):
        """
        Return True if at least one CDS has an overlap with eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the region contains the evaluated region.
        :rtype: bool
        """
        cds = self.children if len(self.children) > 0 else self.getCDSFromTranscript()
        return len(cds.getContainers(eval_region)) > 0

    def hasOverlap(self, eval_region):
        """
        Return True if at least one CDS has an overlap with eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the protein has an overlap with evaluated region.
        :rtype: bool
        """
        has_overlap = False
        if self.reference.name == eval_region.reference.name:
            cds = self.children if len(self.children) > 0 else self.getCDSFromTranscript()
            for curr_cds in cds:
                if not curr_cds.start > eval_region.end and not curr_cds.end < eval_region.start:
                    has_overlap = True
        return has_overlap
