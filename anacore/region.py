# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing regions defines on references sequences."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.8.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


class Region:
    """Class to manage region on a sequence."""

    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None):
        """
        Build and return an instance of Region.

        :param start: The start position on the reference. This position is 1-based and ascending (start <= end).
        :type start: int
        :param end: The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        :type end: int
        :param strand: The strand of the instance ("+" or "-").
        :type strand: str
        :param reference: The region object or the region name of the reference.
        :type reference: Region | str
        :param name: The name of the region.
        :type name: str
        :param annot: The annotations of the region.
        :type annot: dict
        :return: The new instance.
        :rtype: region.Region
        """
        if start is not None and end is not None and int(start) > int(end):
            raise Exception("Start must be inferior to end.")
        self.start = None if start is None else int(start)
        self.end = self.start if end is None else int(end)
        self.strand = strand
        self.setReference(reference)
        self.name = name
        self.annot = dict() if annot is None else annot

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        clean_str = self.name
        if clean_str is None:
            if self.reference is not None:
                clean_str = self.reference.name
            if self.start is not None:
                clean_str = "{}:{}-{}".format(clean_str, self.start, self.end)
        return clean_str

    def getCoordinatesStr(self):
        """
        Return printable string representation of the coordinates.

        :return: The printable representation of the coordinates.
        :rtype: str
        """
        return "{}:{}-{}[{}]".format(
            (None if self.reference is None else self.reference.name),
            self.start,
            self.end,
            self.strand
        )

    def length(self):
        """
        Return length of the region.

        :return: The length of the region.
        :rtype: int
        """
        return(self.end - self.start + 1)

    def setReference(self, reference):
        """
        Change the reference region of the instance.

        :param reference: The region object or the region name of the reference.
        :type regreferenceion_pos: region.Region | str
        """
        if reference is None:
            self.reference = None
        elif isinstance(reference, str):
            self.reference = Region(name=reference)
        elif isinstance(reference, Region):
            self.reference = reference
        else:
            raise Exception("The type " + str(type(reference)) + " is not authorized for reference attribute in Region class.")

    def getPosOnRef(self, region_pos):
        """
        Return coordinate on reference sequence from the coordinate on self.

        :param region_pos: Coordinate on region (1-based).
        :type region_pos: int.
        :return: The coordinate on reference sequence (1-based).
        :rtype: int
        """
        if self.strand is None:
            raise Exception("Cannot return a reference position from the region position because the strand is None in {}.".format(
                self
            ))
        reference_pos = None
        if self.strand == '+':
            reference_pos = region_pos + self.start - 1
        else:
            reference_pos = self.end - region_pos + 1
        return reference_pos

    def getPosOnRegion(self, ref_pos):
        """
        Return coordinate on region from the coordinate on reference sequence.

        :param ref_pos: The coordinate on reference sequence (1-based).
        :type ref_pos: int
        :return: Coordinate on region (1-based).
        :rtype: int
        """
        if self.strand is None:
            raise Exception("Cannot return a region position from the reference position because the strand is None in {}.".format(
                self.name
            ))
        if ref_pos < self.start or ref_pos > self.end:
            raise ValueError("The region {} does not contains the position {}.".format(
                self, ref_pos
            ))
        region_pos = None
        if self.strand == '+':
            region_pos = ref_pos - self.start + 1
        else:
            region_pos = self.end - ref_pos + 1
        return region_pos

    def contains(self, eval_region):
        """
        Return True if the region contains the eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the region contains the evaluated region.
        :rtype: bool
        """
        contains = False
        if self.reference.name == eval_region.reference.name:
            if self.start <= eval_region.start and self.end >= eval_region.end:
                contains = True
        return contains

    def strandedContains(self, eval_region):
        """
        Return True if the region contains the eval_region and their are on the same strand.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the region contains the evaluated region.
        :rtype: bool
        """
        contains = False
        if self.strand == eval_region.strand:
            contains = self.contains(eval_region)
        return contains

    def hasOverlap(self, eval_region):
        """
        Return True if the region has an overlap with eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the region has an overlap with evaluated region.
        :rtype: bool
        """
        has_overlap = False
        if self.reference.name == eval_region.reference.name:
            if not self.start > eval_region.end and not self.end < eval_region.start:
                has_overlap = True
        return has_overlap

    def hasStrandedOverlap(self, eval_region):
        """
        Return True if the region has an overlap with eval_region and their are on the same strand.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: True if the region has an overlap with evaluated region.
        :rtype: bool
        """
        has_overlap = False
        if self.strand == eval_region.strand:
            has_overlap = self.hasOverlap(eval_region)
        return has_overlap

    def getMinDist(self, eval_region):
        """
        Return the distance with the eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: The distance between the instance and the evaluated region.
        :rtype: int
        """
        if self.reference.name != eval_region.reference.name:
            raise Exception('The minimal distance between regions cannot be processed because their are located on diffrents reference ("' + self.reference.name + '" vs "' + eval_region.reference.name + '").')
        min_dist = None
        if self.hasOverlap(eval_region):  # Eval region overlaps region
            min_dist = 0
        elif eval_region.end < self.start:  # Eval region is located before region
            min_dist = self.start - eval_region.end
        else:  # Eval region is located after region
            min_dist = eval_region.start - self.end
        return min_dist

    @staticmethod
    def fromStr(interval_str):
        """
        Return region instance from interval string (example: 5:15746211-5746311).

        :param interval_str: interval string in format <ref_name>:<start_pos>-<end_pos>. Positions are 1-based.
        :type interval_str: str
        :return: Region instance from interval string.
        :rtype: anacore.region.Region
        """
        reference, start_end = interval_str.split(":")
        start, end = start_end.split("-")
        return Region(start, end, None, reference)


class RegionTree(Region):
    """Class to manage region with hierarchical relations (example: gene -> transcript -> exons)."""

    def __init__(self, start=None, end=None, strand=None, reference=None, name=None, annot=None, parent=None, children=None):
        """
        Build and return an instance of RegionTree.

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
        :rtype: region.RegionTree
        """
        self.children = RegionList()
        Region.__init__(self, start, end, strand, reference, name, annot)
        if children is not None:
            for child in children:
                self.addChild(child)
        self.parent = parent
        if parent is not None:
            parent.addChild(self)

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        clean_str = Region.__str__(self)
        clean_str = "{} [{} children]".format(clean_str, len(self.children))
        if "feature" in self.annot and self.annot["feature"] is not None:
            clean_str = "{} {}".format(self.annot["feature"], clean_str)
        return clean_str

    def __getattribute__(self, name):
        """
        Return the value of the selected attribute.

        :param name: The name of the selected attribute.
        :type name: str
        :return: The value of the attribute.
        :rtype: *
        """
        value = object.__getattribute__(self, name)
        if value is None:
            if name == "start":
                if len(self.children) > 0:
                    value = min([child.start for child in self.children])
            elif name == "end":
                if len(self.children) > 0:
                    value = max([child.end for child in self.children])
            elif name == "strand":
                if len(self.children) > 0:
                    value = self.children[0].strand
            elif name == "reference":
                if len(self.children) > 0:
                    value = self.children[0].reference
        return value

    def addChild(self, child):
        """
        Add child region in region.

        :param child: The added region.
        :type child: region.RegionTree
        """
        # Check compatibility
        if self.strand is not None and child.strand is not None:
            if self.strand != child.strand:
                raise ValueError("The sub-region {} cannot be added to {} because their strands are different.".format(child, self))
        if self.reference is not None:
            if self.reference.name != child.reference.name:
                raise ValueError("The sub-region {} cannot be added to {} because their reference region are different.".format(child, self))
        # Process
        child.parent = self
        self.children.insert(0, child)
        if child.start is not None:
            self.sortChildren()

    def sortChildren(self):
        """Sort children in order of their apparition on self strand. Add siblings_idx (1-based) in annotations of children."""
        if self.strand is None:
            raise Exception("Cannot sort sub-regions because the strand is None in {}.".format(self))
        # Sort
        if self.strand == "-":
            self.children = RegionList(sorted(self.children, key=lambda x: (x.end, x.start), reverse=True))
        else:
            self.children = RegionList(sorted(self.children, key=lambda x: (x.start, x.end)))
        # Add siblings_idx
        for siblings_idx, child in enumerate(self.children):
            child.annot["siblings_idx"] = siblings_idx + 1  # Siblings idx is 1-based


class RegionList(list):
    """Class to manage regions list. for example all the gene on a chromosome."""

    def __init__(self, regions=None):
        """
        Build and return an instance of RegionList.

        :param regions: List of regions.
        :type regions: list
        :return: The new instance.
        :rtype: region.RegionList
        """
        if regions is not None:
            for curr_region in regions:
                self.append(curr_region)

    def getContainers(self, eval_region):
        """
        Return all the regions that contains the eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: Regions containing the evaluated region.
        :rtype: list
        """
        containers = list()
        for curr_region in self:
            if curr_region.contains(eval_region):
                containers.append(curr_region)
        return containers

    def getOverlapped(self, eval_region):
        """
        Return all the regions that have an overlap with eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: Regions having an overlap with evaluated region.
        :rtype: list
        """
        overlapped = list()
        for curr_region in self:
            if curr_region.hasOverlap(eval_region):
                overlapped.append(curr_region)
        return overlapped

    def getNearests(self, eval_region, select_fct=None):
        """
        Return the nearest region to eval_region.

        :param eval_region: The evaluated region.
        :type eval_region: Region
        :return: The distance between the nearest region and self.
        :rtype: list
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
                        nearests.append(curr_region)
        return min_dist, nearests


def splittedByRef(region_list):
    """
    Return regions list by reference.

    :param region_list: The list of regions.
    :type region_list: RegionList
    :return: Regions list by reference.
    :rtype: dict
    """
    regions_by_ref = {}
    for curr_region in region_list:
        ref_name = curr_region.reference.name
        if ref_name not in regions_by_ref:
            regions_by_ref[ref_name] = RegionList()
        regions_by_ref[ref_name].append(curr_region)
    return regions_by_ref


def mergedRegion(regions, trace=False):
    """
    Return a region corresponding to the merge of all provided regions.

    :param regions: One list of overlapped regions.
    :param type: list
    :param trace: if true the provided regions are store in merged_region.annot["merge_traceback"].
    :type trace: bool
    :return: The region corresponding to the merge of all provided regions.
    :rtype: Region
    """
    # Check reference
    ref = set()
    for curr_region in regions:
        if curr_region.reference is None:
            ref.add(None)
        else:
            ref.add(curr_region.reference.name)
    if len(ref) > 1:
        raise Exception("All the regions in mergedRegion are not defined on the same reference.")
    # Init merged region
    merged_region = Region(
        min([curr.start for curr in regions]),
        max([curr.end for curr in regions]),
        None,
        regions[0].reference,
        None
    )
    # Strand
    strands = list(set([curr.strand for curr in regions]))
    if len(strands) == 1 and strands[0] is not None:
        merged_region.strand = strands[0]
    # Traceback
    if trace:
        merged_region.annot["merge_traceback"] = RegionList()
        for curr_region in regions:
            if "merge_traceback" in curr_region.annot:
                for curr_trace in curr_region.annot["merge_traceback"]:
                    merged_region.annot["merge_traceback"].append(curr_trace)
            else:
                merged_region.annot["merge_traceback"].append(curr_region)
    return merged_region


def consolidated(regions, merge_contiguous=True, trace=False):
    """
    Return consolidated regions: the overlapping regions (and optionally the contiguous) are merged in one region.

    :param regions: List of regions to consolidate.
    :param type: list
    :param merge_contiguous: If true the contiguous regions are also merged.
    :type merge_contiguous: bool
    :param trace: if true the provided regions are store in merged_region.annot["merge_traceback"].
    :type trace: bool
    :return: Consolidated regions.
    :rtype: RegionList
    """
    padding = 1 if merge_contiguous else 0
    merged_regions = RegionList()
    regions_by_ref = splittedByRef(regions)
    for ref_id, ref_regions in sorted(regions_by_ref.items()):
        sorted_regions = sorted(ref_regions, key=lambda x: (x.start, x.end))
        prev_region = sorted_regions[0]
        for curr_region in sorted_regions[1:]:
            if curr_region.start <= prev_region.end + padding:
                prev_region = mergedRegion([prev_region, curr_region], trace)
            else:
                merged_regions.append(prev_region)
                prev_region = curr_region
        merged_regions.append(prev_region)
    return merged_regions


def iterOverlapped(queries, subjects, check_ref=True):
    """
    Return a generator on query and his overlapped subjects. For example you can provide a list of variant (queries) and search the genes overlapped by each variant in the list of genes (subjects).

    :param queries: The queries. They must be defined on the same reference.
    :type queries: list
    :param subjects: The regions where the overlap is searched. They must be defined on the same reference.
    :type subjects: list
    :param check_ref: If True all the queries and subjects are checked to ensure they care defined on the same reference. Use False to gain computation time if you are certain that all elements are defined on the same reference.
    :type check_ref: bool
    :return: Each iteration return one query and the list of subjects overlapped by it. The queries are ordered by ascending coordinates.
    :rtype: A generator on couple (Region, RegionList)
    """
    if check_ref:
        references = set()
        for elt in subjects + queries:
            if elt.reference is None:
                references.add(None)
            else:
                references.add(elt.reference.name)
        if len(references) > 1:
            raise Exception("All the queries and subjects in iterOverlapped are not defined on the same reference.")
    subjects = sorted(subjects, key=lambda x: (x.start, x.end))
    queries = sorted(queries, key=lambda x: (x.start, x.end))
    subject_idx = 0
    nb_subjects = len(subjects)
    for curr_query in queries:
        # Find overlapping transcripts
        while subject_idx < nb_subjects and curr_query.start > subjects[subject_idx].end:
            subject_idx += 1
        first_subject_idx = None
        overlapping_subjects = RegionList()
        while subject_idx < nb_subjects and curr_query.end >= subjects[subject_idx].start:
            curr_subject = subjects[subject_idx]
            if not curr_query.start > curr_subject.end:
                overlapping_subjects.append(curr_subject)
                if first_subject_idx is None:
                    first_subject_idx = subject_idx
            subject_idx += 1
        subject_idx = (subject_idx if first_subject_idx is None else first_subject_idx)  # Return to the first subject overlapping the query
        # Return current query and his overlaps
        yield(curr_query, overlapping_subjects)


def iterOverlappedByRegion(queries_by_chr, subject_by_chr):
    """
    Return a generator on query and his overlapped subjects. For example you can provide a list of variant (queries) and search the genes overlapped by each variant in the list of genes (subjects).

    :param queries_by_chr: By chromosome the queries.
    :type queries_by_chr: dict
    :param subject_by_chr: By chromosome the regions where the overlap is searched.
    :type subject_by_chr: dict
    :return: Each iteration return the chromosome name, one query and the list of subjects overlapped by it. The queries are ordered by ascending coordinates.
    :rtype: A generator on couple (str, Region, RegionList)
    """
    for chrom in sorted(queries_by_chr):
        if chrom not in subject_by_chr:
            for query in queries_by_chr[chrom]:
                yield(chrom, query, RegionList())
        else:
            for query, overlapping_subjects in iterOverlapped(queries_by_chr[chrom], subject_by_chr[chrom], False):
                yield(chrom, query, overlapping_subjects)
