# -*- coding: utf-8 -*-
"""
Classes and functions for reading/writing/processing BED files.

:Code example:

    Test if a file is a BED file

    .. highlight:: python
    .. code-block:: python

        from anacore.bed import BEDIO

        if BEDIO.isValid("test.txt.gz"):
            print("The file is a BED file")

        # Result>
        # The file is a BED file

    Read BED file by line

    .. highlight:: python
    .. code-block:: python

        from anacore.bed import BEDIO

        with BEDIO("test.bed.gz") as reader:
            print("Region", "Name", sep="\\t")
            for record in reader:
                print(
                    "{}:{}-{}({})".format(
                        record.reference.name,
                        record.start,
                        record.end
                        record.strand
                    ),
                    record.name,
                    sep="\\t"
                )

        # Result>
        # Region\tName
        # chr1:158-367(-)\ttarget_01
        # chr1:54899-54999(+)\ttarget_02

    Read whole BED file and store by region

    .. highlight:: python
    .. code-block:: python

        from anacore.bed import getSortedAreasByChr

        areas_by_chr = getSortedAreasByChr("test.bed.gz")
        print("Region", "Name", sep="\\t")
        for chrom, areas_on_curr_chr in areas_by_chr.items():
            for curr_area in areas_on_curr_chr:
                print(
                    "{}:{}-{}({})".format(
                        curr_area.reference.name,
                        curr_area.start,
                        curr_area.end
                        curr_area.strand
                    ),
                    curr_area.name,
                    sep="\\t"
                )

        # Result>
        # Region\tName
        # chr1:158-367(-)\ttarget_01
        # chr1:54899-54999(+)\ttarget_02

    Write BED file from regions

    .. highlight:: python
    .. code-block:: python

        from anacore.bed import BEDIO, BEDRecord
        from anacore.region import Region

        regions = [
            Region(158, 367, "-", "chr1", "target_01"),
            Region(54899, 54999, "+", "chr1", "target_02"),
        ]

        with BEDIO("test.bed.gz", "w", 6) as writer:  # Write BED file with 6 columns
            for curr_region in regions:
                writer.write(
                    BEDRecord.recFromRegion(curr_region)
                )

        # Result>
        # chr1\t157\t367\ttarget_01\t.\t-
        # chr1\t54898\t54999\ttarget_02\t.\t+
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.5.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.region import Region, RegionList
from anacore.abstractFile import AbstractFile


class BEDRecord(Region):
    """Class to manage a BED record."""

    def __init__(self, chrom=None, start=None, end=None, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        """
        Build and return an instance of BEDRecord.

        :param chrom: The name of the chromosome on which the annotation has been defined.
        :type chrom: str
        :param start: The start position on the reference. This position is 1-based and ascending (start <= end).
        :type start: int
        :param end: The end position on the reference. This position is 1-based and ascending (start <= end). [Default: start]
        :type end: int
        :param name: The name of the annotation.
        :type name: str
        :param score: A score between 0 and 1000
        :type score: int
        :param strand: The strand of the annotation ("+" or "-").
        :type strand: str
        :param thickStart: The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). This position is 1-based and ascending (start <= end).
        :type thickStart: int
        :param thickEnd: The ending position at which the feature is drawn thickly (for example the stop codon in gene displays). This position is 1-based and ascending (start <= end).
        :type thickEnd: int
        :param itemRgb: An RGB value of the form R,G,B (e.g. 255,0,0).
        :type itemRgb: list
        :param blockCount: The number of blocks (exons) in the BED line.
        :type blockCount: int
        :param blockSizes: list of the block sizes
        :type blockSizes: list
        :param blockStarts: list of block starts. All of the blockStart positions should be calculated relative to chromStart.
        :type blockStarts: list
        :return: The new instance.
        :rtype: BEDRecord
        """
        Region.__init__(self, start, end, strand, chrom, name)
        self.score = score
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts

    def __getattribute__(self, attr_name):
        attr_value = None
        if attr_name == "chrom":
            attr_value = object.__getattribute__(self, "reference").name
        else:
            attr_value = object.__getattribute__(self, attr_name)
        return attr_value

    def __setattr__(self, attr_name, attr_value):
        if attr_name == "chrom":
            self.setReference(attr_value)
        else:
            object.__setattr__(self, attr_name, attr_value)

    @staticmethod
    def recFromRegion(region_record):
        """
        Return an instance of BEDRecord corresponding to the region_record. The following fields will be None: score, thickStart, thickEnd, itemRgb, blockCount, blockSizes and blockStarts.

        :param region_record:
        :type region_record: anacore.region.Region
        :return: The BEDRecord corresponding to the region_record. The following fields will be None: score, thickStart, thickEnd, itemRgb, blockCount, blockSizes and blockStarts.
        :rtype: BEDRecord
        """
        return BEDRecord(
            region_record.reference,
            region_record.start,
            region_record.end,
            region_record.name,
            None,
            region_record.strand
        )


class BEDIO(AbstractFile):
    """Class to read and write in BED file (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)."""

    def __init__(self, filepath, mode="r", write_nb_col=3):
        """
        Build and return an instance of BEDIO.

        :param filepath: The filepath.
        :type filepath: str.
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str.
        :param write_nb_col: The number of column in BED file when the mode is in ["a", "w"].
        :type write_nb_col: int.
        :return: The new instance.
        :rtype: BEDIO
        """
        AbstractFile.__init__(self, filepath, mode)
        self._write_nb_col = write_nb_col

    def BEDRecordToBEDLine(self, record):
        """
        Return the record in BED format.

        :param record: The record to process.
        :type record: BEDRecord.
        :return: The BED line corresponding to the record.
        :rtype: str
        """
        str_fmt = "\t".join(["{}" for idx in range(self._write_nb_col)])
        line = str_fmt.format(
            record.chrom,
            (record.start - 1),
            record.end,
            ("." if record.name is None else record.name),
            ("." if record.score is None else str(record.score)),
            ("." if record.strand is None else record.strand),
            ("." if record.thickStart is None else str(record.thickStart - 1)),
            ("." if record.thickEnd is None else str(record.thickEnd)),
            ("." if record.itemRgb is None else ",".join(map(str, record.itemRgb))),
            ("." if record.blockCount is None else str(record.blockCount)),
            ("." if record.blockSizes is None else ",".join(map(str, record.blockSizes))),
            ("." if record.blockStarts is None else ",".join(map(str, record.blockStarts)))
        )
        return line

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record (it is not a comment or an header line).

        :param line: The evaluated line.
        :type line: str.
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        is_record = True
        if line.startswith("browser ") or line.startswith("track ") or line.startswith("#"):
            is_record = False
        return is_record

    @staticmethod
    def isValid(filepath):
        """
        Return True if the file can be a BED file.

        :return: True if the file can be a BED file.
        :rtype: bool
        """
        is_valid = False
        try:
            with BEDIO(filepath) as reader:
                record_idx = 0
                for record in reader:
                    if record_idx >= 10:
                        break
                    # Check strand
                    if record.strand is not None and record.strand not in {"+", "-", "."}:
                        raise IOError(
                            "Invalid strand value in {} from file {}.".format(
                                reader.current_line_nb, reader.filepath
                            )
                        )
                    # Check RGB
                    if record.itemRgb is not None and len(record.itemRgb) != 3:
                        raise IOError(
                            "Invalid RGB code in {} from file {}.".format(
                                reader.current_line_nb, reader.filepath
                            )
                        )
                    record_idx += 1
                is_valid = True
        except Exception:
            pass
        return is_valid

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line.
        :rtype: BEDrecord
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        fields[1] = int(fields[1]) + 1  # Start in BED is 0-based
        fields[2] = int(fields[2])
        nb_col = len(fields)
        if nb_col > 3:
            fields[3] = None if fields[3] == "." else fields[3]
            if nb_col > 4:
                fields[4] = None if fields[4] == "." else int(fields[4])  # A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
                if nb_col > 5:
                    fields[5] = None if fields[5] == "." else fields[5]
                    if nb_col > 6:
                        fields[6] = None if fields[6] == "." else int(fields[6]) + 1  # The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
                        if nb_col > 7:
                            fields[7] = None if fields[7] == "." else int(fields[7])  # The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
                            if nb_col > 8:
                                fields[8] = None if fields[8] in {".", "0"} else [int(color) for color in fields[8].split(",")]  # An RGB value of the form R,G,B (e.g. 255,0,0).
                                if nb_col > 9:
                                    fields[9] = None if fields[9] == "." else int(fields[9])  # The number of blocks (exons) in the BED line.
                                    if nb_col > 10:
                                        fields[10] = None if fields[10] == "." else [int(block) for block in fields[10].split(",")]  # A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
                                        if nb_col > 11:
                                            fields[11] = None if fields[11] == "." else [int(block) for block in fields[11].split(",")]  # A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
        return BEDRecord(*fields)

    def write(self, record):
        """
        Write record line in file.

        :param bed: The record to write.
        :type record: anacore.region.Region
        """
        if not issubclass(record.__class__, BEDRecord):
            record = BEDRecord.recFromRegion(record)
        self.file_handle.write(self.BEDRecordToBEDLine(record) + "\n")
        self.current_line_nb += 1

    @staticmethod
    def getMaxNbCol(in_bed):
        """
        Return the maximum number of fields in BED file. The number of fields per line must be consistent throughout any single set of data in an annotation track.

        :param in_bed: The path to the file (format: BED).
        :type in_bed: str
        :return: The maximum number of fields from the first line of the BED file.
        :rtype: int
        """
        nb_col = -1
        with BEDIO(in_bed) as FH_in:
            for record in FH_in:
                nb_col = max(nb_col, FH_in.current_line.count("\t") + 1)
        return None if nb_col == -1 else nb_col


def getAreas(in_bed):
    """
    Return the list of areas from a BED file.

    :param in_bed: The path to the areas description (format: BED).
    :type in_bed: str
    :return: The list of areas.
    :rtype: region.RegionList
    """
    areas = RegionList()
    with BEDIO(in_bed) as FH_panel:
        areas = RegionList(FH_panel.read())
    return areas


def getAreasByChr(in_bed):
    """
    Return by chromosome the list of areas from a BED file.

    :param in_bed: The path to the areas description (format: BED).
    :type in_bed: str
    :return: The list of areas by chromosome (each list is an instance of region.Regionlist).
    :rtype: dict
    """
    areas_by_chr = dict()
    for curr_area in getAreas(in_bed):
        chrom = curr_area.reference.name
        if chrom not in areas_by_chr:
            areas_by_chr[chrom] = RegionList()
        areas_by_chr[chrom].append(curr_area)
    return areas_by_chr


def getSortedAreasByChr(in_bed):
    """
    Return by chromosome the list of sorted areas from a BED file.

    :param in_bed: The path to the areas description (format: BED).
    :type in_bed: str
    :return: The list of sorted areas by chromosome (each list is an instance of region.Regionlist).
    :rtype: dict
    """
    areas_by_chr = {}
    for chrom, areas in getAreasByChr(in_bed).items():
        areas_by_chr[chrom] = RegionList(
            sorted(areas, key=lambda x: (x.start, x.end))
        )
    return areas_by_chr
