# -*- coding: utf-8 -*-
"""
Classes and functions for reading/writing/processing wiggle file.

:Example:

    Test if a file is a wiggle file

    .. highlight:: python
    .. code-block:: python

        from anacore.wig import WigIO

        if WigIO.isValid("test.txt.gz"):
            print("The file is a wiggle file")

        # Result>
        # The file is a wiggle file

    Read wiggle file by line

    .. highlight:: python
    .. code-block:: python

        from anacore.wig import WigIO

        with WigIO("test.wig.gz") as reader:
            print("region", "mean_depth", sep="\\t")
            for record in reader:
                print(
                    "{}:{}-{}".format(
                        record.reference.name,
                        record.start,
                        record.end
                    ),
                    record.annot["value"],
                    sep="\\t"
                )

        # Result>
        # region\tmean_depth
        # chr1:15800-15899\t12
        # chr1:15900-15999\t40
        # chr1:889100-889399\t22.5
        # chr2:200100-200199\t13.4
        # chr2:280100-280100\t5

    Write wiggle

    .. highlight:: python
    .. code-block:: python

        from anacore.wig import WigIO

        with WigIO("test.wig.gz", "w", step=100, span=100) as writer:
            writer.write("chr1", 15800, 12)
            writer.write("chr1", 15900, 40)
            writer.write("chr1", 889100, 22.5)
            writer.step = None
            writer.write("chr2", 200100, 13.4)
            writer.span = None
            writer.write("chr2", 280100, 5)

        # File content>
        # fixeStep chrom=chr1 start=15800 step=100 span=100
        # 12
        # 40
        # fixeStep chrom=chr1 start=889100 step=100 span=100
        # 22.5
        # variableStep chrom=chr2 span=100
        # 200100 13.4
        # variableStep chrom=chr2 span=1
        # 280100 5
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

from anacore.abstractFile import AbstractFile, isEmpty, isGzip
from anacore.region import Region
import gzip


def getInfoFromFixedDef(def_line):
    """
    Return information from definition line as tuple.

    :param def_line: FixedStep definition line from wiggle.
    :type def_line: str
    :return: Information from definition: chrom, pos, step, span.
    :rtype: (str, int, int, int)
    """
    # fixeStep\schrom=.+\sstart=\d+\sstep=\d+(\sspan=\d+)?
    if " span=" in def_line:
        step_type, chrom_part, start_part, step_part, span_part = def_line.split(" ")
        span = int(span_part.split("=", 1)[1])
    else:
        step_type, chrom_part, start_part, step_part = def_line.split(" ")
        span = 1
    chrom = chrom_part.split("=", 1)[1]
    pos = int(start_part.split("=", 1)[1])
    step = int(step_part.split("=", 1)[1])
    return (chrom, pos, step, span)


def getInfoFromVariableDef(def_line):
    """
    Return information from definition line as tuple.

    :param def_line: VariableStep definition line from wiggle.
    :type def_line: str
    :return: Information from definition: chrom, None, None, span.
    :rtype: (str, None, None, int)
    """
    if " span=" in def_line:
        step_type, chrom_part, span_part = def_line.split(" ")
        span = int(span_part.split("=", 1)[1])
    else:
        step_type, chrom_part = def_line.split(" ")
        chrom = chrom_part.split("=", 1)[1]
        span = 1
    chrom = chrom_part.split("=", 1)[1]
    return chrom, None, None, span


class WigIO(AbstractFile):
    """Class to read and write in wiggle file (https://genome.ucsc.edu/goldenPath/help/wiggle.html)."""

    def __init__(self, filepath, mode="r", step=None, span=None):
        """
        Build and return an instance of WigIO.

        :param filepath: The filepath.
        :type filepath: str.
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str.
        :param step: In write mode: distance between data values. If this value is None data are in varaibleStep otherwise they are in fiwedStep.
        :type step: int.
        :param span: In write mode: number of bases that each data value should cover [Default: 1].
        :type span: int.
        :return: The new instance.
        :rtype: WigIO
        """
        AbstractFile.__init__(self, filepath, mode)
        self.step = step  # In write mode set to None convert writer from fixedStep to variableStep
        self.span = 1 if span is None else span
        if mode != "a":
            self._prev = {"chrom": None, "pos": None, "step": None, "span": None}
        else:  # Append mode
            if not isEmpty(filepath):
                # Get _prev
                open_fct = open
                mode = "r"
                if isGzip(filepath):
                    open_fct = gzip.open
                    mode = "rt"
                with open_fct(filepath, mode) as reader:
                    for line in reader:
                        if line.startswith("fixedStep"):
                            prev_chrom, prev_pos, prev_step, prev_span = getInfoFromFixedDef(line.rstrip())
                            prev_pos = prev_pos - prev_step
                        elif line.startswith("variableStep"):
                            prev_chrom, prev_pos, prev_step, prev_span = getInfoFromVariableDef(line.rstrip())
                        else:
                            if prev_step is not None:
                                prev_pos += prev_step
                    end_with_line_return = line.strip() == ""
                self._prev = {"chrom": prev_chrom, "pos": prev_pos, "step": prev_step, "span": prev_span}
                if step is None:
                    self.step = self._prev["step"]
                if span is None:
                    self.span = self._prev["span"]
                # Add first return
                if not end_with_line_return:
                    self.file_handle.write("\n")

    @staticmethod
    def isValid(filepath):
        """
        Return True if the file can be a WIG file.

        :return: True if the file can be a WIG file.
        :rtype: bool
        """
        # Check is empty
        if isEmpty(filepath):
            return True
        # Check first line
        open_fct = open
        mode = "r"
        if isGzip(filepath):
            open_fct = gzip.open
            mode = "rt"
        with open_fct(filepath, mode) as reader:
            first_line = reader.readline()
            if not first_line.startswith("fixedStep") and not first_line.startswith("variableStep"):
                return False
        # Check begin content
        is_valid = True
        try:
            with WigIO(filepath) as reader:
                record_idx = 0
                for record in reader:
                    if record_idx == 0:
                        record
                    if record_idx >= 10:
                        break
                    record_idx += 1
        except Exception:
            is_valid = False
        return is_valid

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record: it is not a step definition.

        :param line: The evaluated line.
        :type line: str.
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        if self.current_line.startswith("fixedStep"):
            chrom, pos, step, span = getInfoFromFixedDef(self.current_line)
            self.step = step
            self.span = span
            self._prev = {"chrom": chrom, "pos": pos - self.step}
            return False
        if self.current_line.startswith("variableStep"):
            chrom, pos, step, span = getInfoFromVariableDef(self.current_line)
            self.step = None
            self.span = span
            self._prev = {"chrom": chrom}
            return False
        if self.current_line.startswith("#") or self.current_line.startswith("browser") or self.current_line.startswith("track"):
            return False
        return True

    def write(self, chrom, pos, val):
        """
        Write record line in file.

        :param chrom: Chromosome.
        :type chrom: str
        :param pos: Position (1-based).
        :type pos: int
        :param val: Value for the position.
        :type val: float
        """
        if self.step:  # fixedStep
            if chrom != self._prev["chrom"] or self.step != self._prev["step"] or pos != self._prev["pos"] + self.step or self.span != self._prev["span"] or not self._prev["step"]:
                line = "fixedStep chrom={} start={} step={}".format(chrom, pos, self.step)
                if self.span and self.span != 1:
                    line += " span={}".format(self.span)
                self.file_handle.write(
                    "{}\n".format(line)
                )
                self.current_line_nb += 1
            self.file_handle.write("{}\n".format(val))
        else:  # variableStep
            if chrom != self._prev["chrom"] or self.span != self._prev["span"] or self._prev["step"]:
                line = "variableStep chrom={}".format(chrom)
                if self.span and self.span != 1:
                    line += " span={}".format(self.span)
                self.file_handle.write(
                    "{}\n".format(line)
                )
                self.current_line_nb += 1
            self.file_handle.write("{} {}\n".format(pos, val))
        self.current_line_nb += 1
        self._prev = {"chrom": chrom, "pos": pos, "step": self.step, "span": self.span}

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line. Value is in elt.annot["value"]
        :rtype: anacore.region.Region
        """
        if self.step:  # fixedStep
            pos = self._prev["pos"] + self.step
            value = float(self.current_line.rstrip())
            self._prev["pos"] = pos
        else:  # variableStep
            pos_str, val_str = self.current_line.rstrip().split()
            pos = int(pos_str)
            value = float(val_str)
        return Region(
            reference=self._prev["chrom"],
            start=pos,
            end=pos + self.span - 1,
            annot={"value": value}
        )
