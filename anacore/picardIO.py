# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing picard tools outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


def castCol(values, col_type):
    """
    Cast all the values to the column type.

    :param values: Values in the column.
    :type val: list
    :param col_type: Column type among "float" or "int" or None or "str".
    :type col_type: str
    """
    cast_fct = eval(col_type)
    for curr_idx, curr_val in enumerate(values):
        if curr_val == "":
            values[curr_idx] = None
        else:
            values[curr_idx] = cast_fct(curr_val)


def getColType(values):
    """
    Return the column type from these values.

    :param values: Values in the column.
    :type val: list
    :return: Type of the value among "float" or "int" or None or "str".
    :rtype: str
    """
    types = {getValType(curr_val) for curr_val in values if curr_val != "" and curr_val is not None}
    col_type = "int" if len(types) != 0 else None
    if "str" in types:
        col_type = "str"
    elif "float" in types:
        col_type = "float"
    return col_type


def getValType(val):
    """
    Return type of the value.

    :param val: Evaluated value.
    :type val: *
    :return: Type of the value among "float" or "int" or None or "str".
    :rtype: str
    """
    val_type = None
    if val is not None:
        val_type = "str"
        if val.isdigit():
            val_type = "int"
        elif val.replace(".", "", 1).isdigit():
            val_type = "float"
    return val_type


class PicardReader(object):
    """Reader for picard tools outputs."""

    def __init__(self, filepath):
        """
        Build and return an instance of PicardReader.

        :param filepath: The filepath.
        :type filepath: str
        :return: The new instance.
        :rtype: anacore.picardIO.PicardReader
        """
        self.filepath = filepath
        self.command = None
        self.metrics = None
        self.histogram = None
        self._parse()

    def _parse(self):
        """
        Read file and store his information.

        :param FH_in: File handle opens to the first line of the picard's file.
        :type FH_in: TextIOWrapper
        """
        with open(self.filepath) as FH_in:
            for curr_line in FH_in:
                if curr_line.startswith("## METRICS CLASS"):
                    self._parseMetrics(FH_in)
                elif curr_line.startswith("## HISTOGRAM"):
                    self._parseHistogram(FH_in)
                elif curr_line.startswith("## "):
                    self._parseHeader(FH_in)

    def _parseHeader(self, FH_in):
        """
        Read header in file and store his information in self.command.

        :param FH_in: File handle opens to the first line of the picard's file.
        :type FH_in: TextIOWrapper
        """
        self.command = FH_in.readline()[1:].strip()
        curr_line = FH_in.readline()
        while curr_line is not None and curr_line.strip() != "":  # Same section
            curr_line = FH_in.readline()

    def _parseHistogram(self, FH_in):
        """
        Read histogram section in file and store his information in self.histogram.

        :param FH_in: File handle opens to the histogram section of the picard's file: after tag "## HISTOGRAM".
        :type FH_in: TextIOWrapper
        """
        self.histogram = dict()
        # Parse header
        header_titles = [elt.strip() for elt in FH_in.readline().split("\t")]
        titles_by_idx = {}
        for curr_idx, curr_title in enumerate(header_titles):
            titles_by_idx[curr_idx] = curr_title
            self.histogram[curr_title] = list()
        # Get values
        curr_line = FH_in.readline()
        while curr_line is not None and curr_line.strip() != "":  # Same section
            for curr_idx, curr_val in enumerate(curr_line.split("\t")):
                col_name = titles_by_idx[curr_idx]
                col_val = curr_val.strip()
                self.histogram[col_name].append(col_val)
            curr_line = FH_in.readline()
        # Apply types
        for curr_title in header_titles:
            col_type = getColType(self.histogram[curr_title])
            if col_type != "str" and col_type is not None:
                castCol(self.histogram[curr_title], col_type)

    def _parseMetrics(self, FH_in):
        """
        Read metrics section in file and store his information in self.metrics.

        :param FH_in: File handle opens to the metrics section of the picard's file: after tag "## METRICS CLASS".
        :type FH_in: TextIOWrapper
        """
        self.metrics = list()
        # Parse header
        header_titles = [elt.strip() for elt in FH_in.readline().split("\t")]
        titles_by_idx = {curr_idx: curr_title for curr_idx, curr_title in enumerate(header_titles)}
        # Get values
        curr_line = FH_in.readline()
        while curr_line is not None and curr_line.strip() != "":  # Same section
            lib_metrics = {curr_title: None for curr_title in header_titles}
            for curr_idx, curr_val in enumerate(curr_line.split("\t")):
                col_name = titles_by_idx[curr_idx]
                col_val = curr_val.strip()
                lib_metrics[col_name] = col_val
            self.metrics.append(lib_metrics)
            curr_line = FH_in.readline()
        # Apply types
        for curr_title in header_titles:
            col_type = getColType([elt[curr_title] for elt in self.metrics])
            if col_type != "str" and col_type is not None:
                cast_fct = eval(col_type)
                for library in self.metrics:
                    if library[curr_title] == "":
                        library[curr_title] = None
                    else:
                        library[curr_title] = cast_fct(library[curr_title])
