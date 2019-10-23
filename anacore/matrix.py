# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing/processing files containing 2D matrix."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import numpy as np


class DistanceMatrixIO:
    """Class to read and write a distance matrix from or to separate value file."""

    def __init__(self, path, mode="r", separator="\t", cast_fct=float):
        """
        Build and return an instance of DistanceMatrixIO.

        :param path: Path to the distances file.
        :type path: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param separator: Separator used between values.
        :type separator: str
        :param cast_fct: Function used to cast values in matrix.
        :type cast_fct: fct
        :return: The new instance.
        :rtype: matrix.DistanceMatrixIO
        """
        self.cast_fct = cast_fct
        self.dist_matrix = list()
        self.filepath = path
        self.more = mode
        if mode != "r":
            raise NotImplementedError(
                "Only read mode is currently implemented for {}.".format(
                    self.__class__.__name__
                )
            )
        self.names = None
        self.separator = separator
        self._parse()

    def _parse(self):
        """
        Parse file to load data in instance.

        :todo: Select parser by input (_parseDense(1), _parseDense(0), _parseTriangular(), _parseSparse()).
        """
        self._parseDense()

    def _parseDense(self, col_offset=1):
        """
        Parse file in dense format to load data in instance.

        :param col_offset: Skip this number of columns from the left of file. This option is used to skip colnames.
        :param col_offset: int
        """
        with open(self.filepath) as FH_in:
            self.names = [elt.strip() for elt in FH_in.readline().split(self.separator)]
            self.names = self.names[1:]
            for row in FH_in:
                row_fields = list()
                for elt in row.split(self.separator)[col_offset:]:
                    elt = elt.strip()
                    if elt == "":
                        elt = np.nan
                    else:
                        elt = self.cast_fct(elt)
                    row_fields.append(elt)
                self.dist_matrix.append(row_fields)
        self.dist_matrix = np.array(self.dist_matrix)
