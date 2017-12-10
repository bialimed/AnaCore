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
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import numpy as np


class DistanceMatrixIO:
    """
    @summary: Class to read and write a distance matrix from or to separate
    value file.
    """
    def __init__(self, path, mode="r", separator="\t", cast_fct=float):
        """
        @param path: [str] Path to the distances file.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        @param separator: [str] Separator used between values.
        @param cast_fct: [fct] Function used to cast values in matrix.
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
        @summary: Parses file to load data in instance.
        @todo: Select parser by input (_parseDense(1), _parseDense(0),
        _parseTriangular(), _parseSparse()).
        """
        self._parseDense()

    def _parseDense(self, col_offset=1):
        """
        @summary: Parses file in dense format to load data in instance.
        @param col_offset: [int] Skip this number of columns from the left of
        file. This option is used to skip colnames.
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
