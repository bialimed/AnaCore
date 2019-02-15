#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sv import HashedSVIO


class MAFIO(HashedSVIO):
    """Class to read and write a file in Mutation Annotation Format (see https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format)."""

    def __init__(self, filepath, mode="r"):
        """
        Build and return an instance of MAFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: MAFIO
        """
        self.version = "2.4"
        super().__init__(filepath, mode, "\t", "")

    def _parseHeader(self):
        """Parse MAF header to set the attribute titles."""
        if self.current_line_nb == 0:
            self.current_line = self.file_handle.readline()
            self.version = self.current_line.split(" ")[1]
            self.current_line_nb += 1
            while not self.current_line.startswith("Hugo_Symbol"):
                self.current_line = self.file_handle.readline().rstrip()
                self.current_line_nb += 1
        self.titles = [elt.strip() for elt in self.current_line.split(self.separator)]

    def write(self, record):
        """
        Write record line in file.

        :param record: The record.
        :type record: dict
        """
        if self.current_line_nb == 0 and self.titles is not None:
            self.file_handle.write(
                "#version {}\n{}\n".format(
                    self.version,
                    self.separator.join(self.titles)
                )
            )
            self.current_line_nb += 2
        self.file_handle.write(self.recordToLine(record) + "\n")
        self.current_line_nb += 1
