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
__version__ = '1.1.0'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

import gzip


def is_gzip(file):
    """
    @return: [bool] True if the file is gziped.
    @param file: [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open( file )
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip


class AbstractFile(object):
    def __init__( self, filepath, mode="r" ):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode + "t" )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 1
        self.current_line = None

    def close( self ) :
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def closed( self ):
        return self.file_handle.closed()

    def __del__( self ):
        self.close()

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def isRecordLine(self, line):
        """
        @summary: Returns True if the line corresponds to a record (it is not a comment or an header line).
        @param line: [str] The evaluated line.
        @return: [bool] True if the line corresponds to a record.
        """
        return True

    def __iter__( self ):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            if not self.isRecordLine(self.current_line): # Skip non-record line
                continue
            record = None
            try:
                record = self._parseLine()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
            else:
                yield record

    def read(self):
        """
        @summary: Returns all the record from the file.
        @return: [list] The list of records.
        """
        return [record for record in self]
