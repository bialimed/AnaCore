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


class SVIO:
    def __init__(self, filepath, mode="r", separator="\t"):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        @param separator: [str] Separator used between values ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.separator = separator
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 1
        self.current_line = None
        self.titles = None

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __iter__(self):
        try:
            for line in self.file_handle:
                self.current_line_nb += 1
                self.current_line = line
                if line.startswith("#"):
                    self.titles = [field.strip() for field in line[1:].split( self.separator )]
                else:
                    yield( self._parse_line(self.current_line) )
        except:
            raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + "." )

    def next(self):
        """
        @summary: Returns the next record.
        @return: [list] The next record.
        """
        self.current_line_nb += 1
        self.current_line = FH_in.file_handle.next()
        return( self._parse_line(self.current_line) )
    
    #~ def _parse_line(self, line): ###################### A surcharger en version dict
        #~ fields = [field.strip() for field in line.split( self.separator )]
        #~ return fields

    def _parse_line(self, line): ###################### A surcharger en version dict
        fields = [field.strip() for field in line.split( self.separator )]
        record = {self.titles[idx]:val for idx, val in enumerate(fields)}
        return record

    @staticmethod
    def is_valid(filepath, separator="\t"):
        is_valid = False
        FH_in = SVIO(filepath, separator=separator)
        try:
            line_idx = 0
            line = FH_in.file_handle.readline()
            while line_idx < 10 and line:
                line = FH_in.file_handle.readline()
                line_idx += 1
            is_valid = True
        except:
            pass
        finally:
            FH_in.close()
        return is_valid

    def write(self, record):
        self.file_handle.write( self.recordToLine(record) + "\n" )

    def recordToLine(self, record):
        """
        @summary : Returns the record in SV format.
        @param record : [list or dict] The record to process.
        @return : [str] The SV line corresponding to the record.
        """
        #################################check if record is iterable and if has key
        line = self.separator.join(record)
        return( line )
