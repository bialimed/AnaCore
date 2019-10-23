# -*- coding: utf-8 -*-
"""Classes and funcions for reading/writing text files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import gzip


def isEmpty(path):
    """
    Return True if the file does not exists or is empty.

    :param path: Path to processed file.
    :type path: str
    :return: True if the file does not exists or is empty.
    :rtype: bool
    """
    is_empty = True
    if os.path.exists(path):
        if not isGzip(path):
            is_empty = os.path.getsize(path) == 0
        else:
            with gzip.open(path, 'rt') as FH:
                is_empty = FH.seek(0, 2) == 0
    return is_empty


def isGzip(path):
    """
    Return True if the file is gziped.

    :param path: Path to processed file.
    :type path: str
    :return: True if the file is gziped.
    :rtype: bool
    """
    is_gzip = None
    FH_input = gzip.open(path)
    try:
        FH_input.readline()
        is_gzip = True
    except Exception:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip


class AbstractFile(object):
    """Class inherited by classes to read and write text file."""

    def __init__(self, filepath, mode="r"):
        """
        Build and return an instance of AbstractFile.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: abstractFile.AbstractFile
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and isGzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 0
        self.current_line = None

    def close(self):
        """Closes file handle."""
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def closed(self):
        """
        Return True if the file is closed.

        :retrun: True if the file is closed.
        :rtype: bool
        """
        return self.file_handle.closed()

    def __del__(self):
        self.close()

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record (it is not a comment or an header line).

        :param line: The evaluated line.
        :type line: str.
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        return True

    def __iter__(self):
        for line in self.file_handle:
            self.current_line = line.rstrip("\n")
            self.current_line_nb += 1
            if not self.isRecordLine(self.current_line):  # Skip non-record line
                continue
            record = None
            try:
                record = self._parseLine()
            except Exception:
                raise IOError('The line {} in "{}" cannot be parsed by {}. Line content: {}'.format(self.current_line_nb, self.filepath, self.__class__.__name__, self.current_line))
            else:
                yield record

    def read(self):
        """
        Return all the record from the file.

        :return: The list of records.
        :rtype: list
        """
        return [record for record in self]
