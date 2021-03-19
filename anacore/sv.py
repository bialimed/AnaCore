# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing Separated Values files (CSV, TSV, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.abstractFile import isEmpty, AbstractFile


class SVIO(AbstractFile):
    """Class to read and write in separated value files like TSV, CSV, etc. Each rows is view as a list."""

    def __init__(self, filepath, mode="r", separator="\t", title_starter=None, has_title=True, metadata_starter="##"):
        """
        Build and return an instance of SVIO.

        :param filepath: The filepath.
        :type filepath: str.
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str.
        :param separator: Separator used between values.
        :type separator: str.
        :param title_starter: The string used to introduce the title line.
        :type title_starter: str.
        :param has_title: If the first line contains the titles of columns.
        :type has_title: bool.
        :param title_metadata: The string used to introduce a metadata line.
        :type title_metadata: str.
        :return: The new instance.
        :rtype: SVIO
        """
        # Convert mode for append in empty file
        if mode == "a" and isEmpty(filepath):
            mode = "w"
        # Get existing titles and metadata if the file is in append mode
        pre_nb_line = 0
        pre_titles = None
        pre_metadata = []
        if mode == "a":  # Get titles from existing file
            with SVIO(filepath, "r", separator, title_starter, has_title) as reader:
                pre_titles = reader.titles
                pre_metadata = reader.metadata
                pre_nb_line = reader.current_line_nb + sum([1 for line in reader.file_handle])
        # Initialise instance
        AbstractFile.__init__(self, filepath, mode)
        self.separator = separator
        self.metadata = pre_metadata
        self.metadata_starter = metadata_starter
        self._has_title = has_title
        self.titles = pre_titles
        self.title_starter = title_starter
        if mode == "r":
            self._parseHeader()
        elif mode == "a":
            self.current_line_nb = pre_nb_line

    def _parseHeader(self):
        """Parse SV header to set the attributes titles and metadata."""
        if self.current_line_nb == 0 and not isEmpty(self.filepath):
            # Metadata
            curr_line_pointer = self.file_handle.tell()
            self.current_line = self.file_handle.readline().rstrip()
            self.current_line_nb += 1
            while self.current_line.startswith(self.metadata_starter):
                self.metadata.append(self.current_line[len(self.metadata_starter):])
                curr_line_pointer = self.file_handle.tell()
                self.current_line = self.file_handle.readline().rstrip()
                self.current_line_nb += 1
            # Titles
            if self._has_title:
                cleaned_line = self.current_line
                if self.title_starter is not None and self.title_starter != "":
                    if not self.current_line.startswith(self.title_starter):
                        raise Exception('The title line does not starts with "{}".'.format(self.title_starter))
                    cleaned_line = self.current_line[len(self.title_starter):]
                self.titles = [elt.strip() for elt in cleaned_line.split(self.separator)]
            else:
                self.file_handle.seek(curr_line_pointer)
                self.current_line_nb -= 1

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record (it is not a comment or an header line).

        :param line: The evaluated line.
        :type line: str.
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        is_record = True
        if line.startswith("#"):
            is_record = False
        return is_record

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line.
        :rtype: list
        """
        fields = [field.strip() for field in self.current_line.split(self.separator)]
        return fields

    @staticmethod
    def isValid(filepath, separator="\t", metadata_starter="##"):
        """
        Return True is the file can be a SV file.

        :return: True if the file can be a SV file.
        :rtype: bool
        """
        is_valid = False
        pointer_after_metadata = False
        nb_fields = list()  # The number of fields for each row
        reader = SVIO(filepath, separator=separator, title_starter="", has_title=False, metadata_starter="##")
        try:
            # Read the 10 first lines
            line_idx = 1
            line = reader.file_handle.readline()
            if line is None:  # File is empty
                is_valid = True
            else:  # File has content
                if not line.startswith(metadata_starter):
                    pointer_after_metadata = True
                    nb_fields.append(line.count(separator))
                while line_idx < 10 and line:
                    line = reader.file_handle.readline()
                    if line:
                        if pointer_after_metadata or not line.startswith(metadata_starter):
                            pointer_after_metadata = True
                            nb_fields.append(line.count(separator))
                            line_idx += 1
                if len(set(nb_fields)) == 1 and nb_fields[0] > 1:  # All lines have the same number of columns
                    is_valid = True
        except Exception:
            pass
        finally:
            reader.close()
        return is_valid

    def writeHeader(self):
        """Write header lines in file."""
        # Metadata
        for curr_metadata in self.metadata:
            self.file_handle.write(
                self.metadata_starter + curr_metadata + "\n"
            )
            self.current_line_nb += 1
        # Titles
        if self.titles is not None:
            line_starter = "" if self.title_starter is None else self.title_starter
            self.file_handle.write(
                line_starter + self.separator.join(self.titles) + "\n"
            )
            self.current_line_nb += 1

    def write(self, record):
        """
        Write record line in file.

        :param record: The record.
        :type record: list.
        """
        # Write header if the file is empty
        if self.current_line_nb == 0:
            self.writeHeader()
        # Write record
        self.file_handle.write(self.recordToLine(record) + "\n")
        self.current_line_nb += 1

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: The record to process.
        :type record: list.
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        line = self.separator.join([str(elt) for elt in record])
        return(line)


class HashedSVIO(SVIO):
    """Class to read and write in separated value files like TSV, CSV, etc. Each rows is view as a dict indexed by columns titles."""

    def __init__(self, filepath, mode="r", separator="\t", title_starter=None, metadata_starter="##"):
        """
        Build and return an instance of HashedSVIO.

        :param filepath: The filepath.
        :type filepath: str.
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str.
        :param separator: Separator used between values.
        :type separator: str.
        :param title_starter: The string used to introduce the title line.
        :type title_starter: str.
        :param title_metadata: The string used to introduce a metadata line.
        :type title_metadata: str.
        :return: The new instance.
        :rtype: HashedSVIO
        """
        super().__init__(filepath, mode, separator, title_starter, True, metadata_starter)

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line.
        :rtype: dict
        """
        fields = super()._parseLine()
        record = {self.titles[idx]: val for idx, val in enumerate(fields)}
        return record

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: The record to process.
        :type record: dict.
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        fields = [record[curr_title] for curr_title in self.titles]
        line = self.separator.join([str(elt) for elt in fields])
        return(line)
