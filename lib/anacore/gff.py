# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing/processing GFF files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import copy
from anacore.abstractFile import AbstractFile, isEmpty
from anacore.region import Region


class GFF3Record(Region):
    """Class to manage a record comind from a GFF3."""

    def __init__(self, seq_id=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=None, attributes=None):
        """
        Build and return an instance of GFF3Record.

        :param seq_id: The ID of the landmark used to establish the coordinate system for the current feature.
        :type seq_id: str.
        :param source: The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type that is a subclass of the type in the type column.
        :type source: str.
        :param type: The type of the feature (previously called the "method"). This is constrained to be either: (a) a term from the "lite" sequence ontology, SOFA; or (b) a SOFA accession number. The latter alternative is distinguished using the syntax SO:000000.
        :type type: str.
        :param start: The start of the feature, in 1-based integer coordinates, relative to the landmark given in seq_id. Start is always less than or equal to end.
        :type start: int.
        :param end: The end of the feature, in 1-based integer coordinates, relative to the landmark given in seq_id. Start is always less than or equal to end.
        :type end: int.
        :param score: The score of the feature. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features.
        :type score: float.
        :param strand: The strand of the instance ("+" or "-").
        :type strand: str.
        :param phase: For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3. For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field.
        :type phase: int.
        :param attributes: The annotations of the feature.
        :type attributes: dict.
        :return: The new instance.
        :rtype: GFF3Record
        """
        name = None
        cleaned_attributes = attributes
        if attributes is not None and "Name" in attributes:
            name = attributes["Name"]
            cleaned_attributes = copy.deepcopy(attributes)
            del(cleaned_attributes["Name"])
        Region.__init__(self, start, end, strand, seq_id, name, cleaned_attributes)
        self.source = source
        self.type = type
        self.score = score
        self.phase = phase

    def __getattribute__(self, attr_name):
        """
        Return the value of the selected attribute.

        :param attr_name: The name of the selected attribute.
        :type attr_name: str
        :return: The value of the attribute.
        :rtype: *
        """
        attr_value = None
        if attr_name == "seq_id":
            attr_value = object.__getattribute__(self, "reference").name
        else:
            attr_value = object.__getattribute__(self, attr_name)
        return attr_value

    def __setattr__(self, attr_name, attr_value):
        """
        Set the value of the selected attribute.

        :param attr_name: The name of the selected attribute.
        :type attr_name: str
        :param attr_value: The new value for the selected attribute.
        :type attr_value: str
        """
        if attr_name == "seq_id":
            self.setReference(attr_value)
        else:
            object.__setattr__(self, attr_name, attr_value)

    def addToAnnot(self, tag, value):
        """
        Add one value on an existing tag.

        :param tag: tag of the annottation.
        :type tag: str.
        :param value: value to add.
        :type value: str.
        """
        if tag in self.annot:
            self.annot[tag] = self.annot[tag] + "," + value
        else:
            self.annot[tag] = value

    def _annotToGff(self):
        """
        Return a string in GFF3 field format corresponding to the attributes of the record.

        :return: The record.
        :rtype: str
        """
        gff_str = ""
        for tag, value in sorted(self.annot.items()):
            gff_str += "{}={};".format(GFF3IO.encodedValue(tag), GFF3IO.encodedValue(value))
        if gff_str != "":
            gff_str = gff_str[:-1]
        return gff_str

    def toGff(self):
        """
        Return the string corrsponding to the GFF3Record in GFF3 format.

        :return: Line in GFF3 format.
        :rtype: str
        """
        attributes = self._annotToGff()
        if self.name is not None:
            name = "{}={}".format("Name", GFF3IO.encodedValue(self.name))
            attributes = name + ";" + attributes
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.seq_id,
            ("." if self.source is None else self.source),
            self.type,
            self.start,
            self.end,
            ("." if self.score is None else self.score),
            ("." if self.strand is None else self.strand),
            ("." if self.phase is None else self.phase),
            attributes
        )


class GFF3IO(AbstractFile):
    """Class to manage read and write in GFF3 file."""

    def __init__(self, filepath, mode="r", separator="\t", title_starter="#", has_title=True):
        """
        Build and return an instance of GFF3IO.

        :param filepath: The filepath.
        :type filepath: str.
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str.
        :return: The new instance.
        :rtype: GFF3IO
        """
        # Convert mode for append in empty file
        if mode == "a" and isEmpty(filepath):
            mode = "w"
        # Initialise instance
        AbstractFile.__init__(self, filepath, mode)
        # Add header
        if mode == "w":
            self.file_handle.write("##gff-version 3\n")
            self.current_line_nb += 1

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record (it is not a comment or an header line).

        :param line: The evaluated line.
        :type line: str
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        is_record = True
        if line.startswith("#"):
            is_record = False
        return is_record

    def _parseLine(self):
        """
        Return a structured record from the GFF current line.

        :return: The record.
        :rtype: GFF3Record
        """
        fields = [elt.strip() for elt in self.current_line.split("\t")]
        gff_record = GFF3Record(
            fields[0],
            (None if fields[1] == "." else fields[1]),
            fields[2],
            (None if fields[3] == "." else int(fields[3])),
            (None if fields[4] == "." else int(fields[4])),
            (None if fields[5] == "." else fields[5]),
            (None if fields[6] == "." else fields[6]),
            (None if fields[7] == "." else fields[7]),
        )
        # Parse attributes
        attributes_field = "\t".join(fields[8:])
        if attributes_field.strip().endswith(";"):  # if attributes end with ';'
            attributes_field = attributes_field.strip()[:-1]
        attributes_array = attributes_field.split(";")
        cleaned_attributes = list()
        for attribute in attributes_array:
            if "=" not in attribute:
                cleaned_attributes[len(cleaned_attributes) - 1] += ";" + attribute
            else:
                cleaned_attributes.append(attribute)
        for attribute in cleaned_attributes:
            matches = re.match("^([^=]+)=(.*)", attribute)
            tag = GFF3IO.decodedValue(matches.group(1).strip())
            values = [GFF3IO.decodedValue(elt) for elt in matches.group(2).split(',')]
            for current_val in values:
                if tag == "Name":
                    gff_record.name = current_val
                else:
                    gff_record.addToAnnot(tag, current_val)
        return gff_record

    def write(self, record):
        """
        Write one line on GFF file.

        :param record: The object to write.
        :type record: GFF3Record
        """
        self.file_handle.write(record.toGff() + "\n")
        self.current_line_nb += 1

    @staticmethod
    def decodedValue(encoded_value):
        """
        Return the attribute tag/value in human readable format.

        :param encoded_value: The attribute tag or value to decode.
        :type encoded_value: str
        :return: Human readable value.
        :rtype: str
        :see: RFC 3986 Percent-Encoding.
        """
        decoded_value = encoded_value.replace('%3B', ';').replace('%2C', ',').replace('%3D', '=')
        return decoded_value

    @staticmethod
    def encodedValue(dirty_value):
        """
        Return value after GFF3 attribute cleaning. Cleanning:

            - URL escaping rules are used for tags or values containing the following characters: ",=;".
            - Spaces are allowed in this field, but tabs must be replaced by space.
            - Quotes ' and " are deleted.

        :param dirty_value: The attribute tag or value to encode.
        :type dirty_value: str
        :return: GFF3 encoded value.
        :rtype: str
        :see: RFC 3986 Percent-Encoding.
        """
        cleaned_value = dirty_value.replace(';', '%3B')
        cleaned_value = cleaned_value.replace(',', '%2C')
        cleaned_value = cleaned_value.replace('=', '%3D')
        cleaned_value = cleaned_value.replace('\t', ' ')
        cleaned_value = cleaned_value.replace("'", '')
        cleaned_value = cleaned_value.replace('"', '')
        return cleaned_value
