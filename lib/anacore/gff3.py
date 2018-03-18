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
__version__ = '1.2.0'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

import re
from anacore.abstractFile import AbstractFile


class GFF3Record:
    """
    @summary: Record for GFF3.
    """
    def __init__( self, seq_id=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=None, attributes=None ):
        self.seq_id = seq_id
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def setAttribute( self, tag, value ):
        """
        @summary: Create or replace an attribute tag.
        @param tag: tag of the attribute.
        @param value: value of the attribute tag.
        """
        cleaned_tag = GFF3Record._getCleanedAttribute(tag)
        cleaned_value = GFF3Record._getCleanedAttribute(value)
        if self.attributes is not None :
            self.attributes[cleaned_tag] = cleaned_value
        else:
            raise ValueError("The attribute 'Attributes' is not initialized.")

    def addToAttribute( self, tag, value ):
        """
        @summary: Add one value on an existing tag.
        @param tag: tag of the attribute.
        @param value: value to add of the tag.
        """
        cleaned_tag = GFF3Record._getCleanedAttribute(tag)
        cleaned_value = GFF3Record._getCleanedAttribute(value)
        if self.attributes is not None :
            if cleaned_tag in self.attributes:
                self.attributes[cleaned_tag] = self.attributes[cleaned_tag] + "%2C"  + cleaned_value
            else:
                self.attributes[cleaned_tag] = cleaned_value
        else:
            raise ValueError("The attibute 'Attributes' is not initialized.")

    def _attributesToGff( self ):
        """
        @summary: Returns a string in GFF3 format attributes field from the GFF3Record.attributes.
        @return: [str] the attributes in GFF3 format.
        """
        gff_string = ""
        for tag in self.attributes:
            gff_string = gff_string + tag + "=" + str(self.attributes[tag]) + ";"

        return gff_string[:-1]

    def toGff( self ):
        """
        @summary: Returns a string in GFF3 format from the GFF3Record object.
        @return: [str] the line in GFF3 format.
        """
        gff_record = "\t".join( [self.seq_id, self.source, self.type, str(self.start), str(self.end), str(self.score), self.strand, str(self.phase), self._attributesToGff()] )

        return gff_record

    def attributesToStr( self, tag ):
        """
        @summary: Returns the attribute value in human readable format.
        @param tag: [str] the attribute tag.
        @return: [str] the human readable value.
        @see: RFC 3986 Percent-Encoding
        """
        cleaned_tag = GFF3Record._getCleanedAttribute(tag)
        if cleaned_tag in self.attributes:
            readable_value = self.attributes[cleaned_tag].replace('%3B', ';')
            readable_value = readable_value.replace('%2C', ',')
            redable_value = readable_value.replace('%3D', '=')
            return redable_value
        else:
            raise ValueError("The attibute 'Attributes' is not initialized.")

    @staticmethod
    def _getCleanedAttribute( dirty_value ):
        """
        @summary: Returns value after GFF3 attribute cleaning. cleanning :
            - URL escaping rules are used for tags or values containing the following characters: ",=;".
            - Spaces are allowed in this field, but tabs must be replaced with the space.
            - Quotes ' and " are deleted.
        @param dirty_value: [str] value before cleaning.
        @return: [str] the clean value.
        @see: RFC 3986 Percent-Encoding
        """
        cleaned_value = dirty_value.replace(';', '%3B')
        cleaned_value = cleaned_value.replace(',', '%2C')
        cleaned_value = cleaned_value.replace('=', '%3D')
        cleaned_value = cleaned_value.replace('\t', ' ')
        cleaned_value = cleaned_value.replace("'", '')
        cleaned_value = cleaned_value.replace('"', '')

        return cleaned_value


class GFF3IO(AbstractFile):
    def isRecordLine(self, line):
        """
        @summary: Returns True if the line corresponds to a record (it is not a comment or an header line).
        @param line: [str] The evaluated line.
        @return: [bool] True if the line corresponds to a record.
        """
        is_record = True
        if line.startswith("#"):
            is_record = False
        return is_record

    def _parseLine(self):
        """
        @summary: Returns a structured record from the GFF current line.
        @return: [GFF3Record] The record.
        """
        fields = self.current_line.split("\t")
        gff_record = GFF3Record(
            fields[0],
            fields[1],
            fields[2],
            ("." if fields[3] == "." else int(fields[3])),
            ("." if fields[4] == "." else int(fields[4])),
            fields[5],
            fields[6],
            ("." if fields[7] == "." else int(fields[7])),
        )
        # Parse attributes
        gff_record.attributes = dict()
        attributes = "\t".join(fields[8:])
        if attributes.strip().endswith(";"): # if attributes end with ';'
            attributes = attributes.strip()[:-1]
        attributes_array = attributes.split(";")
        cleaned_attributes = list()
        for attribute in attributes_array:
            if not "=" in attribute:
                cleaned_attributes[len(cleaned_attributes)-1] += " %3B " + attribute
            else:
                cleaned_attributes.append(attribute)
        for attribute in cleaned_attributes:
            matches = re.match("^([^=]+)=(.*)", attribute)
            tag = matches.group(1).strip()
            values = matches.group(2).split(',')
            for current_val in values:
                gff_record.addToAttribute(tag, current_val)
        return gff_record

    def write( self, gff_record ):
        """
        @summary: Write one line on gff file.
        @param gff_record: [GFF3Record] the object to write.
        """
        self._handle.write( gff_record.toGff() + "\n" )
        self._line += 1
