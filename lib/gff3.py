#
# Copyright (C) 2015 INRA
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

import sys, re

class GFF3Record:
    """
     @summary : Record for GFF3.
    """
    def __init__( self ):
        self.seq_id = None
        self.source = None
        self.type = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.phase = None
        self.attributes = None

    def setAttribute( self, tag, value ):
        """
         @summary : Create or replace an attribute tag.
          @param tag : tag of the attribute.
          @param value : value of the attribute tag.
        """
        cleaned_tag = GFF3Record._getCleanedAttribute(tag)
        cleaned_value = GFF3Record._getCleanedAttribute(value)
        if self.attributes is not None :
            self.attributes[cleaned_tag] = cleaned_value
        else:
            raise ValueError("The attibute 'Attributes' is not initialized.")

    def addToAttribute( self, tag, value ):
        """
         @summary : Add one value on an existing tag.
          @param tag : tag of the attribute.
          @param value : value to add of the tag.
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
         @summary : Returns a string in GFF3 format attributes field from the GFF3Record.attributes.
         @return : [str] the attributes in GFF3 format.
        """
        gff_string = ""
        for tag in self.attributes:
            gff_string = gff_string + tag + "=" + str(self.attributes[tag]) + ";"

        return gff_string[:-1]

    def toGff( self ):
        """
         @summary : Returns a string in GFF3 format from the GFF3Record object.
         @return : [str] the line in GFF3 format.
        """
        gff_record = "\t".join( [self.seq_id, self.source, self.type, str(self.start), str(self.end), str(self.score), self.strand, str(self.phase), self._attributesToGff()] )

        return gff_record

    def attributesToStr( self, tag ):
        """
         @summary : Returns the attribute value in human readable format.
          @param tag : [str] the attribute tag.
         @return : [str] the human readable value.
         @see : RFC 3986 Percent-Encoding
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
         @summary : Returns value after GFF3 attribute cleaning. cleanning :
            - URL escaping rules are used for tags or values containing the following characters: ",=;".
            - Spaces are allowed in this field, but tabs must be replaced with the space.
            - Quotes ' and " are deleted.
          @param dirty_value : [str] value before cleaning.
         @return : [str] the clean value.
         @see : RFC 3986 Percent-Encoding
        """
        cleaned_value = dirty_value.replace(';', '%3B')
        cleaned_value = cleaned_value.replace(',', '%2C')
        cleaned_value = cleaned_value.replace('=', '%3D')
        cleaned_value = cleaned_value.replace('\t', ' ')
        cleaned_value = cleaned_value.replace("'", '')
        cleaned_value = cleaned_value.replace('"', '')

        return cleaned_value

    @staticmethod
    def fromGff( line ):
        """
         @summary : Returns a GFF3Record from a GFF3 line.
          @param line : line of the GFF.
         @return : [GFF3Record] the record.
        """
        gff_record = GFF3Record()
        line_fields = line.split("\t")
        gff_record.seq_id = line_fields[0]
        gff_record.source = line_fields[1]
        gff_record.type = line_fields[2]
        gff_record.start = int(line_fields[3])
        gff_record.end = int(line_fields[4])
        gff_record.score = line_fields[5]
        gff_record.strand = line_fields[6]
        gff_record.phase = line_fields[7]
        # Parse attributes
        gff_record.attributes = dict()
        attributes = "\t".join(line_fields[8:])
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
            tag = matches.group(1)
            values = matches.group(2).split(',')
            for current_val in values:
                gff_record.addToAttribute(tag, current_val)
        return gff_record


class GFF3IO:
    """
     @summary : Specific handler for GFF3 file.
    """
    def __init__( self, file_path, mode="r" ):
        self._path = file_path
        self._handle = open( file_path, mode )
        self._line = 1

    def __del__( self ):
        self.close()

    def __iter__( self ):
        for line in self._handle:
            line = line.rstrip()
            self._line += 1
            if line.startswith('#') :
                continue
            try:
                gff_record = GFF3Record.fromGff(line)
            except:
                raise IOError( "The line " + str(self._line) + " in '" + self._path + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + line )
            else:
                yield gff_record

    def close( self ) :
        if  hasattr(self, '_handle') and self._handle is not None:
            self._handle.close()
            self._handle = None
            self._line = None

    def write( self, gff_record ):
        """
         @summary : Write one line on gff file.
          @param gff_record : [GFF3Record] the object to write.
        """
        self._handle.write( gff_record.toGff() + "\n" )
        self._line += 1
