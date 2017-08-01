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
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

from abstractFile import AbstractFile


class BEDRecord:
    def __init__(self, chrom=None, start=None, end=None, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        """
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts 


class BEDIO(AbstractFile):
    def BEDRecordToBEDLine(self, bed_record):
        """
        @summary : Returns the record in BED format.
        @param bed_record : [BEDRecord] The sequence to process.
        @return : [str] The BED line.
        """
        line = "{}\t{}\t{}".format(
            self.chrom,
            (self.start - 1),
            self.end
        )
#             ("." if self.name is None else self.name),
#             ("." if self.score is None else str(self.score)),
#             ("." if self.strand is None else self.strand),
#             ("." if self.thickStart is None else str(self.thickStart - 1)),
#             ("." if self.thickEnd is None else str(self.thickEnd)),
#             ("." if self.itemRgb is None else ",".join(self.itemRgb)),
#             ("." if self.blockCount is None else str(self.blockCount)),
#             ("." if self.blockSizes is None else ",".join(map(str, self.blockSizes))),
#             ("." if self.blockStarts is None else ",".join(map(str, self.blockStarts)))
        return line

    def isRecordLine(self, line):
        """
        @summary: Returns True if the line corresponds to a record (it is not a comment or an header line).
        @param line: [str] The evaluated line.
        @return: [bool] True if the line corresponds to a record.
        """
        is_record = True
        if line.startswith("browser ") or line.startswith("track ") or line.startswith("#"):
            is_record = False
        return is_record

    @staticmethod
    def isValid(filepath):
        """
        @summary: Returns True is the file can be a BED.
        @return: [bool] True if the file can be a BED.
        """
        is_valid = False
        try:
            with BEDIO(filepath) as FH_in:
                record_idx = 0
                for record in FH_in:
                    if record_idx >= 10:
                        break
                    if record.strand is not None and record.strand not in ["+", "-", "."]:
                        raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
                    #~ if record.itemRgb is not None and len(record.itemRgb.split(",")) != 3:
                        #~ raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               #~ "Line content : " + self.current_line )                        
                    record_idx += 1
                is_valid = True
        except:
            pass
        return is_valid

    def _parseLine(self):
        """
        @summary: Returns a structured record from the BED current line.
        @return: [BEDrecord] The record.
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        fields[1] = int(fields[1]) + 1 # Start in BED is 0-based
        fields[2] = int(fields[2])
        if len(fields) >= 5:
            fields[4] = None if fields[4] == "." else int(fields[4]) # A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
            if len(fields) >= 7:
                fields[6] = int(fields[6]) + 1 # The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
                if len(fields) >= 8:
                    fields[7] = int(fields[7]) # The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
                    if len(fields) >= 10:
                        fields[9] = int(fields[9]) # The number of blocks (exons) in the BED line.
                        if len(fields) >= 11:
                            fields[10] = [int(block) for block in fields[10].split(",")] # A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
                            if len(fields) >= 12:
                                fields[11] = [int(block) for block in fields[11].split(",")] # A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
        return BEDRecord(*fields)

    def write(self, bed_record):
        """
        @summary: Writes record line in file.
        @param bed_record: [BEDRecord] The record.
        """
        self.file_handle.write( self.BEDRecordToBEDLine(bed_record) + "\n" )
        self.current_line_nb += 1