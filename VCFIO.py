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
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'


import re
import gzip


def is_gzip(file):
    """
    @return: [bool] True if the file is gziped.
    @param file : [str] Path to processed file.
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


class VCFRecord:
    def __init__(self, region=None, position=None, knownSNPId=None, refAllele=None, altAlleles=None, qual=None, pFilter=None, info=None, pFormat=None, samples=None):
        self.chrom   = region
        self.pos     = position
        self.id      = knownSNPId
        self.ref     = refAllele
        self.alt     = altAlleles
        self.qual    = qual
        self.filter  = pFilter
        self.info    = info if info is not None else dict()
        self.format  = pFormat if pFormat is not None else list()
        self.samples = samples if samples is not None else dict()
    
    def isIndel(self):
        """
        @summary : Return True if the variant is an insertion or a deletion.
        @return : [bool]
        """
        isIndel = False
        if self.ref == "." or self.ref == "-":
            isIndel = True
        else:
            for allele in self.alt:
                if len(allele) != len(self.ref) or allele == "." or allele == "-":
                    isIndel = True
        return isIndel
        
    def type(self):
        """
        @summary : Returns the vrariant type.
        @return : [str] 'snp' or 'indel' or 'variation'.
        """
        record_type = "snp"
        if self.isIndel():
            record_type = "indel"
        elif len(self.ref) > 1:
            record_type = "variation"
        return record_type


class VCFIO:
    def __init__( self, filepath, mode="r" ):
        """
        @param filepath : [str] The filepath.
        @param mode : [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode + "t" )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 0
        self.current_line = None
        self.samples = list()
        self.info = dict() # { "IDREP":{"type": int, "number": 1, "description": "Number of times RU is repeated in indel allele."} }
        self.format = dict() # { "GT":{"type": str, "number": 1, "description": "Genotype"} }
    
    def __del__( self ):
        self.close()
    
    def __iter__( self ):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            if self.current_line.startswith('#'):
                type_fct = {
                    "String": str,
                    "Integer": int,
                    "Float": float,
                    "Character": str,
                    "Flag": None
                } ################################################################### Flag
                if self.current_line.startswith("##INFO"):
                    line = self.current_line[8:-1] # Remove ##INFO=< and >
                    match = re.search('ID=([^,]+)', line)
                    id = match.group(1)
                    match = re.search('Type=([^,]+)', line)
                    type = match.group(1)
                    match = re.search('Description="([^"]+)"', line)
                    description = match.group(1)
                    match = re.search('Number=([^,]+)', line)
                    number = None
                    if match.group(1) != "." and match.group(1) != "A":
                        number = int(match.group(1))
                    self.info[id] = { "type": type_fct[type], "number": number, "description": description }
                elif self.current_line.startswith("##FORMAT"):
                    line = self.current_line[10:-1] # Remove ##FORMAT=< and >
                    match = re.search('ID=([^,]+)', line)
                    id = match.group(1)
                    match = re.search('Type=([^,]+)', line)
                    type = match.group(1)
                    match = re.search('Description="([^"]+)"', line)
                    description = match.group(1)
                    match = re.search('Number=([^,]+)', line)
                    number = None
                    if match.group(1) != "." and match.group(1) != "A":
                        number = int(match.group(1))
                    self.format[id] = { "type": type_fct[type], "number": number, "description": description }
                elif line.startswith("#CHROM\tPOS"):
                    self.samples = [spl.strip() for spl in line.split("\t")[9:]]
                continue
            try:
                vcf_record = self._parse_line()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
            else:
                yield vcf_record

    def close( self ) :
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def _parse_line( self ):
        """
        @summary : Returns a structured record from the VCF current line.
        @return : [VCFRecord]
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        variation = VCFRecord()
        variation.chrom  = fields[0]
        variation.pos    = int(fields[1])
        variation.id     = fields[2]
        variation.ref    = fields[3]
        variation.alt    = fields[4].split(',')
        variation.qual   = float(fields[5]) if fields[5] != "." else None
        variation.filter = fields[6]

        if len(fields) >= 8:
            # Field INFO
            if fields[7] != '.':
                info = dict()
                for tag_and_value in fields[7].split(';'):
                    if "=" in tag_and_value:
                        tag, value = tag_and_value.split('=')
                        if self.info[tag]["number"] is None:
                            info[tag] = [self.info[tag]["type"](list_elt) for list_elt in value.split(",")]
                        elif self.info[tag]["number"] == 1:
                            info[tag] = self.info[tag]["type"](value)
                        elif field_format["number"] > 1:
                            info[tag] = [self.info[tag]["type"](list_elt) for list_elt in value.split(",")]
                        else: # Number == 0
                            info[tag] = None
                    else:
                        info[tag_and_value] = True
                variation.info = info

            if len(fields) >= 9:
                # Field FORMAT
                variation.format = fields[8].split(':') if fields[8] != "." else None

                # Fields samples
                data_by_spl = dict()
                for spl_idx, spl_cell in enumerate(fields[9:]):
                    spl_data = dict()
                    for field_idx, field_data in enumerate(spl_cell.split(':')):
                        field_id = variation.format[field_idx]
                        field_format = self.format[field_id]
                        if field_format["number"] is None:
                            spl_data[field_id] = [self.format[field_id]["type"](list_elt) for list_elt in field_data.split(",")]
                        elif field_format["number"] == 1:
                            spl_data[field_id] = self.format[field_id]["type"](field_data)
                        elif field_format["number"] > 1:
                            spl_data[field_id] = [self.format[field_id]["type"](list_elt) for list_elt in field_data.split(",")]
                        else: # Number == 0
                            spl_data[field_id] = None
                    data_by_spl[self.samples[spl_idx]] = spl_data
                variation.samples = data_by_spl

        return variation

"""
def _write_header(self):
    self.file_handle.write( "##fileformat=VCFv4.0\n" )
    for tag in self.format:
        self.file_handle.write( '##FORMAT=<' + \
            'ID=' + tag + ',' + \
            'Number=' + str(self.format[tag]["number"]) + ',' + \
            'Type=' + type(self.format[tag]["type"]) + ',' + \ #######################
            'Description="' + self.format[tag]["description"] + "'" + \
            '>'
    self.file_handle.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([spl for spl in self.samples]) + "\n" )

def _parse_header(self):
    pass
"""

def get_VCF_samples( filepath ):
    samples = list()
    with open( filepath ) as FH:
        for line in FH: #################################################### Limit read
            if line.startswith("#CHROM\tPOS"):
                samples = [spl.strip() for spl in line.split("\t")[9:]]
    return samples
