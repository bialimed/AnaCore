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
__version__ = '1.7.2'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'


import re
import gzip
import warnings
from copy import deepcopy


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


class VCFRecord:
    def __init__(self, region=None, position=None, knownSNPId=None, refAllele=None, altAlleles=None, qual=None, pFilter=None, info=None, pFormat=None, samples=None):
        self.chrom   = region
        self.pos     = position
        self.id      = knownSNPId
        self.ref     = refAllele
        self.alt     = altAlleles
        self.qual    = qual
        self.filter  = pFilter if pFilter is not None else list()
        self.info    = info if info is not None else dict()
        self.format  = pFormat if pFormat is not None else list()
        self.samples = samples if samples is not None else dict()

    def containsIndel(self):
        """
        @summary: Returns True if the variant contains an allele corresponding to an insertion or a deletion.
        @return: [bool] True if the variant contains an allele corresponding to an insertion or a deletion.
        @note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        contains_indel = False
        ref = self.ref.replace(".", "").replace("-", "")
        alt = [curr_alt.replace(".", "").replace("-", "") for curr_alt in self.alt]
        for allele in alt:
            if len(allele) != len(ref):
                contains_indel = True
        return contains_indel

    def isDeletion(self):
        """
        @summary: Returns True if the variant is a deletion.
        @return: [bool] True if the variant is a deletion.
        @warnings: This method can only be used on record with only one alternative allele.
        @note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isDeletion' cannot be used on multi-allelic variant.")
        is_deletion = False
        ref = self.ref.replace(".", "").replace("-", "")
        alt = self.alt[0].replace(".", "").replace("-", "")
        if len(alt) < len(ref):
            is_deletion = True
        return is_deletion

    def isIndel(self):
        """
        @summary: Returns True if the variant is an insertion or a deletion.
        @return: [bool] True if the variant is an insertion or a deletion.
        @warnings: This method can only be used on record with only one alternative allele.
        @note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isIndel' cannot be used on multi-allelic variant.")
        is_indel = False
        ref = self.ref.replace(".", "").replace("-", "")
        alt = self.alt[0].replace(".", "").replace("-", "")
        if len(alt) != len(ref):
            is_indel = True
        return is_indel

    def isInsertion(self):
        """
        @summary: Returns True if the variant is an insertion.
        @return: [bool] True if the variant is an insertion.
        @warnings: This method can only be used on record with only one alternative allele.
        @note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isInsertion' cannot be used on multi-allelic variant.")
        is_insertion = False
        ref = self.ref.replace(".", "").replace("-", "")
        alt = self.alt[0].replace(".", "").replace("-", "")
        if len(alt) > len(ref):
            is_insertion = True
        return is_insertion

    def type(self):
        """
        @summary: Returns the vrariant type.
        @return: [str] 'snp' or 'indel' or 'variation'.
        @warnings: This method can only be used on record with only one alternative allele.
        """
        record_type = "snp"
        if self.isIndel():
            record_type = "indel"
        elif len(self.ref) > 1:
            record_type = "variation"
        return record_type

    def standardizeSingleAllele(self, empty_marker="."):
        """
        @summary: The empty allele marker is replaced by the empty_marker and the alternative and reference allele are reduced to the minimal string. The position of record is also updated. Example: ATG/A becomes TG/. ; AAGC/ATAC becomes AG/TA.
        @param empty_marker: The value used for the empty allele (example: A/. for a deletion of A).
        @warnings: This method can only be used on record with only one alternative allele.
        """
        def twoSideTrimming( record ):
            """
            @summary: Removes identical end and start between reference allele and alternative allele (example: ATGACT/ATCT becomes G/). The end is removed first to manage repeat cases.
            @param record: [VCFRecord] The mono-alternative variant.
            """
            while record.ref != "" and record.alt[0] != "" and record.ref[-1] == record.alt[0][-1]:
                record.ref = record.ref[:-1]
                record.alt[0] = record.alt[0][:-1]
            while record.ref != "" and record.alt[0] != "" and record.ref[0] == record.alt[0][0]:
                record.ref = record.ref[1:]
                record.alt[0] = record.alt[0][1:]
                record.pos += 1

        if len(self.alt) > 1:
            raise Exception("The function 'standardizeSingleAllele' cannot be used on multi-allelic variant.")
        self.ref = self.ref.upper()
        self.alt[0] = self.alt[0].upper()
        # Deletion with marker
        if self.alt[0] in [".", "-", ""]:
            self.alt[0] = empty_marker
        # Deletion without marker
        elif len(self.alt[0]) < len(self.ref):
            if self.ref.startswith(self.alt[0]):
                self.ref = self.ref[len(self.alt[0]):]
                self.pos += len(self.alt[0])
                self.alt[0] = empty_marker
            else:
                twoSideTrimming( self )
                if self.alt[0] != "":
                    warnings.warn( 'The deletion "' + self.ref + '/' + self.alt[0] + '" at location ' + self.chrom + ':' + str(self.pos) + ' cannot be standardized.' )
                else:
                    self.alt[0] = empty_marker
        # Insertion with marker
        elif self.ref in [".", "-", ""]:
            self.ref = empty_marker
        # Insertion without marker
        elif len(self.alt[0]) > len(self.ref):
            if self.alt[0].startswith(self.ref):
                self.alt[0] = self.alt[0][len(self.ref):]
                self.pos += len(self.ref)
                self.ref = empty_marker
            else:
                twoSideTrimming( self )
                if self.ref != "":
                    warnings.warn( 'The insertion "' + self.ref + '/' + self.alt[0] + '" at location ' + self.chrom + ':' + str(self.pos) + ' cannot be standardized.' )
                else:
                    self.ref = empty_marker
        # Substitution
        elif len(self.alt[0]) == len(self.ref) and len(self.ref) != 1:
            twoSideTrimming( self )

    def getMostUpstream(self, ref_seq):
        """
        @summary: Returns the most upstream variant that can have the same alternative sequence of the instance.
        @param ref_seq: [str] The reference sequence where the variant has been identified (example: the sequence of the chromosome).
        @return: [VCFRecord] The standardized most upstream variant (see standardizeSingleAllele).
        @warnings: This method can only be used on record with only one alternative allele.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'getMostUpstream' cannot be used on multi-allelic variant.")
        new_record = deepcopy(self)
        new_record.standardizeSingleAllele()
        if new_record.ref == "." or new_record.alt[0] == ".": # Standardized indel
            uc_ref_seq = ref_seq.upper()
            ref = new_record.ref
            alt = new_record.alt[0]
            # Deletion
            if new_record.isDeletion():
                if uc_ref_seq[(new_record.pos - 1):(new_record.pos - 1 + len(ref))] != ref:
                    raise Exception('The reference on position ' + new_record.chrom + ':' + str(new_record.pos) + ' does not correspond to "' + ref + '".')
                before_var = uc_ref_seq[0:(new_record.pos - 1)]
                while before_var != "" and before_var[-1] == ref[-1]:
                    # shift to upstream
                    before_var = before_var[:-1]
                    ref = ref[-1] + ref[:-1]
                new_record.pos = len(before_var) + 1
                new_record.ref = ref
            # Insertion
            else:
                before_var = uc_ref_seq[0:(new_record.pos - 1)]
                while before_var != "" and before_var[-1] == alt[-1]:
                    # shift to upstream
                    before_var = before_var[:-1]
                    alt = alt[-1] + alt[:-1]
                new_record.pos = len(before_var) + 1
                new_record.alt = [alt]
        return new_record

    def getMostDownstream(self, ref_seq):
        """
        @summary: Returns the most downstream variant that can have the same alternative sequence of the instance.
        @param ref_seq: [str] The reference sequence where the variant has been identified (example: the sequence of the chromosome).
        @return: [VCFRecord] The standardized most downstream variant (see standardizeSingleAllele).
        @warnings: This method can only be used on record with only one alternative allele.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'getMostDownstream' cannot be used on multi-allelic variant.")
        new_record = deepcopy(self)
        new_record.standardizeSingleAllele()
        if new_record.ref == "." or new_record.alt[0] == ".": # Standardized indel
            uc_ref_seq = ref_seq.upper()
            ref = new_record.ref
            alt = new_record.alt[0]
            # Deletion
            if new_record.isDeletion():
                if uc_ref_seq[(new_record.pos - 1):(new_record.pos - 1 + len(ref))] != ref:
                    raise Exception('The reference on position ' + new_record.chrom + ':' + str(new_record.pos) + ' does not correspond to "' + ref + '".')
                after_var = uc_ref_seq[new_record.pos + len(ref) - 1:]
                nb_move = 0
                while after_var != "" and after_var[0] == ref[0]:
                    # shift to downstream
                    after_var = after_var[1:]
                    ref = ref[1:] + ref[0]
                    nb_move += 1
                new_record.pos += nb_move
                new_record.ref = ref
            # Insertion
            else:
                after_var = uc_ref_seq[new_record.pos - 1:]
                nb_move = 0
                while after_var != "" and after_var[0] == alt[0]:
                    # shift to upstream
                    after_var = after_var[1:]
                    alt = alt[1:] + alt[0]
                    nb_move += 1
                new_record.pos += nb_move
                new_record.alt = [alt]
        return new_record

    def getPopDP( self ):############################################### test
        """
        @summary: Returns the depth for the population (it is composed by all samples).
        @return: [int] The depth.
        """
        DP = None
        if "DP" in self.info: # The DP is already processed for the population
            DP = self.info["DP"]
        else:
            spl_wout_DP = [spl_name for spl_name in self.samples if "DP" not in self.samples[spl_name]]
            if len(spl_wout_DP) == 0: # All samples have DP
                DP = sum([self.samples[spl_name]["DP"] for spl_name in self.samples])
            else:
                raise Exception( 'The depth cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )
        return( DP )

    def getAD( self, spl_name ):############################################### test
        AD = None

        # Retrieve AD from self
        if "AD" in self.samples[spl_name]: # The AD is already processed for the sample
            AD = self.samples[spl_name]["AD"]
        elif len(self.samples) == 1 and spl_name in self.samplesend and "AD" in self.info: # Only one sample and AD is already processed for population
            AD = self.info["AD"]
        else: # AD must be calculated
            try:
                AF = self.get_AF( spl_name )
                DP = self.get_DP( spl_name )
                AD = [round(curr_AF * DP, 0) for curr_AF in AF]
            except:
                raise Exception( 'The alternative alleles depths cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )

        # Transform AD to list
        if not isinstance(AD, (list, tuple)):
            AD = [AD]

        # Remove the reference allele frequency
        if len(AD) == len(self.alt) + 1:
            AD = AD[1:]

        return( AD )

    def getPopAF( self ):############################################### test
        """
        @summary: Returns the list of alleles frequencies for the population (it is composed by all samples). The reference frequency is removed from the result if it exists.
        @return: [list] The list of alleles frequencies.
        """
        # Retrieve DP
        DP = None
        try:
            DP = self.getPopDP()
        except:
            DP = None

        # Retrieve AF from self
        AF = None
        if "AF" in self.info: # The AF is already processed for the population
            AF = self.info["AF"]
        elif "AD" in self.info and DP is not None: # The AF can be processed directly from the population information
            AD = self.info["AD"]
            AF = [curr_AD/float(DP) for curr_AD in AD]
        else: # The AF must be calculated from samples
            spl_names = list(self.samples.keys())
            if len(self.samples) == 1 and "AF" in self.samples[spl_names[0]]: # Only one sample and it contains AF
                AF = self.samples[spl_names[0]]["AF"]
            else:
                try:
                    AD = None
                    for idx_spl, spl_name in enumerate(self.samples):
                        if idx_spl == 0:
                            AD = self.getAD( spl_name )
                        else:
                            for idx_allele, curr_AD in self.getAD( spl_name ):
                                AD[idx_allele] += curr_A
                    AF = [curr_AD/float(DP) for curr_AD in AD]
                except:
                    raise Exception( 'The allele frequency cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )

        # Transform AF to list
        if not isinstance(AF, (list, tuple)):
            AF = [AF]

        # Remove the reference allele frequency
        if len(AF) == len(self.alt) + 1:
            AF = AF[1:]

        return( AF )

    def get_AF( self, spl_name ):
        """
        @summary: Returns the list of alleles frequencies for the specified sample. The reference frequency is removed from the result if it exists.
        @param spl_name: [str] The sample name.
        @return: [list] The list of alleles frequencies.
        """
        # Retrieve DP
        DP = None
        try:
            DP = self.get_DP( spl_name )
        except:
            DP = None

        # Retrieve AF from self
        AF = None
        if "AF" in self.samples[spl_name]: # The AF is already processed for the sample
            AF = self.samples[spl_name]["AF"]
            if isinstance(AF, (list, tuple)) and len(AF) == len(self.alt) + 1: # Skip the reference allele frequency
                AF = AF[1:]
        elif "AD" in self.samples[spl_name] and DP is not None: # The AF must be processed for the sample
            AD = self.samples[spl_name]["AD"]
            if len(AD) == len(self.alt) + 1: # Skip the reference allele depth
                AD = AD[1:]
            AF = [curr_AD/float(DP) for curr_AD in AD]
        elif len(self.samples) == 1 and spl_name in self.samples and "AF" in self.info: # Only one sample and AF is already processed for population
            AF = self.info["AF"]
            if isinstance(AF, (list, tuple)) and len(AF) == len(self.alt) + 1: # Skip the reference allele frequency
                AF = AF[1:]
        elif len(self.samples) == 1 and spl_name in self.samples and "AD" in self.info and DP is not None: # Only one sample and AF must be processed for population
            AD = self.info["AD"]
            if len(AD) == len(self.alt) + 1: # Skip the reference allele depth
                AD = AD[1:]
            AF = [curr_AD/float(DP) for curr_AD in AD]
        else:
            raise Exception( 'The allele frequency cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )

        # Transform AF to list
        if not isinstance(AF, (list, tuple)):
            AF = [AF]

        return( AF )

    def get_AD( self, spl_name ):
        """
        @summary: Returns the list of alleles depths for the specified sample. The reference depth is removed from the result if it exists.
        @param spl_name: [str] The sample name.
        @return: [list] The list of alleles depths.
        """
        # Retrieve DP
        DP = None
        try:
            DP = self.get_DP( spl_name )
        except:
            DP = None

        # Retrieve AD from self
        AD = None
        if "AD" in self.samples[spl_name]: # The AD is already processed for the sample
            AD = self.samples[spl_name]["AD"]
            if isinstance(AD, (list, tuple)) and len(AD) == len(self.alt) + 1: # Skip the reference allele depth
                AD = AD[1:]
        elif "AF" in self.samples[spl_name] and DP is not None: # The AD must be processed for the sample
            AF = self.samples[spl_name]["AF"]
            if len(AF) == len(self.alt) + 1: # Skip the reference allele frequency
                AF = AF[1:]
            AD = [int(curr_AF * DP) for curr_AF in AF]
        elif len(self.samples) == 1 and spl_name in self.samples and "AD" in self.info: # Only one sample and AD is already processed for population
            AD = self.info["AD"]
            if isinstance(AD, (list, tuple)) and len(AD) == len(self.alt) + 1: # Skip the reference allele depth
                AD = AD[1:]
        elif len(self.samples) == 1 and spl_name in self.samples and "AF" in self.info and DP is not None: # Only one sample and AD must be processed for population
            AF = self.info["AF"]
            if len(AF) == len(self.alt) + 1: # Skip the reference allele frequency
                AF = AF[1:]
            AD = [int(curr_AF * DP) for curr_AF in AF]
        else:
            raise Exception( 'The allele frequency cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )

        # Transform AF to list
        if not isinstance(AD, (list, tuple)):
            AD = [AD]

        return( AD )

    def get_DP( self, spl_name ):
        """
        @summary: Returns the depth for the specified sample.
        @param spl_name: [str] The sample name.
        @return: [int] The depth.
        """
        DP = None
        if "DP" in self.samples[spl_name]: # The DP is already processed for the sample
            DP = self.samples[spl_name]["DP"]
        elif len(self.samples) == 1 and spl_name in self.samples and "DP" in self.info: # Only one sample and DP is already processed for population
            DP = self.info["DP"]
        else:
            raise Exception( 'The depth cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".' )
        return( DP )


class VCFIO:
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
        self.current_line_nb = 0
        self.current_line = None
        self.samples = list()
        self.info = dict() # { "IDREP":{"type": int, "number": 1, "description": "Number of times RU is repeated in indel allele."} }
        self.format = dict() # { "GT":{"type": str, "number": 1, "description": "Genotype"} }
        if mode == "r":
            self._parse_header()

    def __del__( self ):
        self.close()

    def __iter__( self ):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            if self.current_line.startswith('#'):
                self._parse_header_line()
                continue
            try:
                vcf_record = self._parse_line()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
            else:
                yield vcf_record

    def close( self ):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _parse_header( self ):
        """
        @summary: Parses VCF header to set info, format and samples attributes.
        """
        if self.current_line_nb == 0:
            self.current_line_nb += 1
            self.current_line = self.file_handle.readline().rstrip()
        is_header = False
        if self.current_line.startswith('#'):
            is_header = True
        while is_header:
            self._parse_header_line()
            if self.current_line is None or self.current_line.startswith("#CHROM"): # Last header line
                is_header = False
            else:
                self.current_line_nb += 1
                self.current_line = self.file_handle.readline().rstrip()

    def _parse_header_line( self ):
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
            number_tag = match.group(1)
            number = None
            if match.group(1) not in [".", "A", "R", "G"]:
                number = int(match.group(1))
            self.info[id] = { "type": type_fct[type], "type_tag": type, "number": number, "number_tag": number_tag, "description": description }
        elif self.current_line.startswith("##FORMAT"):
            line = self.current_line[10:-1] # Remove "##FORMAT=<" and ">"
            match = re.search('ID=([^,]+)', line)
            id = match.group(1)
            match = re.search('Type=([^,]+)', line)
            type = match.group(1)
            match = re.search('Description="([^"]+)"', line)
            description = match.group(1)
            match = re.search('Number=([^,]+)', line)
            number_tag = match.group(1)
            number = None
            if match.group(1) not in [".", "A", "R", "G"]:
                number = int(match.group(1))
            self.format[id] = { "type": type_fct[type], "type_tag": type, "number": number, "number_tag": number_tag, "description": description }
        elif self.current_line.startswith("#CHROM\tPOS"):
            self.samples = [spl.strip() for spl in self.current_line.split("\t")[9:]]

    def _parse_line( self ):
        """
        @summary: Returns a structured record from the VCF current line.
        @return: [VCFRecord] The variant described by the current line.
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        variation = VCFRecord()
        variation.chrom  = fields[0]
        variation.pos    = int(fields[1])
        variation.id     = fields[2]
        variation.ref    = fields[3]
        variation.alt    = fields[4].split(',')
        variation.qual   = float(fields[5]) if fields[5] != "." else None
        variation.filter = fields[6].split(";") if fields[6] != "." and fields[6] != "" else None

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
                        if field_format["number"] is None or field_format["number"] > 1:
                            spl_data[field_id] = list()
                            for list_elt in field_data.split(","):
                                if list_elt == ".":
                                    spl_data[field_id].append( None )
                                else:
                                    spl_data[field_id].append( self.format[field_id]["type"](list_elt) )
                        elif field_format["number"] == 1:
                            if field_data == ".":
                                spl_data[field_id] = None
                            else:
                                spl_data[field_id] = self.format[field_id]["type"](field_data)
                        else: # Number == 0
                            spl_data[field_id] = None
                    data_by_spl[self.samples[spl_idx]] = spl_data
                variation.samples = data_by_spl

        return variation

    def write(self, record):
        self.file_handle.write( self.recToVCFLine(record) + "\n" )

    def recToVCFLine(self, record):
        """
        @summary: Returns the record in VCF format.
        @param record: [VCFRecord] The record to process.
        @return: [str] The VCF line.
        """
        # Standard columns
        line = "\t".join([
            record.chrom,
            str(record.pos),
            ("." if record.id is None else record.id),
            record.ref,
            ",".join(record.alt),
            ("." if record.qual is None else str(record.qual)),
            ("." if record.filter is None else ";".join(record.filter))])
        # Info
        if record.info is None:
            line += "\t."
        else:
            info_fields = list()
            for key in sorted(record.info):
                if self.info[key]["number"] is None or self.info[key]["number"] > 1:
                    info_fields.append( key + "=" + ",".join(map(str, record.info[key])) )
                else:
                    if self.info[key]["type"] is None: ############################################## Flag
                        info_fields.append( key )
                    else:
                        info_fields.append( key + "=" + str(record.info[key]) )
            line += "\t" + ";".join(info_fields)
        # Format
        if len(self.format) != 0:
            # Format column
            if record.format is None or record.format == ".":
                line += "\t."
            else:
                line += "\t" + ":".join(record.format)
            # Samples columns
            for spl in self.samples:
                if record.samples is None or len(record.samples) == 0:
                    line += "\t" + ":".join(["." for idx in range(len(record.format))])
                else:
                    spl_fields = list()
                    for key in record.format:
                        if key not in record.samples[spl]:
                            spl_fields.append(".")
                        else:
                            if self.format[key]["number"] is None or self.format[key]["number"] > 1:
                                values = list()
                                for current_val in record.samples[spl][key]:
                                    val = (str(current_val) if current_val is not None else ".")
                                    values.append( val )
                                spl_fields.append( ",".join(values) )
                            else:
                                if self.format[key]["type"] is None: ############################################## Flag
                                    spl_fields.append( key )
                                else:
                                    val = (str(record.samples[spl][key]) if record.samples[spl][key] is not None else ".")
                                    spl_fields.append( val )
                    line += "\t" + ":".join(spl_fields)
        return line

    def _writeHeader( self ):
        """
        @note: Draft
        """
        self.file_handle.write( "##fileformat=VCFv4.0\n" )
        for tag in self.info:
            number_tag = self.info[tag]["number"] if self.info[tag]["number"] is not None else "."
            if "number_tag" in self.info[tag]:
                number_tag = self.info[tag]["number_tag"]
            self.file_handle.write( '##INFO=<' + \
                'ID=' + tag + ',' + \
                'Number=' + str(number_tag) + ',' + \
                'Type=' + self.info[tag]["type_tag"] + ',' + \
                'Description="' + self.info[tag]["description"] + '"' + \
                '>\n'
            )
        for tag in self.format:
            number_tag = self.format[tag]["number"] if self.format[tag]["number"] is not None else "."
            if "number_tag" in self.format[tag]:
                number_tag = self.format[tag]["number_tag"]
            self.file_handle.write( '##FORMAT=<' + \
                'ID=' + tag + ',' + \
                'Number=' + str(number_tag) + ',' + \
                'Type=' + self.format[tag]["type_tag"] + ',' + \
                'Description="' + self.format[tag]["description"] + '"' + \
                '>\n'
            )
        self.file_handle.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([spl for spl in self.samples]) + "\n" )

def getAlleleRecord( FH_vcf, record, idx_alt ):
    new_record = VCFRecord(
        region=record.chrom,
        position=record.pos,
        knownSNPId=record.id,
        refAllele=record.ref.upper(), ########################### pb transfo
        altAlleles=[record.alt[idx_alt].upper()], ########################### pb transfo
        qual=record.qual,
        pFilter=deepcopy(record.filter),
        pFormat=deepcopy(record.format)
    )
    # Info
    for key in sorted(record.info):
        if FH_vcf.info[key]["number_tag"] == "A":
            new_record.info[key] = [record.info[key][idx_alt]]
        elif FH_vcf.info[key]["number_tag"] == "R":
            new_record.info[key] = [record.info[key][0], record.info[key][idx_alt + 1]]
        elif FH_vcf.info[key]["number"] is None or FH_vcf.info[key]["number"] > 1:
            new_record.info[key] = [elt for elt in record.info[key]]
        else:
            new_record.info[key] = record.info[key]
    # Samples
    for spl in record.samples:
        new_record.samples[spl] = dict()
        for key in record.format:
            if key in record.samples[spl]:
                if FH_vcf.format[key]["number_tag"] == "A":
                    new_record.samples[spl][key] = [record.samples[spl][key][idx_alt]]
                elif FH_vcf.format[key]["number_tag"] == "R":
                    new_record.samples[spl][key] = [record.samples[spl][key][0], record.samples[spl][key][idx_alt + 1]]
                elif FH_vcf.format[key]["number"] is None or FH_vcf.format[key]["number"] > 1:
                    new_record.samples[spl][key] = [elt for elt in record.samples[spl][key]]
                else:
                    new_record.samples[spl][key] = record.samples[spl][key]
    return new_record
