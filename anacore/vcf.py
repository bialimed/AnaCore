# -*- coding: utf-8 -*-
"""
Classes and functions for reading/writing/processing VCF.

:Example:

    Read VCF by line

    .. highlight:: python
    .. code-block:: python

        from anacore.vcf import VCFIO

        with VCFIO("test.vcf.gz") as reader:
            print("Variant", "\\t".join(reader.samples), sep="\\t")
            for record in reader:
                alt_freq = [str(record.getAF(curr_spl)) for curr_spl in reader.samples]
                print(record.getName(), "\\t".join(alt_freq), sep="\\t")

        # Result>
        # Variant\tN01\tN02\tN03
        # chr1:35-35=A/T\t0.3\t0.4\t0.1
        # chr1:128-128=G/A\t0.1\t0.6\t0.05

    Read VCF by coordinate

    .. highlight:: python
    .. code-block:: python

        from anacore.vcf import VCFIO

        with VCFIO("test.vcf.gz", "i") as reader:
            for record in reader.getSub("chr1", 10, 100):
                print(record.getName())

        # Result>
        # chr1:35-35=A/T

    Write VCF

    .. highlight:: python
    .. code-block:: python

        from anacore.vcf import VCFIO

        with VCFIO("test.vcf.gz", "w") as writer:
            # Header
            writer.samples = ["my_sample"]
            writer.filter = [
                {"q10": HeaderFilterAttr("q10", "Quality below 10")}
            ]
            self.info = {
                "DB": HeaderInfoAttr("DB", "dbSNP membership, build 129", "Flag", 0)
            }
            self.format = {
                "AF": HeaderFormatAttr("AF", "Allele Frequency", "Float", "A")
            }
            writer.extra_header = [
                "##source=myImputationProgramV3.1",
                "##phasing=partial"
            ]
            writer.writeHeader()
            # Record
            for record in vcf_record_list:
                writer.write(record)

        # Result>
        # ##fileformat=VCFv4.3
        # ##source=myImputationProgramV3.1
        # ##phasing=partial
        # ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
        # ##FILTER=<ID=q10,Description="Quality below 10">
        # ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
        # #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmy_sample
        # chr1\t35\t.\tA\tT.\tPASS\tDB\tAF\t0.1
        # ...
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.32.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.abstractFile import AbstractFile
from copy import deepcopy
from pysam import TabixFile
import sys
import warnings


def decodeInfoValue(val):
    """
    Return text where special characters was unescaped. This method is used to read values coming from INFO or FORMAT.

    :param val: A string value of an item of INFO or FORMAT.
    :type val: str
    :return: Text where special characters was unescaped (%3A, %3B, %3D, %2C).
    :rtype: str
    """
    if "%" not in val:
        return val
    return val.replace("%3A", ":").replace("%3B", ";").replace("%3D", "=").replace("%2C", ",")


def encodeInfoValue(val):
    """
    Return text where special characters was escaped. This method is used to write values coming from INFO or FORMAT.

    :param val: A string value of an item of INFO or FORMAT.
    :type val: str
    :return: Text where special characters was escaped (:, ;, =, ,).
    :rtype: str
    """
    return val.replace(":", "%3A").replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")


def getHeaderAttr(header_line):
    """
    Return an instance of HeaderAttr or a child corresponding to the VCF line.

    :param header_line: Declaration line declaration of an VCF attributes from the VCF header.
    :type header_line: str
    :return: An instance of HeaderAttr or a child corresponding to the VCF line.
    :rtype: HeaderAttr or a child
    """
    # Get content
    header_category, header_content = header_line.split("=", 1)  # ##INFO=<ID=AD,Version="1">
    header_category = header_category[2:]  # ##INFO to INFO
    header_content = header_content[1:-1]  # <ID=AD,Version="1"> to ID=AD,Version="1"
    # Get attributes
    attributes = {}
    attributes_cases = {}
    stack = ""
    opened_quote = False
    for curr_char in header_content:
        if curr_char == "," and not opened_quote:
            key, val = stack.split("=", 1)
            attributes[key.lower()] = val
            if key != key.capitalize():  # Keep case different of capitalize
                attributes_cases[key.lower()] = key
            stack = ""
            opened_quote = False
        elif curr_char == '"':
            if stack[-1] == '\\':
                stack = stack[:-1] + '"'  # replace '\"' by '"'
            else:
                if opened_quote:  # The quote is the second
                    opened_quote = False
                else:  # The quote is the first
                    opened_quote = True
        else:
            stack += curr_char
    if stack != "":
        key, val = stack.split("=", 1)
        attributes[key.lower()] = val
        if key != key.capitalize():  # Keep case different of capitalize
            attributes_cases[key.lower()] = key
    # Return
    header_class_name = "Header{}Attr".format(header_category.capitalize())
    header_class = getattr(sys.modules[__name__], header_class_name)
    attr = header_class(**attributes)
    attr.case_by_attr = attributes_cases
    return attr


class HeaderAttr:
    """Class to manage a dict attribute of the VCF header: one tag FORMAT, one tag INFO, one tag FILTER ..."""

    def __init__(self, id, **kwargs):
        """
        Build and return an instance of HeaderAttr.

        :param id: Identifier of the attribute.
        :type id: str
        :return: The new instance.
        :rtype: HeaderAttr
        """
        super().__setattr__("datastore", kwargs)
        self.id = id
        self._required_attr = ["ID"]  # List of required attributes names.
        self._without_quote = {"ID"}  # List of required attributes without quote in str value.
        self.case_by_attr = {}  # Store case sensitive title for extra attributes if it is different of capitalized attribute name.

    def keys(self):
        """
        Return the list of the attribute’s keys.

        :return: list of the attribute’s keys.
        :rtype: list
        """
        return self.datastore.keys()

    def items(self):
        """
        Return the list of the dictionary’s items ((key, value) pairs).

        :return: List of the dictionary’s items.
        :rtype: list
        """
        return self.datastore.items()

    def __delattr__(self, name):
        if name[0] != "_" and name != "case_by_attr":
            del self.datastore[name]
        else:
            super().delattr(name)

    def __setattr__(self, name, value):
        """
        Assign value to the attribute.

        :param name: The attribute name.
        :type name: str
        :param value: The value to be assigned to the attribute.
        :type value: *
        """
        if name[0] != "_" and name != "case_by_attr":
            self.datastore[name] = value
        else:
            super().__setattr__(name, value)

    def __getattr__(self, name):
        value = None
        if name[0] != "_" and name != "case_by_attr":
            value = self.datastore[name]
        else:
            value = super().__getattr__(name)
        return value

    def __str__(self):
        pre_attr = []
        for key in self._required_attr:
            if key in self._without_quote:
                pre_attr.append('{}={}'.format(key, self.datastore[key.lower()]))
            else:
                pre_attr.append('{}="{}"'.format(key, self.datastore[key.lower()].replace('"', '\\"')))
        extra_keys = {elt.lower() for elt in self.keys()} - {elt.lower() for elt in self._required_attr}
        if len(extra_keys) != 0:
            for key in sorted(extra_keys):
                title = key.capitalize()
                if key in self.case_by_attr:
                    title = self.case_by_attr[key]
                pre_attr.append('{}="{}"'.format(title, self.datastore[key].replace('"', '\\"')))
        return "<{}>".format(",".join(pre_attr))

    def __repr__(self):
        return self.datastore.__repr__()


class HeaderDescAttr(HeaderAttr):
    """Class to manage VCF header attribute where ID and Description are required."""

    def __init__(self, id, description, **kwargs):
        """
        Build and return an instance of HeaderDescAttr.

        :param id: Identifier of the attribute.
        :type id: str
        :param description: Description of the attribute.
        :type description: str
        :return: The new instance.
        :rtype: HeaderDescAttr
        """
        super().__init__(id, **kwargs)
        self.description = description
        self._required_attr = ["ID", "Description"]

    def __setattr__(self, name, value):
        """
        Assign value to the attribute.

        :param name: The attribute name.
        :type name: str
        :param value: The value to be assigned to the attribute.
        :type value: *
        """
        super().__setattr__(name, value)


class HeaderTypedAttr(HeaderDescAttr):
    """Class to manage VCF header attribute where ID, Description, Type and Number are required."""

    def __init__(self, id, description, type="String", number="1", **kwargs):
        """
        Build and return an instance of HeaderTypedAttr.

        :param id: Identifier of the attribute.
        :type id: str
        :param description: Description of the attribute.
        :type description: str
        :param type: Data types. Authorized values: Integer (32-bit, signed), Float (32-bit IEEE-754),Flag, Character, or String.
        :type type: str
        :param number: Number of values that can be included with the field.
        :type number: char
        :return: The new instance.
        :rtype: HeaderTypedAttr
        """
        super().__init__(id, description, **kwargs)
        self.type = type
        self.number = number
        self._required_attr = ["ID", "Number", "Type", "Description"]  # List of required attributes names.
        self._without_quote = {"ID", "Number", "Type"}  # List of required attributes without quote in str value.

    def __setattr__(self, name, value):
        """
        Assign value to the attribute.

        :param name: The attribute name.
        :type name: str
        :param value: The value to be assigned to the attribute.
        :type value: *
        """
        if name == "type":
            type_fct = {
                "String": str,
                "Integer": int,
                "Float": float,
                "Character": str,
                "Flag": None
            }
            self._type = type_fct[value]
        elif name == "number":
            self._number = None
            if value not in [".", "A", "R", "G"]:
                self._number = int(value)
        super().__setattr__(name, value)


class HeaderFilterAttr(HeaderDescAttr):
    """Class to manage a FILTER attribute of VCF header."""

    pass


class HeaderFormatAttr(HeaderTypedAttr):
    """Class to manage a FORMAT attribute of VCF header."""

    pass


class HeaderInfoAttr(HeaderTypedAttr):
    """Class to manage an INFO attribute of VCF header."""

    pass


class HeaderSampleAttr(HeaderAttr):
    """Class to manage an SAMPLE attribute of VCF header."""

    pass


class VCFRecord:
    """Class to manage a variant record."""

    def __init__(
          self, region=None, position=None, knownSNPId=None, refAllele=None,
          altAlleles=None, qual=None, pFilter=None, info=None, pFormat=None,
          samples=None):
        """
        Build and return an instance of VCFRecord.

        :param region: Identifier of the reference region/contig/chromosome containing the variant (example: chr3).
        :type region: str
        :param position: Start position for the variant (1-based).
        :type position: int
        :param knownSNPId: Semi-colon separated list of unique identifiers where available (example rs numbers from dbSNP).
        :type knownSNPId: str
        :param refAllele: Reference allele.
        :type refAllele: char
        :param altAlleles: Alternatives alleles.
        :type altAlleles: list
        :param qual: Variant quality.
        :type qual: int
        :param pFilter: List of filter tags. It can be in three possible states: (1) if no filter was applied, the field contains an empty list ("." in VCF file); (2) if filters were applied but the record passes filters, the field should contain ["PASS"]; (3) if filters were applied and the record does not pass filters, the field should contain ["filter_name", ...].
        :type pFilter: list
        :param info: Variant additionnal information.
        :type info: dict
        :param pFormat: List of format tags.
        :type pFormat: list
        :param samples: Variant format information by sample.
        :type samples: dict
        :return: The new instance.
        :rtype: VCFRecord
        """
        self.chrom = region
        self.pos = position
        self.id = knownSNPId
        self.ref = refAllele
        self.alt = altAlleles
        self.qual = qual
        self.filter = list() if pFilter is None else pFilter
        self.info = dict() if info is None else info
        self.format = list() if pFormat is None else pFormat
        self.samples = dict() if samples is None else samples
        self._normalized = None

    @staticmethod
    def getEmptyAlleleMarker():
        """
        Return the marker internally used to represent an empty allele (reference allele in insertion and alternative allele in deletion).

        :return: The marker.
        :rtype: str
        """
        return "-"

    def __setattr__(self, name, value):
        """
        Assign value to the attribute.

        :param name: The attribute name.
        :type name: str
        :param value: The value to be assigned to the attribute.
        :type value: *
        """
        if value is not None:
            if name == "ref":
                value = (VCFRecord.getEmptyAlleleMarker() if value == "" else value)
                super(VCFRecord, self).__setattr__("_normalized", None)
            elif name == "alt":
                value = [(VCFRecord.getEmptyAlleleMarker() if elt == "" else elt) for elt in value]
                super(VCFRecord, self).__setattr__("_normalized", None)
            elif name == "pos" or name == "chrom":
                super(VCFRecord, self).__setattr__("_normalized", None)
        super(VCFRecord, self).__setattr__(name, value)

    def containsIndel(self):
        """
        Return True if the variant contains an allele corresponding to an insertion or a deletion.

        :return: True if the variant contains an allele corresponding to an insertion or a deletion.
        :rtype: bool
        :note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        contains_indel = False
        ref = self.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = [curr_alt.replace(VCFRecord.getEmptyAlleleMarker(), "") for curr_alt in self.alt]
        for allele in alt:
            if len(allele) != len(ref):
                contains_indel = True
        return contains_indel

    def refStart(self):
        """
        Return the first position on reference affected by the alternative allele.

        :return: The first position on reference affected by the alternative allele. For an insertion between two nucleotids the value will be: first nucleotids pos + 0.5.
        :rtype: float
        :warnings: This method can only be used on record with only one alternative allele.
        :examples:
                # Insertion:
                chr1  12  A    TG    => returns 12
                chr1  12  A    AGT   => chr1  13  .   GT  (std) => returns 12.5
                chr1  11  AA   AGT   => chr1  12  A   GT  (std) => returns 12
                chr1  10  AAA  AAGT  => chr1  12  A   GT  (std) => returns 12
                chr1  10  AATC AACGT => chr1  12  TC  CGT (std) => returns 12

                # Substitution:
                chr1  10  A   T  => returns 10
                chr1  10  AA  TT => returns 10

                # Deletion:
                chr1  12  A     .  => returns 12
                chr1  10  AAA   .  => returns 10
                chr1  10  AAA   A  => chr1  11  AA  . (std) => returns 11
                chr1  10  AGC   A  => chr1  11  GC  . (std) => returns 11
                chr1  10  AAAT  TG => returns 10
        """
        if len(self.alt) > 1:
            raise Exception("The function 'refStart' cannot be used on multi-allelic variant {}.".format(self.getName()))
        if self._normalized is None:
            self._normalized = deepcopy(self)
            self._normalized.normalizeSingleAllele()
        record = self._normalized
        start = record.pos
        if record.ref == VCFRecord.getEmptyAlleleMarker():
            start -= 0.5
        return start

    def refEnd(self):
        """
        Return the last position on reference affected by the alternative allele.

        :return: The last position on reference affected by the alternative allele. For an insertion between two nucleotids the value will be: first nucleotids pos + 0.5.
        :rtype: float
        :warnings: This method can only be used on record with only one alternative allele.
        :examples:
                # Insertion:
                chr1  12  A    TG    => returns 12
                chr1  12  A    AGT   => chr1  13  -   GT  (std) => returns 12.5
                chr1  11  AA   AGT   => chr1  12  A   GT  (std) => returns 12
                chr1  10  AAA  AAGT  => chr1  12  A   GT  (std) => returns 12
                chr1  10  AATC AACGT => chr1  12  TC  CGT (std) => returns 13

                # Substitution:
                chr1  10  A   T  => returns 10
                chr1  10  AA  TT => returns 11

                # Deletion:
                chr1  12  A     -  => returns 12
                chr1  10  AAA   -  => returns 12
                chr1  10  AAA   A  => chr1  11  AA  - (std) => returns 12
                chr1  10  AGC   A  => chr1  11  GC  - (std) => returns 12
                chr1  10  AAAT  TG => returns 13
        """
        if len(self.alt) > 1:
            raise Exception("The function 'refEnd' cannot be used on multi-allelic variant {}.".format(self.getName()))
        if self._normalized is None:
            self._normalized = deepcopy(self)
            self._normalized.normalizeSingleAllele()
        record = self._normalized
        end = record.pos
        if record.ref == VCFRecord.getEmptyAlleleMarker():
            end -= 0.5
        else:
            end += len(record.ref) - 1
        return end

    def getName(self):
        """
        Return an unique name to identified the variant.

        :return: The variant name.
        :rtype: str
        """
        return "{}:{}={}/{}".format(
            self.chrom,
            self.pos,
            self.ref,
            "/".join(self.alt)
        )

    def isDeletion(self):
        """
        Return True if the variant is a deletion.

        :return: True if the variant is a deletion.
        :rtype: bool
        :warnings: This method can only be used on record with only one alternative allele.
        :note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isDeletion' cannot be used on multi-allelic variant {}.".format(self.getName()))
        is_deletion = False
        ref = self.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = self.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
        if len(alt) < len(ref):
            is_deletion = True
        return is_deletion

    def isInsAndDel(self):
        """
        Return True if the variant is an insertion and also a deletion.

        :return: True if the variant is an insertion and also a deletion.
        :rtype: bool
        :warnings: This method can only be used on record with only one alternative allele.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isInsAndDel' cannot be used on multi-allelic variant {}.".format(self.getName()))
        is_ins_and_del = False
        if self._normalized is None:
            self._normalized = deepcopy(self)
            self._normalized.normalizeSingleAllele()
        ref = self._normalized.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = self._normalized.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
        if len(ref) > len(alt):  # Deletion exists and is more longer (eg: AT/ or AT/G)
            if len(alt) != 0:  # Insertion exists (eg: AT/G)
                is_ins_and_del = True
        elif len(ref) < len(alt):  # Insertion exists and is more longer (eg: /AT or G/AT)
            if len(ref) != 0:  # Deletion exists (eg: G/AT)
                is_ins_and_del = True
        else:  # Same length: SNV or MNV
            if len(ref) != 1:  # MNV: AT/GC
                is_ins_and_del = True
        return is_ins_and_del

    def isIndel(self):
        """
        Return True if the variant is an insertion or a deletion.

        :return: True if the variant is an insertion or a deletion.
        :rtype: bool
        :warnings: This method can only be used on record with only one alternative allele.
        :note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isIndel' cannot be used on multi-allelic variant {}.".format(self.getName()))
        is_indel = False
        ref = self.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = self.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
        if len(alt) != len(ref):
            is_indel = True
        return is_indel

    def isInsertion(self):
        """
        Return True if the variant is an insertion.

        :return: True if the variant is an insertion.
        :rtype: bool
        :warnings: This method can only be used on record with only one alternative allele.
        :note: If the alternative allele and the reference alle have the same length the variant is considered as a substitution. For example: AA/GC is considered as double substitution and not as double insertion after double deletion.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'isInsertion' cannot be used on multi-allelic variant {}.".format(self.getName()))
        is_insertion = False
        ref = self.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = self.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
        if len(alt) > len(ref):
            is_insertion = True
        return is_insertion

    def type(self):
        """
        Return the variant type.

        :return: 'snp' or 'indel' or 'variation'.
        :rtype: str
        :warnings: This method can only be used on record with only one alternative allele.
        """
        record_type = "snp"
        if self.isIndel():
            record_type = "indel"
        elif len(self.ref) > 1:
            record_type = "variation"
        return record_type

    def fastDownstreamed(self, seq_handler, padding=500):
        """
        Return a simplified record (CHROM, POS, ALT and REF) moved to the most downstream postition.

        :param seq_handler: File handle to the reference sequences file.
        :type seq_handler: anacore.sequenceIO.IdxFastaIO
        :param padding: Number of nucleotids to inspect after variant. Upstream movement is limited to this number of nucleotids.
        :type padding: int
        :return: The simplified record moved to the most downstream postition.
        :rtype: anacore.vcf.VCFRecord
        """
        if len(self.alt) > 1:
            raise Exception('The function "fastDownstreamed" cannot be used on multi-allelic variant {}.'.format(self.getName()))
        downstream_rec = VCFRecord(self.chrom, self.pos, None, self.ref, [elt for elt in self.alt])
        if self._normalized is not None:
            downstream_rec._normalized = deepcopy(self._normalized)
        ref = self.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
        alt = self.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
        if len(ref) != 1 or len(alt) != 1:
            # Move to downstream
            sub_region = seq_handler.getSub(self.chrom, self.pos, self.pos + len(self.ref) + padding)
            downstream_rec.pos = 1  # Switch position from chromosome to position from subregion
            downstream_rec = downstream_rec.getMostDownstream(sub_region)
            downstream_rec.pos = self.pos + downstream_rec.pos - 1  # Switch position from subregion to position from chromosome
            if len(ref) != len(alt):
                # Add previous nt
                prev_nt = seq_handler.getSub(downstream_rec.chrom, downstream_rec.pos - 1, downstream_rec.pos - 1)
                if downstream_rec.ref == downstream_rec.getEmptyAlleleMarker():  # Insertion
                    downstream_rec.ref = prev_nt
                    downstream_rec.alt[0] = prev_nt + downstream_rec.alt[0]
                    downstream_rec.pos -= 1
                elif downstream_rec.alt[0] == VCFRecord.getEmptyAlleleMarker():  # Deletion
                    downstream_rec.ref = prev_nt + downstream_rec.ref
                    downstream_rec.alt[0] = prev_nt
                    downstream_rec.pos -= 1
        return downstream_rec

    def fastStandardize(self, seq_handler, padding=500):
        """
        Standardize record (move to the most upstream position and remove unecessary nucleotids). This standardization concerns only POS, REF and ALT.

        :param seq_handler: File handle to the reference sequences file.
        :type seq_handler: anacore.sequenceIO.IdxFastaIO
        :param padding: Number of nucleotids to inspect before variant. Upstream movement is limited to this number of nucleotids.
        :type padding: int
        """
        # Move to upstream
        sub_start = max(self.pos - padding, 1)
        real_padding = self.pos - sub_start
        sub_region = seq_handler.getSub(self.chrom, sub_start, self.pos + len(self.ref))
        self.pos = real_padding + 1  # Switch position from chromosome to position from subregion
        upstream_rec = self.getMostUpstream(sub_region)
        self.pos = upstream_rec.pos - 1 + sub_start  # Switch position from subregion to position from chromosome
        self.ref = upstream_rec.ref
        self.alt[0] = upstream_rec.alt[0]
        # Add previous nt
        if self.isIndel():
            prev_nt = seq_handler.getSub(self.chrom, self.pos - 1, self.pos - 1)
            if self.ref == VCFRecord.getEmptyAlleleMarker():  # Insertion
                self.ref = prev_nt
                self.alt[0] = prev_nt + self.alt[0]
                self.pos -= 1
            elif self.alt[0] == VCFRecord.getEmptyAlleleMarker():  # Deletion
                self.ref = prev_nt + self.ref
                self.alt[0] = prev_nt
                self.pos -= 1

    def normalizeSingleAllele(self):
        """
        Replace empty allele by the empty_marker and reduce the alternative and reference allele to the minimal string. The position of record is also updated. Example: ATG/A becomes TG/. ; AAGC/ATAC becomes AG/TA.

        :warnings: This method can only be used on record with only one alternative allele.
        """
        def twoSideTrimming(record):
            """
            Remove identical end and start between reference allele and alternative allele (example: ATGACT/ATCT becomes G/). The end is removed first to manage repeat cases.

            :param record: The mono-alternative variant.
            :type record: VCFRecord
            """
            ref = record.ref
            alt = record.alt[0]
            while ref != "" and alt != "" and ref[-1] == alt[-1]:
                ref = ref[:-1]
                alt = alt[:-1]
            while ref != "" and alt != "" and ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                record.pos += 1
            record.ref = ref
            record.alt = [alt]

        if len(self.alt) > 1:
            raise Exception("The function 'normalizeSingleAllele' cannot be used on multi-allelic variant {}.".format(self.getName()))
        self.ref = self.ref.upper()
        self.alt[0] = self.alt[0].upper()
        # Deletion or insertion with marker
        if self.alt[0] == VCFRecord.getEmptyAlleleMarker() or self.ref == VCFRecord.getEmptyAlleleMarker():
            pass
        # Deletion without marker
        elif len(self.alt[0]) < len(self.ref):
            if self.ref.startswith(self.alt[0]):
                self.ref = self.ref[len(self.alt[0]):]
                self.pos += len(self.alt[0])
                self.alt[0] = VCFRecord.getEmptyAlleleMarker()
            else:
                twoSideTrimming(self)
                if self.alt[0] != VCFRecord.getEmptyAlleleMarker():
                    warnings.warn(
                        'The deletion "{}/{}" at location {}:{} cannot be normalized.'.format(
                            self.ref, self.alt[0], self.chrom, self.pos
                        )
                    )
        # Insertion without marker
        elif len(self.alt[0]) > len(self.ref):
            if self.alt[0].startswith(self.ref):
                self.alt[0] = self.alt[0][len(self.ref):]
                self.pos += len(self.ref)
                self.ref = VCFRecord.getEmptyAlleleMarker()
            else:
                twoSideTrimming(self)
                if self.ref != VCFRecord.getEmptyAlleleMarker():
                    warnings.warn(
                        'The insertion "{}/{}" at location {}:{} cannot be normalized.'.format(
                            self.ref, self.alt[0], self.chrom, self.pos
                        )
                    )
        # Substitution
        elif len(self.alt[0]) == len(self.ref) and len(self.ref) != 1:
            twoSideTrimming(self)

    def getMostUpstream(self, ref_seq):
        """
        Return the most upstream variant that can have the same alternative sequence of the instance.

        :param ref_seq: The reference sequence where the variant has been identified (example: the sequence of the chromosome).
        :type ref_seq: str
        :return: The normalized most upstream variant (see normalizeSingleAllele).
        :rtype: VCFRecord
        :warnings: This method can only be used on record with only one alternative allele.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'getMostUpstream' cannot be used on multi-allelic variant {}.".format(self.getName()))
        if self._normalized is None:
            self._normalized = deepcopy(self)
            self._normalized.normalizeSingleAllele()
        new_record = self._normalized
        if new_record.ref == VCFRecord.getEmptyAlleleMarker() or new_record.alt[0] == VCFRecord.getEmptyAlleleMarker():  # normalized indel
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
        Return the most downstream variant that can have the same alternative sequence of the instance.

        :param ref_seq: The reference sequence where the variant has been identified (example: the sequence of the chromosome).
        :type ref_seq: str
        :return: The normalized most downstream variant (see normalizeSingleAllele).
        :rtype: VCFRecord
        :warnings: This method can only be used on record with only one alternative allele.
        """
        if len(self.alt) > 1:
            raise Exception("The function 'getMostDownstream' cannot be used on multi-allelic variant {}.".format(self.getName()))
        if self._normalized is None:
            self._normalized = deepcopy(self)
            self._normalized.normalizeSingleAllele()
        new_record = self._normalized
        if new_record.ref == VCFRecord.getEmptyAlleleMarker() or new_record.alt[0] == VCFRecord.getEmptyAlleleMarker():  # Normalized indel
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

    def getPopRefAD(self):
        """
        Return the reference allele depth for the population (it is composed by all samples).

        :return: The reference allele population depth.
        :rtype: int
        """
        ref_pop_AD = None
        if "AD" in self.info and len(self.info["AD"]) == len(self.alt) + 1:  # INFO.AD contains ref allele
            ref_pop_AD = self.info["AD"]
        elif "AF" in self.info and "DP" in self.info and len(self.info["AF"]) == len(self.alt) + 1:  # INFO.AF contains ref allele and INFO.DP exists
            ref_pop_AD = int(self.info["AF"][0] * self.info["DP"])
        elif len(self.samples) != 0:  # Must be processed by samples data
            first_spl = self.samples.keys()[0]
            if "AD" in self.format and len(first_spl["AD"]) == len(self.alt) + 1:  # AD contains ref allele
                try:
                    ref_pop_AD = 0
                    for spl_name, spl_data in self.samples.items():  # Sum ref AD for all samples
                        ref_pop_AD += spl_data["AD"][0]
                except Exception:
                    ref_pop_AD = None
            if ref_pop_AD is None and "AF" in self.format and "DP" in self.format and len(first_spl["AF"]) == len(self.alt) + 1:  # AF contains ref allele
                try:
                    ref_pop_AD = 0
                    for spl_name, spl_data in self.samples.items():  # Sum ref AD for all samples
                        ref_pop_AD += spl_data["AF"][0] * spl_data["DP"]
                except Exception:
                    ref_pop_AD = None
        if ref_pop_AD is None:  # Erroneous when several variants exist on position and are line splitted
            if "AD" in self.info:  # INFO.AD does not contain ref
                ref_pop_AD = self.getPopDP() - sum(self.info["AD"])
            elif "AF" in self.info and "DP" in self.info:  # INFO.AF does not contain ref and INFO.DP exists
                ref_pop_AD = int((1 - sum(self.info["AF"])) * self.info["DP"])
            else:  # Must be processed by samples data
                alt_pop_AD = sum(self.getPopAltAD())
                ref_pop_AD = self.getPopDP() - alt_pop_AD
        return ref_pop_AD

    def getPopAltAD(self):
        """
        Return the list of alternative alleles depths for the population (it is composed by all samples). The reference depth is removed from the result if it exists.

        :return: The list of alternative alleles depths.
        :rtype: list
        """
        # Retrieve AD from self
        AD = None
        if "AD" in self.info:  # The AD is already processed for the population
            AD = self.info["AD"]
        else:
            # Get population DP
            DP = None
            try:
                DP = self.getPopDP()
            except Exception:
                pass
            if "AF" in self.info and DP is not None:  # The AD can be processed directly from the population information
                AF = self.info["AF"] if isinstance(self.info["AF"], (list, tuple)) else [self.info["AF"]]
                AD = [int(round(curr_AF * DP, 0)) for curr_AF in AF]
            else:  # The AD must be calculated from samples information
                spl_names = list(self.samples.keys())
                if len(self.samples) == 1 and "AD" in self.samples[spl_names[0]]:  # Only one sample and it contains AD
                    AD = self.samples[spl_names[0]]["AD"]
                else:
                    try:
                        for idx_spl, spl_name in enumerate(self.samples):
                            if idx_spl == 0:
                                AD = [curr_AD for curr_AD in self.getAltAD(spl_name)]
                            else:
                                for idx_allele, curr_AD in enumerate(self.getAltAD(spl_name)):
                                    AD[idx_allele] += curr_AD
                    except Exception:
                        raise Exception('The allele depth cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".')
        # Transform AD to list
        if not isinstance(AD, (list, tuple)):
            AD = [AD]
        # Remove the reference allele depth
        if len(AD) == len(self.alt) + 1:
            AD = AD[1:]
        # Return
        return(AD)

    def getPopRefAF(self):
        """
        Return the reference allele frequency for the population (it is composed by all samples).

        :return: The reference allele population frequency.
        :rtype: int
        """
        ref_pop_AF = None
        if "AF" in self.info and len(self.info["AF"]) == len(self.alt) + 1:  # INFO.AF contains ref allele
            ref_pop_AF = self.info["AF"]
        elif "AD" in self.info and "DP" in self.info and len(self.info["AD"]) == len(self.alt) + 1:  # INFO.AD contains ref allele and INFO.DP exists
            ref_pop_AF = self.info["AD"][0] / self.info["DP"]
        elif len(self.samples) != 0:  # Must be processed by samples data
            first_spl = self.samples.keys()[0]
            if "AD" in self.format and "DP" in self.format and len(first_spl["AD"]) == len(self.alt) + 1:  # AD contains ref allele
                try:
                    ref_pop_AD = 0
                    pop_DP = 0
                    for spl_name, spl_data in self.samples.items():  # Sum ref AD for all samples
                        pop_DP += spl_data["DP"]
                        ref_pop_AD += spl_data["AD"][0]
                    ref_pop_AF = ref_pop_AD / pop_DP
                except Exception:
                    ref_pop_AF = None
            if ref_pop_AF is None and "AF" in self.format and "DP" in self.format and len(first_spl["AF"]) == len(self.alt) + 1:  # AF contains ref allele
                try:
                    ref_pop_AD = 0
                    pop_DP = 0
                    for spl_name, spl_data in self.samples.items():  # Sum ref AF for all samples
                        pop_DP += spl_data["DP"]
                        ref_pop_AD += spl_data["AF"][0] * spl_data["DP"]
                    ref_pop_AF = ref_pop_AD / pop_DP
                except Exception:
                    ref_pop_AF = None
        if ref_pop_AF is None:  # Erroneous when several variants exist on position and are line splitted
            if "AF" in self.info:  # INFO.AF does not contain ref
                ref_pop_AF = 1 - sum(self.info["AF"])
            elif "AD" in self.info and "DP" in self.info:  # INFO.AD does not contain ref and INFO.DP exists
                ref_pop_AF = (self.info["DP"] - sum(self.info["AD"])) / self.info["DP"]
            else:  # Must be processed by samples data
                alt_pop_AF = sum(self.getPopAltAF())
                ref_pop_AF = 1 - alt_pop_AF
        return ref_pop_AF

    def getPopAltAF(self):
        """
        Return the list of alleles frequencies for the population (it is composed by all samples). The reference frequency is removed from the result if it exists.

        :return: The list of alleles frequencies.
        :rtype: list
        """
        # Retrieve AF from self
        AF = None
        if "AF" in self.info:  # The AF is already processed for the population
            AF = self.info["AF"]
        else:
            # Get population DP
            DP = None
            try:
                DP = self.getPopDP()
            except Exception:
                pass
            if "AD" in self.info and DP is not None:  # The AF can be processed directly from the population information
                AD = self.info["AD"] if isinstance(self.info["AD"], (list, tuple)) else [self.info["AD"]]
                AF = [curr_AD / DP for curr_AD in AD]
            else:  # The AF must be calculated from samples information
                spl_names = list(self.samples.keys())
                if len(self.samples) == 1 and "AF" in self.samples[spl_names[0]]:  # Only one sample and it contains AF
                    AF = self.samples[spl_names[0]]["AF"]
                else:
                    try:
                        pop_AD = None
                        for idx_spl, spl_name in enumerate(self.samples):
                            if idx_spl == 0:
                                pop_AD = [curr_AD for curr_AD in self.getAltAD(spl_name)]
                            else:
                                for idx_allele, curr_AD in enumerate(self.getAltAD(spl_name)):
                                    pop_AD[idx_allele] += curr_AD
                        if DP == 0:
                            AF = [0 for curr_AD in pop_AD]
                            if sum(pop_AD) != 0:
                                raise Exception('popAD and popDP are not compatible for variant "{}".'.format(self.getName()))
                        else:
                            AF = [curr_AD / DP for curr_AD in pop_AD]
                    except Exception:
                        raise Exception('The allele frequency cannot be retrieved in variant "{}".'.format(self.getName()))
        # Transform AF to list
        if not isinstance(AF, (list, tuple)):
            AF = [AF]
        # Remove the reference allele frequency
        if len(AF) == len(self.alt) + 1:
            AF = AF[1:]
        # Return
        return AF

    def getPopDP(self):
        """
        Return the depth for the population (it is composed by all samples).

        :return: The depth in population.
        :rtype: int
        """
        DP = None
        if "DP" in self.info:  # The DP is already processed for the population
            DP = self.info["DP"]
        else:
            # Get population AD
            AD = None
            if "AD" in self.info:
                AD = self.info["AD"] if isinstance(self.info["AD"], (list, tuple)) else [self.info["AD"]]
            # Calculate DP
            if AD is not None and len(AD) == len(self.alt) + 1:  # The DP can be processed from INFO's AD (it contains the depth for all alleles and reference)
                DP = sum(AD)
            elif AD is not None and "AF" in self.info:  # The DP can be processed from INFO's AD and AF
                AF = self.info["AF"] if isinstance(self.info["AF"], (list, tuple)) else [self.info["AF"]]
                if len(AF) == len(self.alt) + 1:
                    AF = AF[1:]
                DP = int(round(AD[0] / AF[0], 0))
            elif len(self.samples) != 0:  # The DP must be calculated from samples information
                DP = sum([self.getDP(spl_name) for spl_name in self.samples])
            else:
                raise Exception('The population depth cannot be retrieved in variant "{}".'.format(self.getName()))
        return DP

    def getAD(self, spl_name):
        """
        Return the list of alleles depths for the specified sample.

        :param spl_name: The sample name.
        :type spl_name: str
        :return: The list of alleles depths.
        :rtype: list
        """
        AD = None
        # Retrieve AD from self
        if "AD" in self.samples[spl_name]:  # The AD is already processed for the sample
            AD = self.samples[spl_name]["AD"]
        elif len(self.samples) == 1 and spl_name in self.samples and "AD" in self.info:  # Only one sample and AD is already processed for population
            AD = self.info["AD"]
        else:  # AD must be calculated
            try:
                AF = self.getAF(spl_name)
                DP = self.getDP(spl_name)
                AD = [int(round(curr_AF * DP, 0)) for curr_AF in AF]
            except Exception:
                raise Exception('The alternative alleles depths cannot be retrieved in variant "{}".'.format(self.getName()))
        # Transform AD to list
        if not isinstance(AD, (list, tuple)):
            AD = [AD]
        # Return
        return AD

    def getAltAD(self, spl_name):
        """
        Return the list of alternative alleles depths for the specified sample. The reference depth is removed from the result if it exists.

        :param spl_name: The sample name.
        :type spl_name: str
        :return: The list of alternative alleles depths.
        :rtype: list
        """
        AD = self.getAD(spl_name)
        # Remove the reference allele depth
        if len(AD) == len(self.alt) + 1:
            AD = AD[1:]
        # Return
        return AD

    def getAFBySample(self, missing_replacement=None):
        """
        Return the list of alleles frequencies by sample name. The reference frequency is removed from the result if it exists.

        :param missing_replacement: Value used to replace missing AF.
        :type missing_replacement: float
        :return: The list of alleles frequencies by sample name.
        :rtype: dict
        """
        AF = dict()
        for curr_spl in self.samples:
            AF[curr_spl] = list()
            for curr_AF in self.getAF(curr_spl):
                if curr_AF is None and missing_replacement is not None:
                    curr_AF = float(missing_replacement)
                AF[curr_spl].append(curr_AF)
        return AF

    def getAF(self, spl_name):
        """
        Return the list of alleles frequencies for the specified sample.

        :param spl_name: The sample name.
        :type spl_name: str
        :return: The list of alleles frequencies.
        :rtype: list
        """
        # Retrieve AF from self
        AF = None
        if "AF" in self.samples[spl_name]:  # The AF is already processed for the sample
            AF = self.samples[spl_name]["AF"]
        else:
            # Get sample AD
            AD = None
            if "AD" in self.samples[spl_name]:
                AD = self.samples[spl_name]["AD"] if isinstance(self.samples[spl_name]["AD"], (list, tuple)) else [self.samples[spl_name]["AD"]]
            if AD is not None and len(AD) == len(self.alt) + 1:  # The AF can be processed from sample's AD (it contains the depth for alleles and reference)
                DP = sum(AD)
                if DP == 0:
                    AF = [0 for curr_AD in AD]
                else:
                    AF = [curr_AD / DP for curr_AD in AD]
            else:
                # Get sample DP
                DP = None
                try:
                    DP = self.getDP(spl_name)
                except Exception:
                    pass
                if AD is not None and DP is not None:  # The AF can be processed from sample's AD and DP
                    if DP == 0:
                        AF = [0 for curr_AD in AD]
                    else:
                        AF = [curr_AD / DP for curr_AD in AD]
                elif len(self.samples) == 1 and spl_name in self.samples and "AF" in self.info:  # Only one sample and AF is already processed for population
                    AF = self.info["AF"]
                # elif len(self.samples) == 1 and spl_name in self.samples and "AD" in self.info and DP is not None: # Only one sample and AF must be processed for population
                #     AD = self.info["AD"]
                #     AF = [curr_AD/float(DP) for curr_AD in AD]
                else:
                    raise Exception('The allele frequency cannot be retrieved in variant "{}".'.format(self.getName()))
        # Transform AF to list
        if not isinstance(AF, (list, tuple)):
            AF = [AF]
        # Return
        return AF

    def getAltAF(self, spl_name):
        """
        Return the list of alternative alleles frequencies for the specified sample. The reference frequency is removed from the result if it exists.

        :param spl_name: The sample name.
        :type spl_name: str
        :return: The list of alternative alleles frequencies.
        :rtype: list
        """
        AF = self.getAF(spl_name)
        # Remove the reference allele frequency
        if len(AF) == len(self.alt) + 1:
            AF = AF[1:]
        # Return
        return AF

    def getDP(self, spl_name):
        """
        Return the depth for the specified sample.

        :param spl_name: The sample name.
        :type spl_name: str
        :return: The depth.
        :rtype: int
        """
        DP = None
        if "DP" in self.samples[spl_name]:  # The DP is already processed for the sample
            DP = self.samples[spl_name]["DP"]
        elif len(self.samples) == 1 and spl_name in self.samples and "DP" in self.info:  # Only one sample and DP is already processed for population
            DP = self.info["DP"]
        elif "AD" in self.samples[spl_name]:  # DP can be calculated
            AD = self.samples[spl_name]["AD"] if isinstance(self.samples[spl_name]["AD"], (list, tuple)) else [self.samples[spl_name]["AD"]]
            if len(AD) == len(self.alt) + 1:  # Sample contains AD for all alleles and reference
                DP = sum(AD)
            elif "AF" in self.samples[spl_name]:  # Sample contains AD for all alleles and AF
                AF = self.samples[spl_name]["AF"] if isinstance(self.samples[spl_name]["AF"], (list, tuple)) else [self.samples[spl_name]["AF"]]
                if len(AF) == len(self.alt) + 1:
                    AF = AF[1:]
                DP = int(round(AD[0] / AF[0], 0))
            else:  # AD does not contain reference AD and AF is missing to calculate DP
                raise Exception('The depth cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".')
        else:  # AD is missing to calculate DP
            raise Exception('The depth cannot be retrieved in variant "' + self.chrom + ":" + str(self.pos) + '".')
        return DP


class VCFIO(AbstractFile):
    """Manage VCF file."""

    def __init__(self, filepath, mode="r"):
        """
        Return instance of VCFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a', 'i'). The mode 'i' allow to open file in indexed mode to fetch by region (tabix).
        :type mode: str
        :return: The new instance
        :rtype: anacore.vcf.VCFIO
        """
        standard_mode = mode if mode != "i" else "r"
        AbstractFile.__init__(self, filepath, standard_mode)
        self.samples = list()
        self.filter = dict()  # {"q10": <HeaderFilterAttr>}
        self.info = dict()  # {"IDREP": <HeaderInfoAttr>, "AD": <HeaderInfoAttr>}
        self.format = dict()  # {"GT": <HeaderFormatAttr>, "AD": <HeaderFormatAttr>}
        self.sample_info = dict()  # {"splA": <HeaderSampleAttr>, "splB": <HeaderSampleAttr>}
        self.extra_header = list()  # ["##source=myImputationProgramV3.1", "##phasing=partial"]
        self._index = None
        if mode == "a":  # Load header in append mode
            with VCFIO(filepath, "r") as tmp_reader:
                self.copyHeader(tmp_reader)
        if mode in {"r", "i"}:
            self._parseHeader()
            if mode == "i":
                self._index = TabixFile(filepath)

    def getSub(self, chr, start, end):
        """
        Return generator on records overlapping the specified region.

        .. warning::
            This method can only be used with an instance of anacore.VCFIO opened
            with mode "i".

        :param chr: Chromosome name of the selected region.
        :type chr: str
        :param start: Start of the selected region (1-based).
        :type start: int
        :param end: End of the selected region (1-based).
        :type end: int
        :return: Records overlappind the specified region.
        :rtype: generator for anacore.vcf.VCFRecord
        """
        if self._index is None:
            raise Exception("The file {} must be open in mode 'i' to use {}.getSub() on it.".format(self.file_handle, __class__.__name__))
        iterator = None
        try:
            iterator = self._index.fetch(chr, start - 1, end)
        except ValueError as error:
            missing_region = "could not create iterator for region '{}:{}-{}'".format(chr, start, end)
            if missing_region == str(error):
                return []
            else:
                raise error
        for row in iterator:
            self.current_line = row
            yield self._parseLine()

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

    def _parseHeader(self):
        """Parse VCF header to set info, format and samples attributes."""
        if self.current_line_nb == 0:
            self.current_line = self.file_handle.readline().rstrip()
        is_header = False
        if self.current_line.startswith('#'):
            is_header = True
        while is_header:
            self._parseHeaderLine()
            if self.current_line is None or self.current_line.startswith("#CHROM"):  # Last header line
                is_header = False
            else:
                self.current_line_nb += 1
                self.current_line = self.file_handle.readline().rstrip()

    def _parseHeaderLine(self):
        """Parse one VCF header line to update info, filter, format or samples attributes."""
        if self.current_line.startswith("##INFO"):
            curr_attr = getHeaderAttr(self.current_line)
            self.info[curr_attr.id] = curr_attr
        elif self.current_line.startswith("##FILTER"):
            curr_attr = getHeaderAttr(self.current_line)
            self.filter[curr_attr.id] = curr_attr
        elif self.current_line.startswith("##FORMAT"):
            curr_attr = getHeaderAttr(self.current_line)
            self.format[curr_attr.id] = curr_attr
        elif self.current_line.startswith("##SAMPLE"):
            curr_attr = getHeaderAttr(self.current_line)
            self.sample_info[curr_attr.id] = curr_attr
        elif self.current_line.startswith("#CHROM\tPOS"):
            self.samples = [spl.strip() for spl in self.current_line.split("\t")[9:]]
        elif not self.current_line.startswith("##fileformat="):
            self.extra_header.append(self.current_line)

    def _parseLine(self):
        """
        Return a structured record from the VCF current line.

        :return: The variant described by the current line.
        :rtype: VCFRecord
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        variation = VCFRecord()
        variation.chrom = fields[0]
        variation.pos = int(fields[1])
        variation.id = fields[2]
        variation.ref = fields[3]
        variation.alt = fields[4].split(',')
        variation.qual = float(fields[5]) if fields[5] != "." else None
        if fields[6] != "." and fields[6] != "":
            variation.filter = fields[6].split(";")

        if len(fields) >= 8:
            # Field INFO
            if fields[7] != '.':
                info = dict()
                for tag_and_value in fields[7].split(';'):
                    if "=" not in tag_and_value:  # The field is a flag (self.info[tag]._number == 0)
                        info[tag_and_value] = True
                    else:
                        tag, value = tag_and_value.split('=', 1)
                        if self.info[tag]._number == 1:  # The field contains an unique value
                            if value != ".":  # Exclude key None value
                                info[tag] = self.info[tag]._type(value)
                                if self.info[tag].type == "String":
                                    info[tag] = decodeInfoValue(info[tag])
                        else:  # The field contains a list (self.info[tag]._number is None or self.info[tag]._number > 1)
                            if value == "":
                                info[tag] = []
                            elif value == ".":
                                pass  # Exclude key with None list
                            elif self.info[tag].type == "String":
                                info[tag] = [decodeInfoValue(self.info[tag]._type(list_elt)) for list_elt in value.split(",")]
                            else:
                                info[tag] = [self.info[tag]._type(list_elt) for list_elt in value.split(",")]
                variation.info = info

            if len(fields) >= 9:
                # Field FORMAT
                variation.format = fields[8].split(':') if fields[8] != "." else None
                # Fields samples
                data_by_spl = dict()
                for spl_idx, spl_cell in enumerate(fields[9:]):
                    spl_data = dict()
                    if variation.format is not None:  # Samples cannot have any data if format is None
                        for field_idx, field_data in enumerate(spl_cell.split(':')):
                            field_id = variation.format[field_idx]
                            field_format = self.format[field_id]
                            if field_format._number is None or field_format._number > 1:  # Value is list
                                if field_data == ".":  # Value is None in VCF
                                    if field_format.number in ["R", "A"]:
                                        spl_data[field_id] = [None for alt in variation.alt]
                                        if field_format.number == "R":
                                            spl_data[field_id].append(None)
                                    elif field_format._number is not None:
                                        spl_data[field_id] = [None for elt in range(field_format._number)]
                                    else:
                                        spl_data[field_id] = [None]
                                else:  # Value is not None in VCF
                                    spl_data[field_id] = list()
                                    for list_elt in field_data.split(","):
                                        if list_elt == ".":
                                            spl_data[field_id].append(None)
                                        else:
                                            value = self.format[field_id]._type(list_elt)
                                            if self.format[field_id].type == "String":
                                                value = decodeInfoValue(value)
                                            spl_data[field_id].append(value)
                            elif field_format._number == 1:  # Value is not a list
                                if field_data == ".":
                                    spl_data[field_id] = None
                                else:
                                    spl_data[field_id] = self.format[field_id]._type(field_data)
                                    if self.format[field_id].type == "String":
                                        spl_data[field_id] = decodeInfoValue(spl_data[field_id])
                            else:  # Number == 0
                                spl_data[field_id] = True
                    data_by_spl[self.samples[spl_idx]] = spl_data
                variation.samples = data_by_spl

        return variation

    def write(self, record):
        """
        Write variant record in VCF.

        :param record: The variant record.
        :type record: VCFRecord
        """
        self.file_handle.write(self.recToVCFLine(record) + "\n")

    def recToVCFLine(self, record):
        """
        Return the record in VCF format.

        :param record: The record to process.
        :type record: VCFRecord
        :return: The VCF line.
        :rtype: str
        """
        # Standard columns
        line = "\t".join([
            record.chrom,
            str(record.pos),
            ("." if record.id is None else record.id),
            record.ref,
            ",".join(record.alt),
            ("." if record.qual is None else str(record.qual)),
            ("." if len(record.filter) == 0 else ";".join(record.filter))
        ])
        # Info
        if record.info is None or len(record.info) == 0 or len(self.info) == 0:
            line += "\t."
        else:
            info_fields = list()
            for key in sorted(record.info):
                if self.info[key]._number is None or self.info[key]._number > 1:  # The info may cointain a list of values
                    values = [encodeInfoValue(str(elt)) for elt in record.info[key]]
                    info_fields.append(key + "=" + ",".join(values))
                else:  # The info contains a flag or a uniq value
                    if self.info[key]._type is None:  # Flag
                        info_fields.append(key)
                    else:
                        value = encodeInfoValue(str(record.info[key]))
                        info_fields.append(key + "=" + value)
            line += "\t" + ";".join(info_fields)
        # Format
        if len(self.format) != 0:  # the VCF contains a column format
            if record.format is None or record.format == ".":  # Current record does not contain information on samples
                # Format column
                line += "\t."
                # Samples columns
                for spl in self.samples:
                    line += "\t."
            else:  # Current record contains information on samples
                # Format column
                line += "\t" + ":".join(record.format)
                # Samples columns
                if record.samples is None or len(record.samples) == 0:
                    for spl in self.samples:
                        line += "\t" + ":".join(["." for elt in record.format])
                else:
                    for spl_name in self.samples:
                        record_spl = record.samples[spl_name]
                        spl_fields = list()
                        for key in record.format:
                            if key not in record_spl:
                                spl_fields.append(".")
                            else:
                                if self.format[key]._number is None or self.format[key]._number > 1:  # The info may cointain a list of values
                                    values = list()
                                    for current_val in record_spl[key]:
                                        value = (encodeInfoValue(str(current_val)) if current_val is not None else ".")
                                        values.append(value)
                                    spl_fields.append(",".join(values))
                                else:  # The format contains a uniq value
                                    value = (encodeInfoValue(str(record_spl[key])) if record_spl[key] is not None else ".")
                                    spl_fields.append(value)
                        line += "\t" + ":".join(spl_fields)
        return line

    def copyHeader(self, model):
        """
        Copy header fields from the specified VCF.

        :param model: The VCF source.
        :type model: VCFIO
        """
        self.filter = deepcopy(model.filter)
        self.format = deepcopy(model.format)
        self.info = deepcopy(model.info)
        self.sample_info = deepcopy(model.sample_info)
        self.samples = deepcopy(model.samples)
        self.extra_header = deepcopy(model.extra_header)

    def writeHeader(self):
        """Write VCF header."""
        self.file_handle.write("##fileformat=VCFv4.3\n")
        for curr in self.extra_header:
            self.file_handle.write(curr + "\n")
        for tag in sorted(self.info):
            self.file_handle.write('##INFO={}\n'.format(self.info[tag]))
        for tag in sorted(self.filter):
            self.file_handle.write('##FILTER={}\n'.format(self.filter[tag]))
        for tag in sorted(self.format):
            self.file_handle.write('##FORMAT={}\n'.format(self.format[tag]))
        for spl in sorted(self.sample_info):
            self.file_handle.write('##SAMPLE={}\n'.format(self.sample_info[spl]))
        last_header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        if len(self.format) != 0 or len(self.samples) != 0:
            last_header_line += "\tFORMAT\t" + "\t".join([spl for spl in self.samples])
        self.file_handle.write(last_header_line + "\n")


def getAlleleRecord(FH_vcf, record, idx_alt):
    """
    Return the record corresponding to the specified allele in variant.

    :param FH_vcf: Handler on record source file.
    :type FH_vcf: anacore.vcf.VCFIO
    :param record: The variant record.
    :type record: anacore.vcf.VCFRecord
    :param idx_alt: The index of the allele in alt attribute.
    :type idx_alt: int
    :return: The record corresponding to the specified allele in variant.
    :rtype: anacore.vcf.VCFRecord
    """
    new_record = VCFRecord(
        region=record.chrom,
        position=record.pos,
        knownSNPId=record.id,
        refAllele=record.ref.upper(),  ########################### pb transfo
        altAlleles=[record.alt[idx_alt].upper()],  ########################### pb transfo
        qual=record.qual,
        pFilter=deepcopy(record.filter),
        pFormat=deepcopy(record.format)
    )
    # Info
    for key in sorted(record.info):
        if FH_vcf.info[key].number == "A":
            new_record.info[key] = [record.info[key][idx_alt]]
        elif FH_vcf.info[key].number == "R":
            new_record.info[key] = [record.info[key][0], record.info[key][idx_alt + 1]]
        elif FH_vcf.info[key]._number is None or FH_vcf.info[key]._number > 1:
            new_record.info[key] = [elt for elt in record.info[key]]
        else:
            new_record.info[key] = record.info[key]
    # Samples
    for spl in record.samples:
        new_record.samples[spl] = dict()
        for key in record.format:
            if key in record.samples[spl]:
                if FH_vcf.format[key].number == "A":
                    new_record.samples[spl][key] = [record.samples[spl][key][idx_alt]]
                elif FH_vcf.format[key].number == "R":
                    new_record.samples[spl][key] = [record.samples[spl][key][0], record.samples[spl][key][idx_alt + 1]]
                elif FH_vcf.format[key]._number is None or FH_vcf.format[key]._number > 1:
                    new_record.samples[spl][key] = [elt for elt in record.samples[spl][key]]
                else:
                    new_record.samples[spl][key] = record.samples[spl][key]
    return new_record


def getFreqMatrix(vcf_path, missing_replacement=0.0, accept_missing=True):
    """
    Return a 2D matrix representing variants frequencies by sample.

    :param vcf_path: The pass to the VCF processed.
    :type vcf_path: str
    :param missing_replacement: The value used to replace missing alleles frequencies (example: np.nan).
    :type missing_replacement: *
    :param accept_missing: If false and a missing AF is found an exception is raised. If true the missing AF are replaced by missing_replacement.
    :type accept_missing: bool
    :return: The list of samples, the list of variants names and the 2D matrix of alleles frequencies.
    :rtype: list
    """
    samples = list()
    variants = list()
    AF_matrix = list()
    with VCFIO(vcf_path) as FH_vcf:
        samples = FH_vcf.samples
        for record in FH_vcf:
            for alt_idx, curr_alt in enumerate(record.alt):  # For each alternative allele in variant
                record_allele = getAlleleRecord(FH_vcf, record, alt_idx)
                variants.append(record_allele.getName())
                row_array = list()
                for curr_spl in samples:
                    AF = record_allele.getAF(curr_spl)[0]
                    if AF is None:
                        AF = missing_replacement
                        if not accept_missing:
                            raise Exception("The AF is missing for variant {} in sample {}.".format(curr_spl, record_allele.getName()))
                    row_array.append(AF)
                AF_matrix.append(row_array)
    return samples, variants, AF_matrix
