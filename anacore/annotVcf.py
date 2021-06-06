# -*- coding: utf-8 -*-
"""
Classes and functions for reading/writing/processing annotated VCF.

:Example:

    Read annotated VCF by line

    .. highlight:: python
    .. code-block:: python

        from anacore.annotVcf import AnnotVCFIO

        # File content>
        # ##fileformat=VCFv4.3
        # ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|HGVSp">
        # ##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
        # ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        # #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_test
        # 17\t3675200\t.\tC\tG\t4845.9\tPASS\t\tDP:AD\t4560:1500
        # 17\t7675200\t.\tC\tG\t31799.6\tPASS\tCSQ=G|missense_variant|MODERATE|TP53|NP_000537.3%3Ap.(Ala138Pro),G|missense_variant|MODERATE|TP53|NP_001119584.1%3Ap.(Ala138Pro)\tDP:AD\t4560:1500

        print("Alteration", "Symbol", "HGVSp", sep="\\t")
        with AnnotVCFIO("test.vcf.gz", annot_field="CSQ") as reader:
            for record in reader:
                if "CSQ" not in record.info:
                    print(record.getName(), "", "", sep="\\t")
                else:
                    for annot in record.info["CSQ"]:
                        print(
                            record.getName(),
                            annot["Symbol"],
                            annot["HGVSp"],
                            sep="\\t"
                        )

        # Result>
        # Alteration\tSymbol\tHGVSp
        # 17:3675200=C/G\t\t
        # 17:7675200=C/G\tTP53\tNP_000537.3%3Ap.(Ala138Pro)
        # 17:7675200=C/G\tTP53\tNP_001119584.1:p.(Ala138Pro)

    Write annotated VCF

    .. highlight:: python
    .. code-block:: python

        from anacore.annotVcf import AnnotVCFIO

        annot_field = "ANN"
        with AnnotVCFIO("test.vcf.gz", "w", annot_field=annot_field) as writer:
            # Header
            writer.samples = ["my_sample"]
            self.info = {
                annot_field: HeaderInfoAttr(
                    annot_field,
                    "Consequence annotations. Format: SYMBOL|HGVSp",
                    "String",
                    "."
                )
            }
            self.format = {
                "AF": HeaderFormatAttr("AF", "Allele Frequency", "Float", "A")
            }
            writer.writeHeader()
            # Record
            for record in vcf_record_list:
                record.info[annot_field] = [
                    {"SYMBOL": "TP53", "HGVSp": "NP_000537.3:p.(Ala138Pro)"},
                    {"SYMBOL": "TP53", "HGVSp": "NP_001119584.1%3Ap.(Ala138Pro)"},
                ]
                writer.write(record)

        # Result>
        # ##fileformat=VCFv4.3
        # ##INFO=<ID=ANN,Number=.,Type=String,Description="Allele Frequency">
        # ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
        # #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmy_sample
        # 17\t7675200\t.\tC\tG\t.\tPASS\tANN=TP53|NP_000537.3%3Ap.(Ala138Pro),TP53|NP_001119584.1%3Ap.(Ala138Pro)\tAF\t0.33
        # ...
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.5.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.vcf import *
from copy import deepcopy
import re


class AnnotVCFIO(VCFIO):
    """Manage VCF file containing variants annotations."""

    def __init__(self, filepath, mode="r", annot_field="ANN"):
        """
        Return instance of AnnotVCFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a', 'i'). The mode 'i' allow to open file in indexed mode to fetch by region (tabix).
        :type mode: str
        :param annot_field: The tag for the field used to store annotations. [Default: ANN]
        :type annot_field: str
        :return: The new instance.
        :rtype: AnnotVCFIO
        """
        self.annot_field = annot_field
        self.ANN_titles = list()
        super().__init__(filepath, mode)

    def copyHeader(self, model):
        """
        Copy header fields from the specified VCF.

        :param model: The source.
        :type model: vcf.VCFIO
        """
        super().copyHeader(model)
        if hasattr(model, 'annot_field') and model.annot_field is not None:
            self.annot_field = model.annot_field
            if hasattr(model, 'ANN_titles') and len(model.ANN_titles) != 0:
                self.ANN_titles = model.ANN_titles

    def _parseAnnot(self, record):
        """
        Transform annotations stored as string in record to a structured annotations (list of dict).

        :param record: The record containing the annotations.
        :type record: vcf.VCFRecord
        """
        if self.annot_field in self.info:
            if self.annot_field not in record.info or record.info[self.annot_field] is None:  # ANN is empty
                record.info[self.annot_field] = list()
            else:  # The variant has consequence(s)
                for annot_idx, annot in enumerate(record.info[self.annot_field]):
                    annot_val = {key: None for key in self.ANN_titles}
                    for field_idx, field_value in enumerate(annot.split("|")):
                        csq_key = self.ANN_titles[field_idx]
                        if field_value.strip() != "":
                            annot_val[csq_key] = field_value.strip()
                        else:
                            annot_val[csq_key] = None
                    record.info[self.annot_field][annot_idx] = annot_val

    def _parseHeader(self):
        """Retrieve and store annotations subfields titles and order from the annotation description in INFO header."""
        super()._parseHeader()
        if self.annot_field in self.info:
            # Get ANN fields
            ann_description = self.info[self.annot_field].description
            match = re.search("Format: ([^ ]+)", ann_description)  # VEP
            if match is None:
                match = re.search("Functional annotations: ?'([^']+)'", ann_description)  # SnpEff and jannovar
                if match is None:
                    raise Exception("The {} description cannot be parsed in file {}.".format(self.annot_field, self.filepath))
            titles_str = match.group(1)
            if titles_str.endswith("."):
                titles_str = titles_str[:-1]
            self.ANN_titles = [elt.strip() for elt in titles_str.split("|")]

    def _parseLine(self):
        """
        Return the record corresponding to the current line.

        :return: The record corresponding to the current line.
        :rtype: vcf.VCFRecord
        """
        record = super()._parseLine()
        self._parseAnnot(record)
        return record

    def recToVCFLine(self, record):
        """
        Return the record in VCF format.

        :param record: The record to process.
        :type record: VCFRecord
        :return: The VCF line.
        :rtype: str
        """
        std_record = record
        if self.annot_field in record.info:
            std_record = deepcopy(record)
            if record.info[self.annot_field] is None or len(record.info[self.annot_field]) == 0:  # ANN is empty
                std_record.info.pop(self.annot_field, None)
            else:  # The variant has consequence(s)
                csq_fields = list()
                for annot in record.info[self.annot_field]:
                    annot_fields = list()
                    for csq_key in self.ANN_titles:
                        if annot[csq_key] is None:
                            annot_fields.append("")
                        else:
                            annot_fields.append(str(annot[csq_key]))
                    csq_fields.append("|".join(annot_fields))
                std_record.info[self.annot_field] = csq_fields
        return super().recToVCFLine(std_record)

    def writeHeader(self):
        """Write VCF header."""
        # Manage declaration of ANN in header
        if self.annot_field in self.info:
            match = re.search("(Format: [^ ]+)", self.info[self.annot_field].description)
            old_titles = match.group(1)
            if old_titles.endswith("."):
                old_titles = old_titles[:-1]
            self.info[self.annot_field].description = self.info[self.annot_field].description.replace(
                old_titles,
                "Format: " + "|".join(self.ANN_titles)
            )
        # Write header
        super().writeHeader()


class VEPVCFIO(AnnotVCFIO):
    """Manage VCF file containing variants annotations produced by VEP and stored in CSQ field."""

    def __init__(self, filepath, mode="r"):
        """
        Return instance of VEPVCFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a', 'i'). The mode 'i' allow to open file in indexed mode to fetch by region (tabix).
        :type mode: str
        :return: The new instance.
        :rtype: VEPVCFIO
        """
        super().__init__(filepath, mode, "CSQ")
