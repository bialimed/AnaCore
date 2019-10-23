# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing/processing annotated VCF."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
from copy import deepcopy
from anacore.vcf import *


class AnnotVCFIO(VCFIO):
    """Manage VCF file containing variants annotations."""

    def __init__(self, filepath, mode="r", annot_field="ANN"):
        """
        Return instance of AnnotVCFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
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
        if hasattr(model, 'annot_field'):
            self.annot_field = model.annot_field
        if hasattr(model, 'ANN_titles'):
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
            match = re.search("Format: ([^ ]+)", self.info[self.annot_field].description)
            if match is None:
                raise Exception("The {} description cannot be parsed in file {}.".format(self.annot_field, self.filepath))
            titles_str = match.group(1)
            if titles_str.endswith("."):
                titles_str = titles_str[:-1]
            self.ANN_titles = titles_str.split("|")

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
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: VEPVCFIO
        """
        super().__init__(filepath, mode, "CSQ")
