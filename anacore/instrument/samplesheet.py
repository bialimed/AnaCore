# -*- coding: utf-8 -*-
"""Classes and functions required for samplesheet management."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'


def cleanedEnd(value):
    """
    Return line without ending empty CSV fields. This method can be used to removes invalid comma used in CSV to obtain the same number of columns in whole file.

    :param value: Line to clean.
    :type value: str
    :return: Line without ending empty CSV fields.
    :rtype: str
    """
    value = value.strip()
    while value.endswith(","):
        value = value[:-1].strip()
    return value


class Sample:
    """Class to manage a sample from sample sheet."""

    def __init__(self, id, barcodes=None, description=None, metadata=None):
        """
        Build and return an instance of Sample.

        :param id: Sample ID (example: splA).
        :type id: str
        :param barcodes: Sample molecular index used in demultiplex. Example:
            * For single-end: {"index": "barcode01"}
            * For paired-end: {"index": "barcode01_FWD", "index2": "barcode01_REV"}
        :type barcodes: dict
        :param description: Sample description.
        :type description: str
        :param metadata: Additional data on sample (example: {"panel": "Enrich_Solid-Tumor_V3", "project": "Rimo"}). These data come from additional column in section Data.
        :type metadata: dict
        :return: The new instance.
        :rtype: Sample
        """
        self.barcodes = dict() if barcodes is None else barcodes
        self.description = description
        self.id = id
        self.metadata = dict() if metadata is None else metadata

    @property
    def basename(self):
        """
        Return the start of file basename corresponding to the sample (example: splA).

        :return: Sample basename.
        :rtype: str
        """
        return self.id

    def toDict(self):
        """
        Return dict representation of the instance.

        :return: Dict representation of the instance.
        :rtype: dict
        """
        res = self.__dict__
        res["basename"] = self.basename
        return res
