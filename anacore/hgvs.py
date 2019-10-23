# -*- coding: utf-8 -*-
"""Classes and functions for parsing and managing HGVS and mutalyzer results."""


__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re


class Accession:
    """Class to manipulate accession number."""

    def __init__(self, id, version=None, source=None):
        """
        Build and return an instance of Accession.

        :param id: ID of the element (from NCBI accession: NP_001245144.1 => id.version).
        :type id: str
        :param version: Version of the element (from NCBI accession: NP_001245144.1 => id.version).
        :type version: int
        :param source: Source of the element. Example: NCBI or ENSEMBL
        :type source: str
        :return: The new instance.
        :rtype: Accession
        """
        self.id = id
        self.version = int(version) if version is not None else None
        self.source = source

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        repr = self.id
        if self.version is not None:
            repr = "{}.{}".format(self.id, self.version)
        return repr


class HGVS:
    """Class to manipulate HGVS information."""

    def __init__(self, accession, type, change):
        """
        Build and return an instance of HGVS.

        :param accession: Accesion of the alterated element.
        :type accession: Accession
        :param type: HGVS type: "g" or "c" or "p".
        :type type: str
        :param change: Change on element alterated (from HGVS string: "accession:type.change").
        :type change: str
        :return: The new instance.
        :rtype: HGVS
        """
        self.accession = accession
        self.type = type
        self.change = change

    def isPredicted(self):
        """
        Return True if the change come from prediction (opposite to experimetal evidence).

        :return: True if the change come from prediction (opposite to experimetal evidence).
        :rtype: bool
        """
        return "(" in self.change or ")" in self.change

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        return "{}:{}.{}".format(self.accession, self.type, self.change)

    @staticmethod
    def fromStr(hgvs_str):
        """
        Return an instance of HGVS from an HGVS string.

        :param hgvs_str: The HGVS string to parse.
        :type hgvs_str: str
        :return: The instance of HGVS.
        :rtype: HGVS
        """
        hgvs = None
        match = re.search("^(.+)\.(.+):(.)\.(.+)$", hgvs_str)
        if match is not None:
            acc_id, acc_version, hgvs_type, hgvs_change = match.groups()
            acc_obj = Accession(acc_id, acc_version)
            hgvs = HGVS(acc_obj, hgvs_type, hgvs_change)
        else:
            match = re.search("^(.+):(.)\.(.+)$", hgvs_str)
            if match is not None:
                acc_id, hgvs_type, hgvs_change = match.groups()
                acc_obj = Accession(acc_id)
                hgvs = HGVS(acc_obj, hgvs_type, hgvs_change)
            else:
                raise Exception('The HGVS "{}" has an invalid format.'.format(hgvs_str))
        return hgvs


class RunMutalyzerLegend:
    """Class to manipulate field legend from runMutalyzer[Light]."""

    def __init__(self, legend_data):
        """
        Build and return an instance of RunMutalyzerLegend.

        :param legend_data: Legend field from runMutalyzer[Light].
        :type legend_data: list
        :return: The new instance.
        :rtype: RunMutalyzerLegend
        """
        self.data = legend_data

    def getIdByName(self):
        """
        Return element accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.

        :return: element accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.
        :rtype: dict
        """
        id_by_name = {}
        for elt in self.data:
            if "id" in elt and "name" in elt:
                id_by_name[elt["name"]] = elt["id"]
        return id_by_name

    def getProtBytr(self):
        """
        Return protein accession by transcript accession.

        :return: Protein accession by transcript accession.
        :rtype: dict
        """
        prot_by_tr = {}
        id_by_name = self.getIdByName()
        for elt_name, elt_id in id_by_name.items():
            if "_i" in elt_name:
                prot_acc = elt_id
                tr_acc = id_by_name[elt_name.replace("_i", "_v")]
                prot_by_tr[tr_acc] = prot_acc
        return prot_by_tr


class RunMutalyzerDescription:
    """Class to manipulate field proteinDescriptions or transcriptDescriptions from runMutalyzer[Light]."""

    def __init__(self, description_data):
        """
        Build and return an instance of RunMutalyzerLegend.

        :param description_data: proteinDescriptions or transcriptDescriptions field from runMutalyzer[Light].
        :type description_data: list
        :return: The new instance.
        :rtype: RunMutalyzerLegend
        """
        self.data = description_data

    def getByAccession(self, id_by_name):
        """
        Return HGVS by RefSeq accession from proteinDescriptions or transcriptDescriptions of runMutalyzer[Light].

        :param id_by_name: Protein/Transcript RefSeq accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.
        :type id_by_name: dict
        :return: HGVS by RefSeq accession.
        :rtype: dict
        """
        HGVS_by_elt = {}
        for elt in self.data:
            match_renamed = re.search("(.+)\((.+)\):(.)\.(.+)", elt)
            if match_renamed:
                ref_acc, elt_name, hgvs_type, change = match_renamed.groups()
                elt_acc = elt_name
                if elt_name in id_by_name:
                    elt_acc = id_by_name[elt_name]
                if change == "?":
                    HGVS_by_elt[elt_acc] = ""
                else:
                    HGVS_by_elt[elt_acc] = "{}:{}.{}".format(elt_acc, hgvs_type, change)
            else:
                match_std = re.search("(.+):(.)\.(.+)", elt)
                if match_std:
                    elt_acc, hgvs_type, change = match_std.groups()
                    if change == "?":
                        HGVS_by_elt[elt_acc] = ""
                    else:
                        HGVS_by_elt[elt_acc] = "{}:{}.{}".format(elt_acc, hgvs_type, change)
                else:
                    raise Exception("The feature description {} cannot be parsed.".format(elt))
        return HGVS_by_elt
