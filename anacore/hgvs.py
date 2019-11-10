# -*- coding: utf-8 -*-
"""Classes and functions for parsing and managing HGVS and mutalyzer results."""


__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import time
import base64
import requests
import urllib.parse


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


class MutalyzerBatch:
    """
    Manage submission of query batch on a mutalyzer server.

    :Code example:
        .. highlight:: python
        .. code-block:: python

            mutalyzer_url = "https://mutalyzer.nl"
            queries = ["NC_000012.12:g.25245351C>T", "NG_012337.1(SDHD_v001):c.274G>T"]
            batch = MutalyzerBatch(mutalyzer_url, queries, "NameChecker")
            for curr_res in batch.syncRequest():
                print(res)
    """

    def __init__(self, queries, url="https://mutalyzer.nl", command="NameChecker", proxy_url=None):
        """
        Build and return an instance of MutalyzerBatch.

        :param url: URL for the mutalyzer server.
        :type url: str
        :param queries: List of queries (example: ["NC_000012.12:g.25245351C>T", "NG_012337.1(SDHD_v001):c.274G>T"]).
        :type queries: list
        :param command: Type of the batch job, choose from: NameChecker, SyntaxChecker, PositionConverter and SnpConverter.
        :type command: str
        :param proxy_url: URL for the proxy server if necessary.
        :type proxy_url: str
        :return: The new instance.
        :rtype: MutalyzerBatch
        """
        self._url = url
        self._command = command
        if command not in {"NameChecker", "SyntaxChecker", "PositionConverter", "SnpConverter"}:
            raise Exception("{} is not in authorized commands: {}".format(command, ["NameChecker", "SyntaxChecker", "PositionConverter", "SnpConverter"]))
        self._queries = queries
        self._proxy_url = proxy_url
        self._status = None
        self._progress = None
        self._id = None
        self._submit_time = None
        self._complete_time = None

    def getExecTime(self):
        """
        Return server execution time for the batch.

        :return: The server execution time for the batch.
        :rtype: float
        """
        return self._complete_time - self._submit_time

    def getRequestURL(self):
        """
        Return the request URL to submit the batch.

        :return: The request URL to submit the batch.
        :rtype: str
        """
        request_data = base64.urlsafe_b64encode("\n".join(self._queries).encode("utf-8"))
        request_data = urllib.parse.quote(request_data, safe='')
        url = '{}/json/submitBatchJob?process={};data={}'.format(self._url, self._command, request_data)
        return url

    def submit(self):
        """Submit the batch."""
        url_request = self.getRequestURL()
        self._submit_time = time.time()
        response = requests.get(
            url_request,
            proxies=(None if self._proxy_url is None else {"https": self._proxy_url, "http": self.self._proxy_url})
        )
        self._status = "submitted"
        if response.status_code != 200:
            self._status = "error"
            raise Exception("Request {} has failed.".format(url_request))
        self._id = response.json()

    def updateStatus(self):
        """Update progress and status of the submitted batch."""
        url_request = '{}/json/monitorBatchJob?job_id={}'.format(self._url, self._id)
        response = requests.get(
            url_request,
            proxies=(None if self._proxy_url is None else {"https": self._proxy_url, "http": self.self._proxy_url})
        )
        if response.status_code != 200:
            self._status = "error"
            raise Exception("Request {} has failed.".format(url_request))
        nb_left = response.json()
        self._progress = (len(self._queries) - nb_left) / len(self._queries)
        if nb_left == 0:
            self._status = "complete"
            self._complete_time = time.time()

    def getResponse(self):
        """
        Return sever result for the submitted batch.

        :return: The sever result for the submitted batch.
        :rtype: str
        """
        url_request = '{}/json/getBatchJob?job_id={}'.format(self._url, self._id)
        response = requests.get(
            url_request,
            proxies=(None if self._proxy_url is None else {"https": self._proxy_url, "http": self.self._proxy_url})
        )
        if response.status_code != 200:
            self._status = "error"
            raise Exception("Request {} has failed.".format(url_request))
        self._status = "success"
        return str(base64.b64decode(response.json()), "utf-8")

    def getParsedResponse(self, time_to_update=2):
        """
        Return sever result for the submitted batch in list format.

        :param time_to_update: Number of seconds to wait before retry to get response from server. After the status becomes "complete", the server take few seconds to update last results.
        :type time_to_update: int
        :return: The sever result for the submitted batch in list format.
        :rtype: list
        """
        response_lines = self.getResponse().strip().split("\n")
        while len(response_lines) - 1 < len(self._queries):  # After complete the result must be update on server
            time.sleep(time_to_update)
            response_lines = self.getResponse().strip().split("\n")
        titles = response_lines.pop(0).split("\t")
        results = []
        for curr_res in response_lines:
            curr_res = curr_res.split("\t")
            results.append({key: val for key, val in zip(titles, curr_res)})
        return results

    def syncRequest(self, time_to_update=3, first_time_to_update=1):
        """
        Return sever result for the batch in list format.

        :param time_to_update: Number of seconds to wait between each update of status from server.
        :type time_to_update: int
        :param first_time_to_update: Number of seconds to wait before the first status update after submit.
        :type first_time_to_update: int
        :return: The sever result for the batch in list format.
        :rtype: list
        """
        self.submit()
        time.sleep(first_time_to_update)
        self.updateStatus()
        while self._status not in {"error", "complete"}:
            time.sleep(time_to_update)
            self.updateStatus()
        return self.getParsedResponse()


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
