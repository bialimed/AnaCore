# -*- coding: utf-8 -*-
"""Classes and functions for parsing and managing HGVS and mutalyzer results."""


__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sequence import AA3LettersAlphabet
import base64
from collections import deque
import re
import requests
import time
import urllib.parse


AA_ONE_BY_THREE = AA3LettersAlphabet.one_by_three

AA_THREE_BY_ONE = {value: key for key, value in AA_ONE_BY_THREE.items()}

ONE_LETTER_AA_LEXIC = set(AA_ONE_BY_THREE.values())

THREE_LETTERS_AA_LEXIC = AA3LettersAlphabet.words | {"*", "X"}


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
        match = re.search(r"^(.+)\.(.+):(.)\.(.+)$", hgvs_str)
        if match is not None:
            acc_id, acc_version, hgvs_type, hgvs_change = match.groups()
            acc_obj = Accession(acc_id, acc_version)
            hgvs = HGVS(acc_obj, hgvs_type, hgvs_change)
        else:
            match = re.search(r"^(.+):(.)\.(.+)$", hgvs_str)
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
            match_renamed = re.search(r"(.+)\((.+)\):(.)\.(.+)", elt)
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
                match_std = re.search(r"(.+):(.)\.(.+)", elt)
                if match_std:
                    elt_acc, hgvs_type, change = match_std.groups()
                    if change == "?":
                        HGVS_by_elt[elt_acc] = ""
                    else:
                        HGVS_by_elt[elt_acc] = "{}:{}.{}".format(elt_acc, hgvs_type, change)
                else:
                    raise Exception("The feature description {} cannot be parsed.".format(elt))
        return HGVS_by_elt


class HGVSProtChange:
    """Class to manipulate change part of proteic HGVS (ex: "Val600Glu")."""

    def __init__(self, start_aa, start_pos, end_aa, end_pos, new_aa, evt, new_elts, predicted=True):
        """
        Build and return an instance of HGVSProtChange.

        :param start_aa: First amino acid impacted by the variant. Example: "Val" in "Val600Glu" or "Gln" in "Gln746_Lys747ins*63" or "Ile" in "Ile327Argfs*?".
        :type start_aa: str
        :param start_pos: Proteic position of the first amino acid impacted by the variant. Example: 600 in "Val600Glu" or 746 in "Gln746_Lys747ins*63" or 327 in "Ile327Argfs*?".
        :type start_pos: int
        :param end_aa: Last amino acid impacted by the variant. Example: None in "Val600Glu" or "Lys" in "Gln746_Lys747ins*63" or None in "Ile327Argfs*?".
        :type end_aa: str
        :param end_pos: Proteic position of the last amino acid impacted by the variant. Example: None in "Val600Glu" or 747 in "Gln746_Lys747ins*63" or None in "Ile327Argfs*?".
        :type end_pos: int
        :param new_aa: Amino acid changed by "fs" or "ext". Example: None in "Val600Glu" or None in "Gln746_Lys747ins*63" or "Arg" in "Ile327Argfs*?".
        :type new_aa: str
        :param evt: Type of evenment in [None, "del", "dup", "ext", "fs", "ins", "insdel"]. Example: None in "Val600Glu" or "ins" in "Gln746_Lys747ins*63" or "fs" in "Ile327Argfs*?".
        :type evt: str
        :param new_elts: New elements list after evenment type in HGVS change. Example ["Glu"] from "Val600Glu" or ["*", "63"] from "Gln746_Lys747ins*63" or ["*", "?"] from "Ile327Argfs*?".
        :type new_elts: list
        :param predicted: True if HGVS consequence is predicted, i.e. without experimental evidence (no RNA or protein sequence analysed).
        :type predicted: bool
        :return: The new instance.
        :rtype: anacore.hgvs.HGVSProtChange
        """
        self.start_aa = start_aa
        self.start_pos = start_pos
        self.end_aa = end_aa
        self.end_pos = end_pos
        self.new_aa = new_aa
        self.evt = evt
        self.new_elts = [] if new_elts is None else new_elts
        self.predicted = predicted

    @property
    def new_aa(self):
        return self._new_aa

    @new_aa.setter
    def new_aa(self, new_val):
        self._new_aa = None if new_val is None else HGVSProtChange._cleanedAA(new_val)

    @property
    def start_aa(self):
        return self._start_aa

    @start_aa.setter
    def start_aa(self, new_val):
        self._start_aa = None if new_val is None else HGVSProtChange._cleanedAA(new_val)

    @property
    def end_aa(self):
        return self._end_aa

    @end_aa.setter
    def end_aa(self, new_val):
        self._end_aa = None if new_val is None else HGVSProtChange._cleanedAA(new_val)

    def __eq__(self, other):
        """
        Return true if the object represent the same HGVS change.

        :param other: HGVS change to compare.
        :type other: anacore.hgvs.HGVSProtChange
        :return: True if the object represent the same HGVS change.
        :rtype: bool
        """
        return repr(self) == repr(other)

    def __repr__(self):
        """
        Return the object representation in string.

        :return: The object representation in string.
        :rtype: str
        """
        rep_str = '"start_aa": {}, "start_pos": {}, ' \
            '"end_aa": {}, "end_pos": {}, ' \
            '"new_aa": {}, "evt": {}, ' \
            '"new_elts": {}, "predicted": {}'.format(
                self.start_aa,
                self.start_pos,
                self.end_aa,
                self.end_pos,
                self.new_aa,
                self.evt,
                self.new_elts,
                self.predicted
            )
        return '{' + rep_str + '}'

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        if self.start_aa is None:
            if self.start_pos == 0:
                repr = "0"
            else:
                repr = "?"
        else:
            repr = "{}{}".format(self.start_aa, self.start_pos)
            if self.end_pos:
                repr += "_{}{}".format(self.end_aa, self.end_pos)
            if self.new_aa:
                repr += self.new_aa
            if self.evt:
                repr += self.evt
            repr += "".join(self.new_elts)
        if self.predicted:
            repr = "({})".format(repr)
        return repr

    def isInFrameIns(self):
        """
        Return True if the HGVS is insertion, duplication or repeat. Frameshift or extension resulting from insertion return False.

        :return: Return True if the HGVS is insertion, duplication or repeat.
        :rtype: bool
        """
        is_ins = False
        if self.evt in {"ins", "dup"}:
            is_ins = True
        elif self.evt is None and len(self.new_elts) != 0 and self.new_elts[0] == "[":  # Reapeat
            is_ins = True
        return is_ins

    @staticmethod
    def _cleanedAA(aa):
        """
        Validate and return an amino acid representation for HGVS.

        :type aa: Amino acid in one or three letter representation.
        :type aa: str
        :return: Return three letter representation.
        :rtype: str
        """
        aa = aa.capitalize()
        if len(aa) == 1:
            if aa not in ONE_LETTER_AA_LEXIC:
                raise Exception("{} is not a valid one letter amino acid. Authorized: {}.".change(aa, sorted(ONE_LETTER_AA_LEXIC)))
            aa = AA_THREE_BY_ONE[aa]
        elif len(aa) == 3:
            if aa not in THREE_LETTERS_AA_LEXIC:
                raise Exception("{} is not a valid three letter amino acid. Authorized: {}.".change(aa, sorted(THREE_LETTERS_AA_LEXIC)))
        else:
            raise Exception('The amino acid "{}" it is not in one or three letter representation.'.format(aa))
        return aa

    @staticmethod
    def _getAANbLetter(start_aa, follow):
        """
        Return the number of letters used to encode amino acids.

        :param start_aa: First amino acid impacted by the variant.
        :type start_aa: str
        :param follow: Part after first impacted positions in HGVS change (ex: "Glu" in Val600Glu or "_Lys747ins*63" in Gln746_Lys747ins*63 or "Argfs*?" in Ile327Argfs*?).
        :type follow: str
        :return: Number of letters used to encode amino acids.
        :rtype: int
        """
        aa_nb_letter = len(start_aa)
        if start_aa == "*":
            aa_nb_letter = None
            if follow.startswith("_"):
                idx = 1
                while not follow[idx].isdigit():
                    idx += 1
                aa_nb_letter = idx - 1
            elif len(follow) == 1 and follow not in {"?", "="}:
                aa_nb_letter = 1
            elif len(follow) >= 3:
                if follow[0:3] in THREE_LETTERS_AA_LEXIC:
                    aa_nb_letter = 3
                elif follow[0].isupper():
                    if follow[1:3] == "fs" and follow[0] in ONE_LETTER_AA_LEXIC:
                        aa_nb_letter = 1
                    elif len(follow) >= 4 and follow[1:4] == "ext" and follow[0] in ONE_LETTER_AA_LEXIC:
                        aa_nb_letter = 1
        return aa_nb_letter

    @staticmethod
    def _splittedOnEvt(change, aa_nb_letter, start_pos, end_pos, follow):
        """
        Return (new_aa, evenement, new_elts_string) from part after modified positions in HGVS change (ex: "Glu" in Val600Glu or "ins*63" in Gln746_Lys747ins*63 or "Argfs*?" in Ile327Argfs*?).

        :param change: Change part of HGVS string (ex: "Val600Glu").
        :type change: str
        :param aa_nb_letter: Number of letters to represent an amino acid: 1 or 3.
        :type aa_nb_letter: int
        :param start_pos: First amino acid changed by variant.
        :type start_pos: int
        :param end_pos: Last amino acid changed by variant.
        :type end_pos: int
        :param follow: Part after modified positions in HGVS change (ex: "Glu" in Val600Glu or "ins*63" in Gln746_Lys747ins*63 or "Argfs*?" in Ile327Argfs*?).
        :type follow: str
        :return: The tuple (new_aa, evenement, new_elts_string) from part after modified positions in HGVS change. Examples: (None, None, "Glu") from Glu or (None, "ins", "*63") from ins*63 or ("Arg", "fs", "*?") from Argfs*?.
        :rtype: (str, str, str)
        """
        curr_aa_lexic = THREE_LETTERS_AA_LEXIC if aa_nb_letter == 3 else ONE_LETTER_AA_LEXIC
        re_aa = "|".join(["({})".format(elt.lower()) for elt in curr_aa_lexic])
        re_aa = re_aa.replace("*", r"\*")
        new_aa = None
        evt = None
        lc_follow = follow.lower()
        if lc_follow.startswith("delins"):
            evt = "delins"
            follow = follow[6:]
            if len(follow) == 0:
                raise Exception('"delins" must be followed by the amino acids inserted: {}'.format(change))
        elif lc_follow.startswith("del"):
            evt = "del"
            follow = follow[3:]
            if len(follow) != 0:
                # Follow must be the deleted aa
                if len(follow) % aa_nb_letter != 0:
                    raise Exception('Deletion must be ended by "del" or the amino acids deleted: {}'.format(change))
                del_len = (start_pos if end_pos is None else end_pos) - start_pos + 1
                if len(follow) / aa_nb_letter != del_len:
                    raise Exception('Deletion must be ended by "del" or the amino acids deleted: {}'.format(change))
                del_aa = [follow[idx:idx+aa_nb_letter].capitalize() for idx in range(0, len(follow), aa_nb_letter)]
                for curr_del_aa in del_aa:
                    if curr_del_aa not in curr_aa_lexic:
                        raise Exception('Deletion must be ended by "del" or the amino acids deleted: {}'.format(change))
                follow = ""
        elif lc_follow.startswith("ins"):
            evt = "ins"
            follow = follow[3:]
            if len(follow) == 0:
                raise Exception("Insertion must be followed by the number or the amino acids inserted: {}".format(change))
            if end_pos - start_pos != 1:
                raise Exception("In insertion the end position must be the next amino acid after the start: {}".format(change))
        elif lc_follow.startswith("dup"):
            evt = "dup"
            follow = follow[3:]
            if len(follow) != 0:
                # Follow must be the deleted aa
                if len(follow) % aa_nb_letter != 0:
                    raise Exception('Duplication must be ended by "dup" or the amino acids duplicated: {}'.format(change))
                del_len = (start_pos if end_pos is None else end_pos) - start_pos + 1
                if len(follow) / aa_nb_letter != del_len:
                    raise Exception('Duplication must be ended by "dup" or the amino acids duplicated: {}'.format(change))
                del_aa = [follow[idx:idx+aa_nb_letter].capitalize() for idx in range(0, len(follow), aa_nb_letter)]
                for curr_del_aa in del_aa:
                    if curr_del_aa not in curr_aa_lexic:
                        raise Exception('Duplication must be ended by "dup" or the amino acids duplicated: {}'.format(change))
                follow = ""
        elif lc_follow == "ext" or lc_follow.startswith("ext-"):
            evt = "ext"
            follow = follow[3:]
        elif re.match(r"^(" + re_aa + ")ext", lc_follow):
            evt = "ext"
            new_aa = follow[:aa_nb_letter]
            follow = follow[aa_nb_letter + 3:]
        elif lc_follow == "fs":
            evt = "fs"
            follow = ""
        elif re.match(r"^(" + re_aa + ")fs", lc_follow):
            evt = "fs"
            new_aa = follow[:aa_nb_letter]
            follow = follow[aa_nb_letter + 2:]
        elif re.match(r"^fs(\*|ter)\d*$", lc_follow):
            evt = "fs"
            follow = follow[2:]
        return new_aa, evt, follow

    @staticmethod
    def _parsedNewElements(change, aa_nb_letter, evt, new_elements_str):
        """
        Return new elements list from part after evenment type in HGVS change (ex: "Glu" in Val600Glu or "*63" in Gln746_Lys747ins*63 or "*?" in Ile327Argfs*?).

        :param change: Change part of HGVS string (ex: "Val600Glu").
        :type change: str
        :param aa_nb_letter: Number of letters to represent an amino acid: 1 or 3.
        :type aa_nb_letter: int
        :param evt: Type of evenment in [None, "del", "dup", "ext", "fs", "ins", "insdel"].
        :type evt: str
        :param new_elements_str: Part after evenment type in HGVS change (ex: "Glu" in Val600Glu or "*63" in Gln746_Lys747ins*63 or "*?" in Ile327Argfs*?).
        :type new_elements_str: str
        :return: New elements list. Example ["Glu"] from "Glu" or ["*", "63"] from "*63" or ["*", "?"] from "*?".
        :rtype: list
        """
        new_elts = []
        remaining = deque(new_elements_str)
        curr_aa = ""
        while len(remaining):
            elt = remaining.popleft()
            if elt == "=":
                if len(remaining) > 0 and remaining[0] not in {"/", "^"}:
                    raise Exception('"=" it is not the last character in an alternative: {}'.format(change))
            elif elt == "*":
                elt = "Ter"
                if len(remaining) != 0:
                    if remaining[0] not in {"?", "^", "/"} and not remaining[0].isdigit():
                        raise Exception('"*" it is not the last element in change and the next element it is not in ["?", "/", "^", digit]: {}'.format(change))
            elif elt == "X":
                elt = "Xaa"
                if aa_nb_letter == 3:
                    if len(remaining) >= 2 and remaining[0] == "a" and remaining[1] == "a":
                        remaining.popleft()
                        remaining.popleft()
            elif elt in {"^", "/"}:
                if len(remaining) == 0:
                    raise Exception('"{}" cannot ends the change string: {}'.format(elt, change))
                ############################## TODO: parse evt in remaining: Val245=/del
            elif elt == "(":
                if evt != "ins":
                    raise Exception("Only insertions can contain brackets: {}.".format(change))
                while len(remaining) != 0 and remaining[0].isdigit():
                    elt += remaining.popleft()
                if elt == "(":
                    raise Exception("Only numbers are authorized in brackets: {}.".format(change))
                if len(remaining) == 0 or remaining[0] != ")":
                    raise Exception("Open bracket in insertion must be closed: {}.".format(change))
                elt += remaining.popleft()
            elif elt == "-":
                if evt != "ext":
                    raise Exception("Only extensions can contain minus: {}.".format(change))
                while len(remaining) != 0 and remaining[0].isdigit():
                    elt += remaining.popleft()
                if elt == "-":
                    raise Exception("Only numbers are authorized after minus: {}.".format(change))
            elif elt.isdigit():
                while len(remaining) != 0 and remaining[0].isdigit():
                    elt += remaining.popleft()
                if evt not in {"ins", "fs", "ext"}:
                    raise Exception('Digits are only authorized in change of type "ext", "fs" or "ins": {}.'.format(change))
                if evt != "ins":
                    if len(remaining) != 0 and remaining[0] != "/":
                        raise Exception('Digit are only authorized at the end of the change or followed by "/" for "fs" and "ext": {}'.format(change))
            elif elt == "?":
                if len(remaining) != 0:
                    raise Exception('"?" must be the last character in change: {}'.format(change))
            else:  # Is amino acid
                while len(elt) != aa_nb_letter:
                    elt += remaining.popleft()
                elt = elt.capitalize()
                if aa_nb_letter == 1:
                    elt = AA_THREE_BY_ONE[elt]
                else:
                    elt = elt.capitalize()
                    if elt not in THREE_LETTERS_AA_LEXIC:
                        raise Exception('The element "{}" in {} is not an amino acid.'.format(elt, change))
                if elt == "Ter":
                    if len(remaining) != 0:
                        if remaining[0] not in {"?", "^", "/"} and not remaining[0].isdigit():
                            raise Exception('"Ter" it is not the last alement in change and the next element it is not in ["?", "/", "^", digit]: {}'.format(change))
            new_elts.append(elt)
        if curr_aa != "":
            raise Exception("One of the amino acids is incomplete: {} in {}".format(curr_aa, change))
        return new_elts

    @staticmethod
    def isValid(change):
        """
        Return True if the string is a valid HGVS change (ex: "Val600Glu" or "(Val582_Asn583insXXXXX)").

        :param change: String to validate.
        :type change: str
        :return: True if the string is a valid HGVS change (ex: "Val600Glu" or "(Val582_Asn583insXXXXX)").
        :rtype: bool
        """
        is_hgvs = None
        try:
            HGVSProtChange.fromStr(change)
            is_hgvs = True
        except Exception:
            is_hgvs = False
        return is_hgvs

    @staticmethod
    def fromStr(change, aa_nb_letter=None):
        """
        Return HGVSProtChange from change part of HGVS string (ex: "Val600Glu").

        :param change: Change part of HGVS string (ex: "Val600Glu").
        :type change: str
        :param aa_nb_letter: Number of letters to represent an amino acid: 1 or 3. By default, the function tries to determine this number automatically.
        :type aa_nb_letter: int
        :return: HGVSProtChange from change part of HGVS string (ex: "Val600Glu").
        :rtype: anacore.hgvs.HGVSProtChange
        """
        predicted = None
        if change[0] == "(" and change[-1] == ")":
            change = change[1:-1]
            predicted = True
        else:
            predicted = False
        if change == "?":
            return HGVSProtChange(None, None, None, None, None, None, None, predicted)
        if change == "0":
            return HGVSProtChange(None, 0, None, None, None, None, None, predicted)
        match = re.search(r"^([^\d]+)(\d+)(.+)$", change)
        if not match:
            raise Exception('"{}" has not a valid syntax for HGVS protein change.'.format(change))
        else:
            start_aa, start_pos, follow = match.groups()
            start_pos = int(start_pos)
            start_aa.capitalize()
            if aa_nb_letter is None:
                aa_nb_letter = HGVSProtChange._getAANbLetter(start_aa, follow)
            end_aa = None
            end_pos = None
            if follow.startswith("_"):
                end_aa = follow[1:aa_nb_letter + 1]
                follow = follow[1 + aa_nb_letter:]
                end_pos = ""
                while follow[0].isdigit():
                    end_pos += follow[0]
                    follow = follow[1:]
                end_pos = int(end_pos)
            new_aa, evt, new_elts_str = HGVSProtChange._splittedOnEvt(change, aa_nb_letter, start_pos, end_pos, follow)
            new_elts = [elt for elt in new_elts_str]
            if len(new_elts_str) != 0:
                if not re.match(r"^\[\d+\]$", new_elts_str):  # The alteration is not a repeat
                    new_elts = HGVSProtChange._parsedNewElements(change, aa_nb_letter, evt, new_elts_str)
                else:  # The alteration is a repeat
                    nb_repeats = new_elts_str[1:-1]
                    if nb_repeats == "1":
                        new_elts = []
                        evt = "dup"
                    else:
                        new_elts = ["[", nb_repeats, "]"]
            return HGVSProtChange(start_aa, start_pos, end_aa, end_pos, new_aa, evt, new_elts, predicted)

    @staticmethod
    def insCouldBeIdentical(hgvs_ins, hgvs_repeat):
        """
        Return True if the insertion could be identical to the repeat (example: L16_I17insGTTL vs G13_L16dup). In this case, the difference bewtenn the two HGVS comes only from a difference of notation.

        :param hgvs_ins: HGVSp change for an insertion (example: L16_I17insGTTL).
        :type hgvs_ins: anacore.hgvs.HGVSProtChange
        :param hgvs_repeat: HGVSp change for a repeat (example: G13_L16dup or G13_L16[1]).
        :type hgvs_repeat: anacore.hgvs.HGVSProtChange
        :return: True if the insertion could be identical to the repeat.
        :rtype: bool
        """
        could_be_identical = False
        if hgvs_ins.isInFrameIns() and hgvs_repeat.isInFrameIns():
            repeat_end_pos = hgvs_repeat.end_pos if hgvs_repeat.end_pos is not None else hgvs_repeat.start_pos
            repeat_end_aa = hgvs_repeat.end_aa if hgvs_repeat.end_aa is not None else hgvs_repeat.start_aa
            len_dup = repeat_end_pos - hgvs_repeat.start_pos + 1
            nb_repeat = 1
            if hgvs_repeat.evt is None:  # Repeat
                nb_repeat = int(hgvs_repeat.new_elts[1])
            len_ins_in_repeat = len_dup * nb_repeat
            if len(hgvs_ins.new_elts) == len_ins_in_repeat:  # same length
                if repeat_end_aa + str(repeat_end_pos) == hgvs_ins.start_aa + str(hgvs_ins.start_pos):  # Same positions
                    if hgvs_repeat.start_aa == hgvs_ins.new_elts[0] and repeat_end_aa == hgvs_ins.new_elts[-1]:  # Same start and end amino acids
                        if nb_repeat == 1:
                            could_be_identical = True
                        else:
                            insertions = ["".join(hgvs_ins.new_elts[start:start + len_dup]) for start in range(0, len_ins_in_repeat, len_dup)]
                            if len(set(insertions)) == 1:  # Insertion is a valid repeat
                                could_be_identical = True
        return could_be_identical
