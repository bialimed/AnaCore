# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing MSI locus information (status, data, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sequence import getShortestRepeatUnit
from copy import deepcopy
import re


def getRefSeqInfo(ref_fh, microsat, flank_size=5):
    """
    Return location, flanks bases, repeat sequences and number of occurences from reference sequence and target.

    :param ref_fh: File handle to reference sequences file.
    :type ref_fh: anacore.sequenceIO.IdxFastaIO
    :param microsat: Microsatellite region.
    :type microsat: anacore.region.Region
    :param flank_size: Size of flanking sequences.
    :type flank_size: int
    :return: Location, flanks bases, repeat sequences and number of occurences from reference sequence and target.
    :rtype: dict
    """
    left_flank = ref_fh.getSub(microsat.reference.name, microsat.start - flank_size, microsat.start - 1)
    right_flank = ref_fh.getSub(microsat.reference.name, microsat.end + 1, microsat.end + flank_size)
    target_seq = ref_fh.getSub(microsat.reference.name, microsat.start, microsat.end)
    repeat_unit = getShortestRepeatUnit(target_seq)
    if repeat_unit is None:
        raise Exception("The region {} doe not contain any repeat: {}.".format(microsat, target_seq))
    return {
        "chromosome": microsat.reference.name,
        "location": microsat.start - 1,
        "repeat_times": int(microsat.length() / len(repeat_unit)),
        "repeat_unit_bases": repeat_unit,
        "left_flank_bases": left_flank,
        "right_flank_bases": right_flank
    }


class Locus:
    """Manage a locus of a microsatellite (name, position and stability status)."""

    def __init__(self, position, name=None, results=None):
        """
        Build and return an instance of Locus.

        :param position: The position of the locus with format chr:start-end (the start is 0-based).
        :type position: str
        :param name: The name of the locus (for example: NR21).
        :type name: str
        :param results: The results of the locus indexed by the name of the method used to produce these results.
        :type results: dict
        :return: The new instance.
        :rtype: Locus
        """
        self.name = name
        if not re.search(r"^.+:\d+-\d+$", position):
            raise ValueError("{} is invalid for locus position. Locus position must be in format: chr:start-end.".format(position))
        self.position = position
        self.results = {} if results is None else results

    def delResult(self, result_id):
        """
        Delete result from locus.

        :param result_id: The ID used for index the results to delete.
        :type result_id: str
        """
        self.results.pop(result_id, None)

    @property
    def end(self):
        """
        Return end position on reference (1-based).

        :return: End position on reference (1-based).
        :rtype: int
        """
        return int(self.position.rsplit("-", 1)[1])

    @staticmethod
    def fromDict(hash):
        """
        Build and return an instance of Locus from a dict. This method is used for load instance from JSON.

        :param hash: The locus information.
        :type hash: dict
        :return: The new instance.
        :rtype: anacore.msi.base.Locus
        """
        cleaned_hash = deepcopy(hash)
        if "results" in hash:
            for method, result in hash["results"].items():
                cleaned_hash["results"][method] = LocusRes.fromDict(result)
        return Locus(**cleaned_hash)

    @property
    def length(self):
        """
        Return length of locus in reference sequence.

        :return: Length of locus in reference sequence.
        :rtype: int
        """
        return self.end - self.start + 1

    @property
    def start(self):
        """
        Return start position on reference (1-based).

        :return: Start position on reference (1-based).
        :rtype: int
        """
        return int(self.position.rsplit("-", 1)[0].rsplit(":", 1)[1]) + 1


class LocusDataDistrib:
    """Microsatellites lengths distribution."""

    def __init__(self, ct_by_len=None, mode="reads"):
        """
        Build and return an instance of LocusDataDistrib.

        :param ct_by_len: Depth by length.
        :type ct_by_len: dict
        :param mode: Depth calculation mode ("reads" or "fragments").
        :type mode: str
        :return: The new instance.
        :rtype: LocusDataDistrib
        """
        self.ct_by_len = {}
        if ct_by_len is not None:
            self.ct_by_len = {int(len): ct for len, ct in ct_by_len.items()}
        if mode not in {"reads", "fragments"}:
            raise Exception('Mode for distribution must be "reads" or "fragments" not: {}.'.format(mode))
        self.mode = mode

    @staticmethod
    def fromDense(dense, mode="reads", start=1):
        """
        Return LocusDataDistrib from dense list of counts (example: [0, 5, 4, 1] => {2: 5, 3: 4, 4: 1}).

        :param dense: List of counts. Each value correspond to a length and each length is represented between start and max length.
        :type dense: list
        :param mode: Depth calculation mode ("reads" or "fragments").
        :type mode: str
        :param start: Lis of counts start form this value. First length is more often 1 but can be 0.
        :type start: int
        :return: LocusDataDistrib from dense list of counts (example: [0, 5, 4, 1] => {2: 5, 3: 4, 4: 1}).
        :rtype: LocusDataDistrib
        """
        ct_by_len = {idx + start: ct for idx, ct in enumerate(dense) if ct != 0}
        return LocusDataDistrib(ct_by_len, mode)

    def getCount(self):
        """
        Return the number of elements in size distribution.

        :return: The number of elements in size distribution.
        :rtype: int
        """
        return sum(list(self.ct_by_len.values()))

    def getDenseCount(self, start=None, end=None):
        """
        Return the number of elements by locus length for absolutely all lengths betwen start and end. The length with 0 are also indicated.

        :param start: The first length of the returned distribution. [Default: the minimum length in the distribution]
        :type start: int
        :param end: The last length of the returned distribution. [Default: the maximum length in the distribution]
        :type end: int
        :return: The number of elements in each length.
        :rtype: list
        """
        if start is None:
            start = self.getMinLength()
        if end is None:
            end = self.getMaxLength()
        dense_count = list()
        for curr_length in range(start, end + 1):
            count = 0
            if curr_length in self.ct_by_len:
                count = self.ct_by_len[curr_length]
            dense_count.append(count)
        return dense_count

    def getDensePrct(self, start=None, end=None):
        """
        Return the percentage of elements by locus length for absolutely all lengths betwen start and end. The lengths with 0 are also indicated. The percentage is based on all the distribution not reduced at start and end.

        :param start: The first length of the returned distribution. [Default: the minimum length in the distribution]
        :type start: int
        :param end: The last length of the returned distribution. [Default: the maximum length in the distribution]
        :type end: int
        :return: The number of elements in each length.
        :rtype: list
        """
        nb_pairs = self.getCount()
        dense_prct = list()
        for curr_count in self.getDenseCount(start, end):
            prct = None
            if nb_pairs != 0:
                prct = (curr_count * 100) / nb_pairs
            dense_prct.append(prct)
        return dense_prct

    def getMaxLength(self):
        """
        Return the maximum length of the locus.

        :return: The maximum length of the locus.
        :rtype: int
        """
        return max([len for len, ct in self.ct_by_len.items() if ct != 0])

    def getMinLength(self):
        """
        Return the minimum length of the locus.

        :return: The minimum length of the locus.
        :rtype: int
        """
        return min([len for len, ct in self.ct_by_len.items() if ct != 0])

    def getMostRepresented(self):
        """
        Return length and count of the most represented length of the distribution. If two lengths have same count, the larger one is returned.

        :return: Length and count of the most represented length of the distribution. If two lengths have same count, the larger one is returned.
        :rtype: dict
        """
        highest_peak = {"length": None, "count": 0}
        for length, count in sorted(self.items()):
            if count >= highest_peak["count"]:
                highest_peak["length"] = length
                highest_peak["count"] = count
        return highest_peak

    def getNbPeaks(self, min_count=1):
        """
        Return number of peaks in the distribution.

        :param min_count: Length with lower count than this value are not taken into account.
        :type min_count: int
        :return: Number of peaks in the distribution.
        :rtype: int
        """
        nb_peak = 0
        for length, count in self.items():
            if count >= min_count:
                nb_peak += 1
        return nb_peak

    def items(self):
        """
        Return a generator on tuples (lengths, count). Lengths are sparse.

        :return: Generator on tuples (lengths, count). Lengths are sparse.
        :rtype: Generator
        """
        for len, ct in self.ct_by_len.items():
            yield (len, ct)


class LocusRes:
    """Manage the stability status for an anlysis of a locus."""

    def __init__(self, status, score=None, data=None):
        """
        Build and return an instance of LocusRes.

        :param status: The status of the locus.
        :type status: msi.Status
        :param score: The confidence score on the status prediction (0 to 1).
        :type score: float
        :param data: The data used to predict the status (example: the size distribution).
        :type results: dict
        :return: The new instance.
        :rtype: LocusRes
        """
        self.status = status
        self.score = score
        self.data = {} if data is None else data

    @staticmethod
    def fromDict(hash):
        """
        Build and return an instance of LocusRes from a dict. This method is used for load instance from JSON.

        :param hash: The locus result information.
        :type hash: dict
        :return: The new instance.
        :rtype: LocusRes
        """
        cleaned_hash = deepcopy(hash)
        if "data" in hash and "lengths" in hash["data"]:
            if not issubclass(hash["data"]["lengths"].__class__, LocusDataDistrib):
                cleaned_hash["data"]["lengths"] = LocusDataDistrib(**hash["data"]["lengths"])
        return LocusRes(**cleaned_hash)
