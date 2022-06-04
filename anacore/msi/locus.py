# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing MSI locus information (status, data, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from copy import deepcopy


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
        self.position = position
        self.name = name
        self.results = {} if results is None else results

    def delResult(self, result_id):
        """
        Delete result from locus.

        :param result_id: The ID used for index the results to delete.
        :type result_id: str
        """
        self.results.pop(result_id, None)

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
