# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing MSI information (status, data, locus, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.7.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import math
import json
import numpy as np
from copy import deepcopy


def toDict(msi_object):
    """
    Retuns a dictionary representing the object. This method is used in json.dump in argument named "default" for recursively convert an object to json.

    :param msi_object: The object to convert.
    :type msi_object: a class of anacore.msi library
    :return: The dictionary representing the object.
    :rtype: dict
    """
    return msi_object.__dict__


class Status:
    """Status for samples (MSISplRes) and loci (LocusRes) in anacore.msi library."""

    stable = "MSS"
    unstable = "MSI"
    undetermined = "Undetermined"  # Cannot be determined
    none = None  # Not evaluated

    @staticmethod
    def authorizedValues():
        """Return the values authorized for stability status."""
        return [attr_value for attr_name, attr_value in Status.__dict__.items() if attr_name != "authorizedValues" and not attr_name.startswith("__")]


class MSILocus:
    """Manage a locus of a microsatellite (name, position and stability status)."""

    def __init__(self, position, name=None, results=None):
        """
        Build and return an instance of MSILocus.

        :param position: The position of the locus with format chr:start-end (the start is 0-based).
        :type position: str
        :param name: The name of the locus (for example: NR21).
        :type name: str
        :param results: The results of the locus indexed by the name of the method used to produce these results.
        :type results: dict
        :return: The new instance.
        :rtype: MSILocus
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
    def fromDict(data):
        """
        Build and return an instance of MSILocus from a dict. This method is used for load instance from JSON.

        :param data: The locus information.
        :type data: dict
        :return: The new instance.
        :rtype: MSILocus
        """
        cleaned_data = deepcopy(data)
        if "results" in data:
            for method, result in data["results"].items():
                if "_class" in result and result["_class"] == "LocusResPairsCombi":
                    cleaned_data["results"][method] = LocusResPairsCombi.fromDict(result)
                elif "_class" in result and result["_class"] == "LocusResDistrib":
                    cleaned_data["results"][method] = LocusResDistrib.fromDict(result)
                else:
                    cleaned_data["results"][method] = LocusRes.fromDict(result)
        return MSILocus(**cleaned_data)


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
        self._class = "LocusRes"
        self.status = status
        self.score = score
        self.data = {} if data is None else data

    @staticmethod
    def fromDict(data):
        """
        Build and return an instance of LocusRes from a dict. This method is used for load instance from JSON.

        :param data: The locus result information.
        :type data: dict
        :return: The new instance.
        :rtype: LocusRes
        """
        cleaned_data = deepcopy(data)
        if "_class" in cleaned_data:
            cleaned_data.pop("_class", None)
        return LocusRes(**cleaned_data)


class LocusResDistrib(LocusRes):
    """Manage the stability status for an anlysis of a locus containing the count by length."""

    def __init__(self, status, score=None, data=None):
        """
        Build and return an instance of LocusResDistrib.

        :param status: The status of the locus.
        :type status: msi.Status
        :param score: The confidence score on the status prediction (0 to 1).
        :type score: float
        :param data: The data used to predict the status (example: the size distribution).
        :type results: dict
        :return: The new instance.
        :rtype: LocusResDistrib
        """
        super().__init__(status, score, data)
        self._class = "LocusResDistrib"

    def getCount(self):
        """
        Return the number of elements in size distribution.

        :return: The number of elements in size distribution.
        :rtype: int
        """
        return sum(list(self.data["nb_by_length"].values()))

    def getMinLength(self):
        """
        Return the minimum length of the locus.

        :return: The minimum length of the locus.
        :rtype: int
        """
        return min([int(elt) for elt in self.data["nb_by_length"]])

    def getMaxLength(self):
        """
        Return the maximum length of the locus.

        :return: The maximum length of the locus.
        :rtype: int
        """
        return max([int(elt) for elt in self.data["nb_by_length"]])

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
            if str(curr_length) in self.data["nb_by_length"]:
                count = self.data["nb_by_length"][str(curr_length)]
            dense_count.append(count)
        return dense_count

    @staticmethod
    def fromDict(data):
        """
        Build and return an instance of LocusResDistrib from a dict. This method is used for load instance from JSON.

        :param data: The locus result information.
        :type data: dict
        :return: The new instance.
        :rtype: LocusResDistrib
        """
        cleaned_data = deepcopy(data)
        if "_class" in cleaned_data:
            cleaned_data.pop("_class", None)
        return LocusResDistrib(**cleaned_data)


class LocusResPairsCombi(LocusResDistrib):
    """Manage the stability status for an anlysis after reads pair combination of a locus."""

    def __init__(self, status, score=None, data=None):
        """
        Build and return an instance of LocusResPairsCombi.

        :param status: The status of the locus.
        :type status: msi.Status
        :param score: The confidence score on the status prediction (0 to 1).
        :type score: float
        :param data: The data used to predict the status (example: the size distribution).
        :type results: dict
        :return: The new instance.
        :rtype: LocusResPairsCombi
        """
        super().__init__(status, score, data)
        self._class = "LocusResPairsCombi"

    def getNbFrag(self):
        """
        Return the number of fragments in size distribution.

        :return: The number of fragments in size distribution.
        :rtype: int
        """
        return self.getCount()

    @staticmethod
    def fromDict(data):
        """
        Build and return an instance of LocusResPairsCombi from a dict. This method is used for load instance from JSON.

        :param data: The locus result information.
        :type data: dict
        :return: The new instance.
        :rtype: LocusResPairsCombi
        """
        cleaned_data = deepcopy(data)
        if "_class" in cleaned_data:
            cleaned_data.pop("_class", None)
        return LocusResPairsCombi(**cleaned_data)


class MSISplRes:
    """Manage the stability status for an anlysis of a sample."""

    def __init__(self, status, score=None, method=None, param=None, version=None):
        """
        Build and return an instance of MSISplRes.

        :param status: The status of the sample.
        :type status: msi.Status
        :param score: The confidence score on the status prediction (0 to 1).
        :type score: float
        :param method: The name of the method used to predict status.
        :type method: str
        :param param: The parameters of the status prediction.
        :type param: dict
        :param version: The version of the method used to predict status.
        :type version: dict
        :return: The new instance.
        :rtype: MSISplRes
        """
        self.status = status
        self.method = method
        self.score = score
        self.param = param
        self.version = version

    @staticmethod
    def fromDict(data):
        """
        Build and return an instance of MSISplRes from a dict. This method is used for load instance from JSON.

        :param data: The sample result information.
        :type data: dict
        :return: The new instance.
        :rtype: MSISplRes
        """
        cleaned_data = deepcopy(data)
        return MSISplRes(**cleaned_data)


class MSISample:
    """Manage a sample in context of microsatellite instability analysis. This object contains information on loci and on analyses (at sample and loci levels)."""

    def __init__(self, name, loci=None, results=None):
        """
        Build and return an instance of MSISample.

        :param name: The sample name.
        :type name: str
        :param loci: By ID (format: chr:start-end with start 0-based) the loci targeted in analysis.
        :type loci: dict
        :param results: By method name the analyses of sample status.
        :type results: dict
        :return: The new instance.
        :rtype: MSISample
        """
        self.name = name
        self.loci = {} if loci is None else loci
        self.results = {} if results is None else results

    def getLociMethods(self):
        """
        Return the different methods used in loci status analysis.

        :return: The names of the methods used in analyses of loci.
        :rtype: set
        """
        methods = set()
        for curr_locus in self.loci:
            for curr_method in self.results.keys():
                methods.add(curr_method)
        return methods

    def addLocus(self, locus):
        """
        Add one locus to the sample.

        :param locus: The locus to add.
        :type locus: MSILocus
        """
        if locus.__class__.__name__ != "MSILocus":
            raise Exception('The class "{}" cannot be used as locus for MSISample.'.format(locus.__class__.__name__))
        self.loci[locus.position] = locus

    def delLoci(self, locus_ids):
        """
        Delete loci from the sample.

        :param locus_ids: The IDs of the loci to delete.
        :type locus_ids: list
        """
        for locus_id in locus_ids:
            self.delLocus(locus_id)

    def delLocus(self, locus_id):
        """
        Delete one locus from the sample.

        :param locus_id: The ID of the locus to delete.
        :type locus_id: str
        """
        self.loci.pop(locus_id, None)

    def _getStatusByMethod(self, method):
        """
        Return the list of loci status predicted by the selected method.

        :param method: The selected method.
        :type method: str
        :return: The list of status.
        :rtype: list
        """
        status = list()
        for locus_id, locus in self.loci.items():
            if method in locus.results:
                status.append(locus.results[method].status)
        return status

    def getNbUnstable(self, method):
        """
        Return the number of loci predicted as unstable with the selected method.

        :param method: The selected method.
        :type method: str
        :return: The number of unstable loci.
        :rtype: int
        """
        nb_unstable = 0
        status = self._getStatusByMethod(method)
        for curr_status in status:
            if curr_status is not None and curr_status == Status.unstable:
                nb_unstable += 1
        return nb_unstable

    def getNbStable(self, method):
        """
        Return the number of loci predicted as stable with the selected method.

        :param method: The selected method.
        :type method: str
        :return: The number of stable loci.
        :rtype: int
        """
        nb_stable = 0
        status = self._getStatusByMethod(method)
        for curr_status in status:
            if curr_status is not None and curr_status == Status.stable:
                nb_stable += 1
        return nb_stable

    def getNbUndetermined(self, method):
        """
        Return the number of loci where the status is undetermined after prediction with the selected method.

        :param method: The selected method.
        :type method: str
        :return: The number of undetermined loci.
        :rtype: int
        """
        nb_undetermined = 0
        status = self._getStatusByMethod(method)
        for curr_status in status:
            if curr_status is not None and curr_status == Status.undetermined:
                nb_undetermined += 1
        return nb_undetermined

    def getNbDetermined(self, method):
        """
        Return the number of loci with a prediction of status with the selected method.

        :param method: The selected method.
        :type method: str
        :return: The number of determined loci.
        :rtype: int
        """
        nb_determined = 0
        status = self._getStatusByMethod(method)
        for curr_status in status:
            if curr_status is not None and curr_status != Status.undetermined:
                nb_determined += 1
        return nb_determined

    def getNbProcessed(self, method):
        """
        Return the number of loci where the status has been evaluated for the sample with the selected method.

        :param method: The selected method.
        :type method: str
        :return: The number of processed loci.
        :rtype: int
        """
        nb_processed = 0
        status = self._getStatusByMethod(method)
        for curr_status in status:
            if curr_status is not None:
                nb_processed += 1
        return nb_processed

    def getNbLoci(self):
        """
        Return the number of loci in the sample.

        :param method: The selected method.
        :type method: str
        :return: The number of loci.
        :rtype: int
        """
        return len(self.loci)

    @staticmethod
    def fromDict(data):
        """
        Build and return an instance of MSISample from a dict. This method is used for load instance from JSON.

        :param data: The sample information.
        :type data: dict
        :return: The new instance.
        :rtype: MSISample
        """
        cleaned_data = deepcopy(data)
        # Loci
        if "loci" in data:
            cleaned_data["loci"] = dict()
            for locus_id, locus in data["loci"].items():
                cleaned_data["loci"][locus_id] = MSILocus.fromDict(locus)
        # Results
        if "results" in data:
            cleaned_data["results"] = dict()
            for method, result in data["results"].items():
                cleaned_data["results"][method] = MSISplRes.fromDict(result)
        # Name
        return MSISample(**cleaned_data)

    def _getScoreCalculation(self, eval_status, method, undetermined_weight=0.5, locus_weight_is_score=True):
        """
        Calculate and return a confidence score for the sample status prediction. This score is calculation take into account each locus and his score with the following formula: sum(scores) / (len(scores) + nb_loci_undetermined * undetermined_weight).

        :param eval_status: The score is calculated for this status (msi.Status.stable or msi.Status.unstable). Except for undetermined_weight equals to 0 the score of the complementary status CANNOT be calculated by 1 - complementary score.
        :type eval_status: msi.Status
        :param method: The status of the loci are extracted from the results of this method.
        :type method: str
        :param undetermined_weight: The weight of the undetermined loci in score calculation.
        :type undetermined_weight: float
        :param locus_weight_is_score: Use the prediction score of each locus as wheight of this locus in score calculation.
        :type locus_weight_is_score: bool
        :return: The prediction score. This score has a value between 0 and 1.
        :rtype: float
        """
        scores = list()
        nb_loci_undetermined = 0
        for locus_id, locus in self.loci.items():
            if locus.results[method].status is not None:
                if locus.results[method].status == Status.undetermined:
                    nb_loci_undetermined += 1
                else:
                    if locus.results[method].status == eval_status:
                        if locus_weight_is_score and locus.results[method].score is not None:
                            scores.append(locus.results[method].score)
                        else:
                            scores.append(1)
                    elif locus.results[method].status != Status.undetermined:
                        if locus_weight_is_score and locus.results[method].score is not None:
                            scores.append(1 - locus.results[method].score)
                        else:
                            scores.append(0)
        score = None
        if len(scores) != 0:
            score = sum(scores) / (len(scores) + nb_loci_undetermined * undetermined_weight)
        return round(score, 5)

    def setScore(self, method, undetermined_weight=0.5, locus_weight_is_score=True):
        """
        Calculate and set a confidence score for the sample status prediction. This score is calculation take into account each locus and his score with the following formula: sum(scores) / (len(scores) + nb_loci_undetermined * undetermined_weight).

        :param method: Name of the selected method.
        :type method: str
        :param undetermined_weight: The weight of the undetermined loci in score calculation.
        :type undetermined_weight: float
        :param locus_weight_is_score: Use the prediction score of each locus as wheight of this locus in score calculation.
        :type locus_weight_is_score: bool
        """
        spl_res = self.results[method]
        if spl_res.status == Status.undetermined or spl_res.status == Status.none:
            spl_res.score = None
        else:
            spl_res.score = self._getScoreCalculation(spl_res.status, method, undetermined_weight)

    def setStatusByInstabilityCount(self, method, min_voting_loci=3, instability_threshold=3):
        """
        Add status and score for the sample by number of unstable loci for the selected method.

        :param method: Name of the selected method.
        :type method: str
        :param min_voting_loci: Minimum number of voting loci (stable + unstable) to determine the sample status. If the number of voting loci is lower than this value the status for the sample will be undetermined.
        :type min_voting_loci: int
        :param instability_threshold: If the number of unstable is superior or equal than this value the status of the sample will be unstable.
        :type instability_threshold: int
        """
        result = MSISplRes.fromDict({
            'status': Status.undetermined,
            'method': method,
            'score': None,
            'param': {"aggregation_method": "instability count", "min_voting_loci": min_voting_loci, "instability_threshold": instability_threshold},
            'version': "1.0.0"
        })
        nb_stable = self.getNbStable(method)
        nb_unstable = self.getNbUnstable(method)
        result.status = Status.undetermined
        if nb_stable + nb_unstable >= min_voting_loci:
            if nb_unstable >= instability_threshold:
                result.status = Status.unstable
            else:
                result.status = Status.stable
        self.results[method] = result
        self.setScore(method)

    def setStatusByInstabilityRatio(self, method, min_voting_loci=1, instability_threshold=0.4):
        """
        Add status and score for the sample by ratio of unstable loci for the selected method.

        :param method: Name of the selected method.
        :type method: str
        :param min_voting_loci: Minimum number of voting loci (stable + unstable) to determine the sample status. If the number of voting loci is lower than this value the status for the sample will be undetermined.
        :type min_voting_loci: int
        :param instability_threshold: If the ratio unstable/(stable + unstable) is superior or equal than this value the status of the sample will be unstable.
        :type instability_threshold: float
        """
        result = MSISplRes.fromDict({
            'status': Status.undetermined,
            'method': method,
            'score': None,
            'param': {"aggregation_method": "instability ratio", "min_voting_loci": min_voting_loci, "instability_threshold": instability_threshold},
            'version': "1.0.0"
        })
        nb_stable = self.getNbStable(method)
        nb_unstable = self.getNbUnstable(method)
        result.status = Status.undetermined
        if nb_stable + nb_unstable >= min_voting_loci:
            if nb_unstable / (nb_stable + nb_unstable) >= instability_threshold:
                result.status = Status.unstable
            else:
                result.status = Status.stable
        self.results[method] = result
        self.setScore(method)

    def setStatusByMajority(self, method, min_voting_loci=1):
        """
        Add status and score for the sample by consensus on loci results for the selected method.

        :param method: Name of the selected method.
        :type method: str
        :param min_voting_loci: Minimum number of voting loci (stable + unstable) to determine the sample status. If the number of voting loci is lower than this value the status for the sample will be undetermined.
        :type min_voting_loci: int
        """
        result = MSISplRes.fromDict({
            'status': Status.undetermined,
            'method': method,
            'score': None,
            'param': {"aggregation_method": "majority", "min_voting_loci": min_voting_loci},
            'version': "1.0.0"
        })
        nb_stable = self.getNbStable(method)
        nb_unstable = self.getNbUnstable(method)
        if nb_stable + nb_unstable >= min_voting_loci:
            result.status = Status.undetermined
            if nb_stable > nb_unstable:
                result.status = Status.stable
            elif nb_stable < nb_unstable:
                result.status = Status.unstable
        self.results[method] = result
        self.setScore(method)


class LocusClassifier:
    """
    Classifier for locus using MSISample objects.

    Synopsis:
        clf = LocusClassifier(locus_id, method_name, classifier)
        clf.fit(train_dataset)
        clf.predict(test_dataset)
        clf.predict_proba(test_dataset)


        clf = LocusClassifier(locus_id, method_name, classifier)
        clf.fit(train_dataset)
        clf.set_status(test_dataset)
    """

    def __init__(self, locus_id, method_name, classifier, model_method_name="model", data_method_name=None):
        """
        Build and return an instance of LocusClassifier.

        :param locus_id: The ID (format: chr:start-end with start 0-based) of the loci targeted in analysis.
        :type locus_id: str
        :param method_name: By method name the analyses of sample status.
        :type method_name: dict
        :param classifier: The classifier used to predict status.
        :type classifier: A sklearn classifier object
        :param model_method_name: The name of the method used in models loci to store expected status and lengths distributions.
        :type model_method_name: str
        :param data_method_name: The data used for the prediction are extracted from the results of this method. By default the selected method is the same of the classifier method_name.
        :type data_method_name: str
        :return: The new instance.
        :rtype: LocusClassifier
        """
        self.classifier = classifier
        self.data_method_name = method_name if data_method_name is None else data_method_name
        self.locus_id = locus_id
        self.method_name = method_name
        self.model_method_name = model_method_name
        self._max_len = None  # This value is used to uniformise lengths distributions in the comparison
        self._min_len = None  # This value is used to uniformise lengths distributions in the comparison
        self._test_dataset = []
        self._train_dataset = []  # The MSISample provided for fitting. These samples are filtered before fit (see _usable_train_dataset)
        self._usable_train_dataset = []  # The MSISample containing status information for the locus

    def _get_min_max_len(self, dataset, method):
        """
        Return the minimum and maximum length for the locus in the dataset.

        :param dataset: The list of MSISample.
        :type dataset: list
        :param method: The lengths distribution of the loci are extracted from the results of this method.
        :type method: str
        :return: The list of MSISample.
        :rtype: list
        """
        min_len = math.inf
        max_len = -1
        for curr_spl in dataset:
            locus_res = curr_spl.loci[self.locus_id].results[method]
            min_len = min(min_len, locus_res.getMinLength())
            max_len = max(max_len, locus_res.getMaxLength())
        return min_len, max_len

    def _set_min_max_len(self):
        """Set the minimum and maximum length for the locus in the dataset (train + test)."""
        train_min, train_max = self._get_min_max_len(self._usable_train_dataset, self.model_method_name)
        test_min, test_max = self._get_min_max_len(self._test_dataset, self.data_method_name)
        self._min_len = min(train_min, test_min)
        self._max_len = max(train_max, test_max)

    def _get_data(self, dataset, method):
        """
        Return uniformised lengths distribution from the dataset.

        :param dataset: The list of MSISample.
        :type dataset: list
        :param method: The lengths distribution of the loci are extracted from the results of this method.
        :type method: str
        :return: The uniformised lengths distribution. Rows are samples, columns are lengths and values are percentages of counts in the sample.
        :rtype: np.matrix
        """
        prct_matrix = []  # rows are samples, columns are lengths and values are percentages of counts in the sample.
        if self._min_len is None and self._max_len is None:
            self._set_min_max_len()
        for curr_spl in dataset:
            locus_res = curr_spl.loci[self.locus_id].results[method]
            prct_matrix.append(
                locus_res.getDensePrct(self._min_len, self._max_len)
            )
        return np.matrix(prct_matrix)

    def _get_test_data(self):
        """
        Return uniformised lengths distribution for samples in test dataset.

        :return: The uniformised lengths distribution. Rows are samples, columns are lengths and values are percentages of counts in the sample.
        :rtype: np.matrix
        """
        return self._get_data(self._test_dataset, self.data_method_name)

    def _get_train_data(self):
        """
        Return uniformised lengths distribution for samples in usable train dataset.

        :return: The uniformised lengths distribution. Rows are samples, columns are lengths and values are percentages of counts in the sample.
        :rtype: np.matrix
        """
        return self._get_data(self._usable_train_dataset, self.model_method_name)

    def _get_train_labels(self):
        """
        Return the list of labels for samples in usable train dataset.

        :return: The list of labels for samples in usable train dataset.
        :rtype: np.array
        """
        labels = []
        for curr_spl in self._usable_train_dataset:
            locus_res = curr_spl.loci[self.locus_id].results[self.model_method_name]
            labels.append(locus_res.status)
        return np.array(labels)

    def fit(self, train_dataset):
        """
        Fit the model using train_dataset as training data and their status as target values.

        :param train_dataset: The list of MSISample containing the locus to classify and in LocusRes the data of the selected method and the status ecpected.
        :type test_dataset: list
        """
        self._train_dataset = train_dataset
        self._usable_train_dataset = [spl for spl in train_dataset if self.model_method_name in spl.loci[self.locus_id].results]
        self.classifier.fit(self._get_train_data(), self._get_train_labels())

    def predict(self, test_dataset):
        """
        Predict the status of the locus in the list of samples.

        :param test_dataset: The list of MSISample containing the locus to classify and in LocusRes the data of the selected method. These samples must already be filtered, for eaxmple on count of elements in length distribution.
        :type test_dataset: list
        :return: The predicted class for the selected locus in each sample.
        :rtype: np.array
        """
        self._test_dataset = test_dataset
        test_min_len, test_max_len = self._get_min_max_len(self._test_dataset, self.data_method_name)
        if test_min_len < self._min_len or test_max_len > self._max_len:
            self.fit(self._train_dataset)
        return self.classifier.predict(self._get_test_data())

    def predict_proba(self, test_dataset):
        """
        Return estimated probabilities for the test_dataset.

        :param test_dataset: The list of MSISample containing the locus to classify and in LocusRes the data of the selected method. These samples must already be filtered, for eaxmple on count of elements in length distribution.
        :type test_dataset: list
        :return: The estimated probabilities for the test_dataset.
        :rtype: np.array
        """
        self._test_dataset = test_dataset
        return self.classifier.predict_proba(self._get_test_data())

    def _get_scores(self, pred_labels):
        """
        Return confidence scores for the predictions of the test_dataset.

        :param pred_labels: The list of predicted status for the test_dataset.
        :type pred_labels: list
        :return: The scores for the predictions. The values are between 0 and 1.
        :rtype: list
        """
        scores = None
        proba_idx_by_label = {label: idx for idx, label in enumerate(self.classifier.classes_)}
        try:
            proba = self.predict_proba(self._test_dataset)
            scores = [round(spl_proba[proba_idx_by_label[spl_label]], 6) for spl_proba, spl_label in zip(proba, pred_labels)]
        except Exception:
            scores = [None for spl in self._test_dataset]
        return scores

    def set_status(self, test_dataset):
        """
        Set status and score for the selected locus in the list of samples.

        :param test_dataset: The list of MSISample containing the locus to classify and in LocusRes the data of the selected method. These samples must already be filtered, for eaxmple on count of elements in length distribution.
        :type test_dataset: list
        """
        self._test_dataset = test_dataset
        pred_labels = self.predict(test_dataset)
        pred_scores = self._get_scores(pred_labels)
        for label, score, sample in zip(pred_labels, pred_scores, self._test_dataset):
            if self.method_name not in sample.loci[self.locus_id].results:  # If the method does not exist in locus results
                sample.loci[self.locus_id].results[self.method_name] = LocusRes(Status.none)
            locus_res = sample.loci[self.locus_id].results[self.method_name]
            if self.method_name != self.data_method_name:  # If data used in prediction come from an other method
                locus_res.data["data_source"] = self.data_method_name
            locus_res.status = label
            locus_res.score = score


class MSIReport:
    """Read/write the JSON file used to store a list of MSISample."""

    @staticmethod
    def parse(in_path):
        """
        Return the list of MSISample stored in a MSIReport file.

        :param in_path: Path to the file storing MSISamples (format: MSIReport).
        :type in_path: str
        :return: The list of MSISample.
        :rtype: list
        """
        spl_data = []
        with open(in_path) as FH_in:
            spl_data = json.load(FH_in)
        return [MSISample.fromDict(curr_spl) for curr_spl in spl_data]

    @staticmethod
    def write(msi_samples, out_path):
        """
        Write the list of MSISample in a MSIReport file.

        :param msi_samples: The list of MSISample.
        :type msi_samples: list
        :param out_path: Path to the output file storing MSISamples (format: MSIReport).
        :type in_path: str
        """
        with open(out_path, "w") as FH_out:
            json.dump(msi_samples, FH_out, sort_keys=True, default=toDict)


def getNbSupporting(report, method="model", loci=None):
    """
    Return the number of samples by locus class.

    :param in_report: List of MSISample.
    :type in_report: list()
    :param method: Evaluated method.
    :type method: str
    :param loci: List of evaluated loci.
    :type loci: list
    :return: The number of samples by locus class. Each item is dict with following keys: locus_name, locus_id, status and support.
    :rtype: list
    """
    # Init loci from report if the argument is None
    if loci is None:
        loci = set()
        for spl in report:
            for locus_id in spl.loci:
                loci.add(locus_id)
    # Get count by status by locus
    nb_by_locus = {}
    for spl in report:
        for locus_id in loci:
            locus = spl.loci[locus_id]
            if locus_id not in nb_by_locus:
                nb_by_locus[locus_id] = {
                    "locus_name": locus.name,
                    "supp_by_status": {elt: 0 for elt in Status.authorizedValues()}
                }
            if method in locus.results:
                locus_status = locus.results[method].status
                nb_by_locus[locus_id]["supp_by_status"][locus_status] += 1
    # To list
    counts = []
    for locus_id, locus_info in nb_by_locus.items():
        for status, support in locus_info["supp_by_status"].items():
            counts.append({
                "locus_id": locus_id,
                "locus_name": locus_info["locus_name"],
                "status": status,
                "support": support
            })
    return counts


def getIncompleteModels(in_report, min_support_model=20):
    """
    Return the list of locus class model with an unsufficient number of supporting samples.

    :param in_report: Path to a MSIReport file.
    :type in_report: str
    :param min_support_model: The minimum number of sample in locus class to validate the model.
    :type min_support_model: int
    :return: The list of locus class model with with an unsufficient number of supporting samples. Each item is dict with following keys: locus_name, locus_id, status and support.
    :rtype: list
    """
    incomplete_models = []
    report = MSIReport.parse(in_report)
    for curr_model in getNbSupporting(report, method="model"):
        if curr_model["status"] in [Status.stable, Status.unstable] and curr_model["support"] < min_support_model:
            incomplete_models.append(curr_model)
    return incomplete_models
