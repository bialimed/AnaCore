# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing MSI sample information (status, data, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.1.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from copy import deepcopy
from anacore.msi.base import Status
from anacore.msi.locus import Locus


class MSISplRes:
    """Manage the stability status for an anlysis of a sample."""

    def __init__(self, status, score=None, method=None, param=None, version=None):
        """
        Build and return an instance of MSISplRes.

        :param status: The status of the sample.
        :type status: anacore.msi.base.Status
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
        :type locus: anacore.msi.base.Locus
        """
        if locus.__class__.__name__ != "Locus":
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
                cleaned_data["loci"][locus_id] = Locus.fromDict(locus)
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
            spl_res.score = self._getScoreCalculation(spl_res.status, method, undetermined_weight, locus_weight_is_score)

    def setStatusByInstabilityRatio(self, method, min_voting_loci=0.5, instability_threshold=0.4):
        """
        Add status and score for the sample by ratio of unstable loci for the selected method.

        :param method: Name of the selected method.
        :type method: str
        :param min_voting_loci: Minimum ratio (stable + unstable)/all to determine the sample status. If the number of voting loci is lower than this value the status for the sample will be undetermined.
        :type min_voting_loci: float
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
        if (nb_stable + nb_unstable) / self.getNbLoci() >= min_voting_loci:
            if nb_unstable / (nb_stable + nb_unstable) >= instability_threshold:
                result.status = Status.unstable
            else:
                result.status = Status.stable
        self.results[method] = result
        self.setScore(method)
