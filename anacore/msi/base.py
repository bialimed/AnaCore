# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing MSI information (status, data, ...)."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.msi.locus import LocusRes
import math
import numpy as np


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
        """
        Return authorized values for stability status.

        :return: Authorized for stability status.
        :rtype: set
        """
        return {attr_value for attr_name, attr_value in Status.__dict__.items() if attr_name != "authorizedValues" and not attr_name.startswith("__")}


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
            min_len = min(min_len, locus_res.data["lengths"].getMinLength())
            max_len = max(max_len, locus_res.data["lengths"].getMaxLength())
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
                locus_res.data["lengths"].getDensePrct(self._min_len, self._max_len)
            )
        return np.asarray(prct_matrix)

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


def getNbSupporting(report, method="model", loci=None):
    """
    Return the number of samples by locus class.

    :param in_report: List of MSI samples.
    :type in_report: list
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
