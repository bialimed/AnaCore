# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing/processing a list of MSI samples."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.msi.base import getNbSupporting, Status, toDict
from anacore.msi.sample import MSISample
import json


class ReportIO:
    """Read/write the JSON file used to store a list of anacore.msi.sample.Sample."""

    @staticmethod
    def parse(in_path):
        """
        Return the list of MSI samples stored in a report file.

        :param in_path: Path to the file storing MSI samples (format: MSIReport).
        :type in_path: str
        :return: The list of MSI samples.
        :rtype: list
        """
        spl_data = []
        with open(in_path) as FH_in:
            spl_data = json.load(FH_in)
        return [MSISample.fromDict(curr_spl) for curr_spl in spl_data]

    @staticmethod
    def write(msi_samples, out_path):
        """
        Write the list of MSI samples in a MSI report file.

        :param msi_samples: The list of MSI samples.
        :type msi_samples: list
        :param out_path: Path to the output file storing MSI samples (format: MSIReport).
        :type in_path: str
        """
        with open(out_path, "w") as FH_out:
            json.dump(msi_samples, FH_out, sort_keys=True, default=toDict)


def getIncompleteModels(in_report, min_support_model=20):
    """
    Return the list of locus class model with an unsufficient number of supporting samples.

    :param in_report: Path to a MSI report file.
    :type in_report: str
    :param min_support_model: The minimum number of sample in locus class to validate the model.
    :type min_support_model: int
    :return: The list of locus class model with with an unsufficient number of supporting samples. Each item is dict with following keys: locus_name, locus_id, status and support.
    :rtype: list
    """
    incomplete_models = []
    report = ReportIO.parse(in_report)
    for curr_model in getNbSupporting(report, method="model"):
        if curr_model["status"] in [Status.stable, Status.unstable] and curr_model["support"] < min_support_model:
            incomplete_models.append(curr_model)
    return incomplete_models
