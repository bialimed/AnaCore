# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing Hubble's outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sv import HashedSVIO
from anacore.msi.base import Status
from anacore.msi.locus import Locus
from anacore.msi.sample import MSISample, MSISplRes
import os
import json


class HubbleDiff(HashedSVIO):
    """Manage differential analysis results file coming from Hubble."""

    def __init__(self, filepath, mode="r"):
        """
        Return the new instance of HubbleDiff.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: New instance of HubbleDiff.
        :rtype: HubbleDiff
        """
        super().__init__(filepath, mode, "\t", "#")
        if mode == "w":
            self.titles = ["Chromosome", "Start", "RepeatUnit", "Assessed", "Distance", "PValue"]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line.
        :rtype: anacore.msi.base.Locus
        """
        record = super()._parseLine()
        return Locus.fromDict({
            "position": "{}:{}".format(record["Chromosome"], record["Start"]),
            "results": {
                "Hubble": {
                    "status": Status.undetermined if record["Assessed"] == "False" else Status.none,
                    "score": None if record["Assessed"] == "False" else 1 - float(record["PValue"]),
                    "data": {
                        "nt": record["RepeatUnit"],
                        "distance": None if record["Assessed"] == "False" else float(record["Distance"]),
                        "p_value": None if record["Assessed"] == "False" else float(record["PValue"])
                    }
                }
            }
        })

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: Record containing differential analysis results coming from Hubble.
        :type record: anacore.msi.base.Locus
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        res = record.results["Hubble"]
        formatted_record = {
            "Chromosome": record.position.split(":")[0],
            "Start": record.position.split(":")[1],
            "RepeatUnit": res.data["nt"],
            "Assessed": "False" if res.data["distance"] is None else "True",
            "Distance": "" if res.data["distance"] is None else res.data["distance"],
            "PValue": "" if res.data["p_value"] is None else res.data["p_value"]
        }
        return super().recordToLine(formatted_record)


class HubbleDist(HashedSVIO):
    """Manage distributions file coming from hubble."""

    def __init__(self, filepath, mode="r", len_limit=100):
        """
        Return the new instance of HubbleDist.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param len_limit: Max length recorded in length_distribution.
        :type len_limit: int
        :return: New instance of HubbleDist.
        :rtype: HubbleDist
        """
        super().__init__(filepath, mode, "\t", "#")
        self.len_limit = len_limit
        if mode == "w":
            self.titles = ["chromosome", "location", "repeat_unit_bases", "covered", "length_distribution"]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line.
        :rtype: anacore.mis.base.Locus
        """
        record = super()._parseLine()
        # Parse distribution
        nb_by_length = {}
        for idx, count in enumerate(record["length_distribution"].split(",")):
            if count != "0":
                nb_by_length[idx + 1] = int(count)  # Lengths count start from 1
        # Locus
        return Locus.fromDict({
            "position": "{}:{}".format(record["chromosome"], record["location"]),
            "results": {
                "Hubble": {
                    "status": Status.none if record["covered"] == "True" else Status.undetermined,
                    "data": {
                        "nt": record["repeat_unit_bases"],
                        "lengths": {"ct_by_len": nb_by_length}
                    }
                }
            }
        })

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: Record containing distributions coming from Hubble.
        :type record: anacore.msi.base.Locus
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        res = record.results["Hubble"]
        formatted_record = {
            "chromosome": record.position.split(":")[0],
            "location": record.position.split(":")[1],
            "repeat_unit_bases": res.data["nt"],
            "covered": "False" if res.data["lengths"].getCount() == 0 else "True",
            "length_distribution": ",".join(
                map(str, res.data["lengths"].getDenseCount(1, self.len_limit))  # Lengths count start from 1
            )
        }
        return super().recordToLine(formatted_record)


def parseHubbleResults(summary_path, differential_path, distributions_path=None, unstable_threshold=0.2):
    """
    Return information on loci and analyses about microsatellites for sample.

    :param summary_path: Path to microsatellites analysis summary file (format: JSON).
    :type summary_path: str
    :param differential_path: Path to microsatellites diffrential analysis results file (format: JSON).
    :type differential_path: str
    :param distributions_path: Path to microsatellites sizes distributions file (format: JSON).
    :type distributions_path: str
    :param unstable_threshold: Over this rate of unstable loci the sample is tagged as unstable.
    :type unstable_threshold: float
    :return: Information on loci and analyses about microsatellites for sample.
    :rtype: anacore.msi.MSISpl
    """
    # Parse summary
    sum_data = None
    with open(summary_path) as reader:
        sum_data = json.load(reader)
    spl_name = os.path.basename(sum_data["Settings"]["DetectionOptions"]["QueryBam"])
    if spl_name.endswith(".bam"):
        spl_name = spl_name[:-4]
    msi_spl = MSISample(spl_name)
    params = sum_data["Settings"]["TumorOnlyComparisonOptions"]
    if params is None:
        params = sum_data["Settings"]["TumorNormalComparisonOptions"]
    status = Status.undetermined
    if sum_data["ResultIsValid"]:
        status = Status.none
        if unstable_threshold:
            if sum_data["PercentageUnstableSites"] > unstable_threshold * 100:
                status = Status.unstable
            else:
                status = Status.stable
    msi_spl.results["Hubble"] = MSISplRes(
        status,
        None,
        "Hubble",
        params
    )
    # Parse differential analysis results
    for curr_locus in HubbleDiff(differential_path):
        res_locus = curr_locus.results["Hubble"]
        if res_locus.status != Status.undetermined:  # Status can be assessed
            res_locus.status = Status.stable
            if res_locus.data["distance"] >= msi_spl.results["Hubble"].param["DistanceThreshold"]:
                if res_locus.data["p_value"] <= msi_spl.results["Hubble"].param["PValueThreshold"]:
                    res_locus.status = Status.unstable
        msi_spl.addLocus(curr_locus)
    # Check parsing
    if sum_data["TotalMicrosatelliteSitesAssessed"] != msi_spl.getNbDetermined("Hubble"):
        raise Exception(
            "Assessed microsat number is discordant betseen {} ({}) and {} ({}).".format(
                summary_path,
                sum_data["TotalMicrosatelliteSitesAssessed"],
                differential_path,
                msi_spl.getNbDetermined("Hubble")
            )
        )
    if sum_data["TotalMicrosatelliteSitesUnstable"] != msi_spl.getNbUnstable("Hubble"):
        raise Exception(
            "Unstable microsat number is discordant betseen {} ({}) and {} ({}).".format(
                summary_path,
                sum_data["TotalMicrosatelliteSitesUnstable"],
                differential_path,
                msi_spl.getNbUnstable("Hubble")
            )
        )
    # Parse distributions
    if distributions_path is not None:
        for curr_locus in HubbleDist(distributions_path):
            old_res = msi_spl.loci[curr_locus.position].results["Hubble"]
            curr_locus.results["Hubble"].status = old_res.status
            curr_locus.results["Hubble"].data["distance"] = old_res.data["distance"]
            curr_locus.results["Hubble"].data["p_value"] = old_res.data["p_value"]
            msi_spl.loci[curr_locus.position] = curr_locus
    return msi_spl
