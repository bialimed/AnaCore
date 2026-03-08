# -*- coding: utf-8 -*-
"""
Classes and functions for reading Sequel's files on run.

:Code example:

    Check if run is ended

    .. highlight:: python
    .. code-block:: python

        from anacore.instrument.sequel.run import Run

        run = Run("my_run/"):
        print("Run steps:")
        print("", "1- sequencing ended:", run.isSequenced())
        print("", "2- copy ended:", run.isCopied())

        # Result>
        # Run steps:
        #   1- sequencing ended: True
        #   2- copy ended: False
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import datetime
import glob
import os
import re
import xml.etree.ElementTree as ET
# https://pacbiofileformats.readthedocs.io/en/11.0/RunDesignCsv.html
# https://www.umassmed.edu/globalassets/deep-sequencing-core/pbce/sequel_ii_and_iie_data_files.pdf


def _consReadsetDateToDatetime(date_str):
    """
    Return datetime object from date string coming from consensus readset file.

    :param date_str: Date string coming from consensus readset file (example: 2026-01-01T10:15:30.3Z).
    :type date_str: str
    :return: Datetime object from date string coming from consensus readset file.
    :rtypr: datetime.datetime
    """
    iso_date_str = date_str.replace("Z", "+00:00")
    iso_date_str = re.sub(r"\.(\d*)\+", r"+", iso_date_str)
    return datetime.datetime.fromisoformat(iso_date_str)


def etreeToDict(node):
    """
    Return dict from XML node.

    :param node: Node to convert.
    :type node: xml.etree.ElementTree.Element
    :return: Dict corresponding to XML node.
    :rtype: dict
    """
    data = {"key": node.tag.split("}")[-1]}
    # Attributes
    if len(node.attrib):
        data["attributes"] = {
            key: val
            for key, val in node.attrib.items()
            if val != "unknown" and val != "0001-01-01T00:00:00"
        }
    # Children
    if len(node):
        for child in node:
            child_dict = etreeToDict(child)
            data[child_dict["key"]] = etreeToDict(child)
    else:
        data["val"] = node.text if node.text != "unknown" and node.text != "0001-01-01T00:00:00" else None
    # Return
    return data


def getMetadataFromConsensusReadset(in_readset_info_path):
    """
    Return run information tree from consensusreadset.xml.

    :param in_readset_info_path: Path to consensus readset file.
    :type in_readset_info_path: str
    :return: Return run information tree ffrom consensusreadset.xml.
    :rtype: dict
    """
    root = ET.parse(in_readset_info_path).getroot()
    return etreeToDict(root)


def getRunInfo(in_smart_folder):
    """
    Return run information (StartDate, EndDate, Operator, InstrumentID, ...) from smart cell folder.

    :param in_smart_folder: Path to smart cell folder.
    :type in_smart_folder: str
    :return: Return run information (StartDate, EndDate, Operator, InstrumentID, ...) from smart cell folder.
    :rtype: dict
    """
    readset_info_path = glob.glob(os.path.join(in_smart_folder, "*.consensusreadset.xml"))[0]
    readset_info = getMetadataFromConsensusReadset(readset_info_path)
    collection_metadata = readset_info["DataSetMetadata"]["Collections"]["CollectionMetadata"]
    run_dict = dict()
    run_dict["InstrumentID"] = collection_metadata["attributes"]["InstrumentId"]
    # Run
    run_dict["RunID"] = collection_metadata["attributes"]["Context"]
    run_details = collection_metadata["RunDetails"]
    if run_details.get("Name") and run_details["Name"]["val"]:
        run_dict["RunName"] = run_details["Name"]["val"]
    if run_details.get("CreatedBy") and run_details["CreatedBy"]["val"]:
        run_dict["CreatedBy"] = run_details["CreatedBy"]["val"]
    if run_details.get("CreatedBy") and run_details["StartedBy"]["val"]:
        run_dict["OperatorName"] = run_details["StartedBy"]["val"]
    # Run dates
    run_dict["StartDate"] = _consReadsetDateToDatetime(collection_metadata["attributes"]["CreatedAt"])
    run_dict["EndDate"] = None
    if readset_info["ExternalResources"].get("ExternalResource"):
        run_dict["EndDate"] = _consReadsetDateToDatetime(
            readset_info["ExternalResources"]["ExternalResource"]["attributes"]["CreatedAt"]
        )
    # Well
    well_info = collection_metadata["WellSample"]
    if well_info.get("WellName") and well_info["WellName"]["val"]:
        run_dict["WellName"] = well_info["WellName"]["val"]
    if well_info.get("Application") and well_info["Application"]["val"]:
        run_dict["Application"] = well_info["Application"]["val"]
    if well_info.get("IsoSeq") and well_info["IsoSeq"]["val"] is not None:
        run_dict["IsoSeq"] = well_info["IsoSeq"]["val"] == "true"
    if well_info.get("IsoSeq") and well_info["IsCCS"]["val"] is not None:
        run_dict["IsCCS"] = well_info["IsCCS"]["val"] == "true"
    return run_dict


def getRunNameFromSequelHifi(in_smart_folder):
    """
    Return run name from SMRT Link file in run folder.

    :param in_smart_folder: Path to run folder.
    :type in_smart_folder: str
    :return: Run name
    :rtype: str
    """
    return getRunInfo(in_smart_folder).get("RunName")


class Run:
    """Read and provide getters on Run from run folder."""

    def __init__(self, path):
        """
        Build and return an instance of Run.

        :param path: Path to the run folder.
        :type path: str
        :return: The new instance.
        :rtype: Run
        """
        self.path = path

    def isCopied(self):
        """
        Return True if copy from temp to output folder is ended by instrument.

        :return: True if copy from temp to output folder is ended by instrument.
        :rtype: boolean
        """
        is_completed = False
        if self.isSequenced():
            copy_complete_file = os.path.join(self.path, "*.transferdone")
            is_completed = len(glob.glob(copy_complete_file)) != 0
        return is_completed

    def isEnded(self):
        """
        Return True if all instrument steps are ended (sequencing and copy).

        :return: True if all instrument steps are ended.
        :rtype: boolean
        """
        return self.isSequenced() and self.isCopied()

    def isSequenced(self):
        """
        Return True if sequencing step is ended on instrument.

        :return: True if sequencing step is ended on instrument.
        :rtype: boolean
        """
        run_complete_file = os.path.join(self.path, "*.zmw_metrics.json.gz")
        return len(glob.glob(run_complete_file)) != 0
