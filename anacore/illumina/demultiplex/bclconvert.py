# -*- coding: utf-8 -*-
"""Classes to read bcl-convert logs and stats files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.illumina.demultiplex.base import AbstractDemultStat
import csv
import datetime
import os


class DemultLog(object):
    """Reader for bcl-convert log."""

    def __init__(self, folder_path):
        """
        Build and return an instance of DemultLog.

        :param folder_path: Path to bcl-convert log folder (${out_folder}/Logs).
        :type folder_path: str
        :return: The new instance.
        :rtype: DemultLog
        """
        self.command = None
        self.end_time = None
        self.folder_path = folder_path
        self.start_time = None
        self.status = None  # COMPLETED or FAILED or RUNNING
        self.version = None
        self._parse()

    def _parse(self):
        """Read logs and store information on the instance's attributes."""
        self.command = None
        self.end_time = None
        self.start_time = None
        self.status = "RUNNING"
        self.version = None
        self._parseInfo()
        self._parseError()

    def _parseError(self):
        """Read error log file and store information on the instance's attributes."""
        error_path = os.path.join(self.folder_path, "Errors.log")
        if os.path.getsize(error_path) > 0:
            self.status = "FAILED"
            with open(error_path) as reader:
                line = reader.readline()
                date, trash, thread, info = line.strip().split(" ", 3)
                self.end_time = datetime.datetime.fromisoformat(date.replace('Z', '+00:00'))
            if self.start_time is None:
                self.start_time = self.end_time

    def _parseInfo(self):
        """Read info log file and store information on the instance's attributes."""
        last = (None, None)  # (date, info)
        info_path = os.path.join(self.folder_path, "Info.log")
        with open(info_path) as reader:
            for line in reader:
                if line.strip().count(" ") > 2:  # info message is not empty
                    date, trash, thread, info = line.strip().split(" ", 3)
                    if self.start_time is None:
                        self.start_time = datetime.datetime.fromisoformat(date.replace('Z', '+00:00'))
                    if self.version is None and info.startswith("bcl-convert Version "):
                        self.version = info.rsplit(" ", 1)[1].strip()
                        self.version = self.version.replace("00.000.000.", "")
                    elif self.command is None and info.startswith("Command Line: "):
                        arguments = info.split(":", 1)[1].strip()
                        self.command = "bcl-convert {}".format(arguments)
                    last = (date, info)
        if last[1] is not None and last[1].startswith("Conversion Complete."):
            self.end_time = datetime.datetime.fromisoformat(last[0].replace('Z', '+00:00'))
            self.status = "COMPLETED"


class DemultStat(AbstractDemultStat):
    """Reader for demultipexing statistics file (Demultiplex_Stats.csv and Top_Unknown_Barcodes.csv) from bcl-convert."""

    def __init__(self, demult_path, undet_path=None):
        """
        Build and return an instance of DemultStat from bcl-convert.

        :param demult_path: Path to the samples demultiplexing statistics file.
        :type demult_path: str
        :param undeter_path: Path to the undetermined demultiplexing statistics file.
        :type undeter_path: str
        :return: The new instance.
        :rtype: DemultStat
        """
        super().__init__(demult_path)
        self._undet_path = undet_path

    def _parse(self):
        """Read self._path content and store information in self.samples."""
        tmp_data = dict()
        with open(self._path) as handle:
            reader = csv.DictReader(handle, delimiter=',')
            for row in reader:
                # Sample
                spl_id = row["SampleID"]
                if spl_id not in tmp_data:
                    tmp_data[spl_id] = {
                        "id": spl_id,
                        "barcodes": dict()
                    }
                spl = tmp_data[spl_id]
                # Barcode
                idx_seq = row["Index"].replace("-", "+")
                if idx_seq not in spl["barcodes"]:
                    spl["barcodes"][idx_seq] = {
                        "seq": idx_seq,
                        "lanes": dict()
                    }
                barcode = spl["barcodes"][idx_seq]
                # Lane
                lane_id = int(row["Lane"])
                if lane_id not in barcode["lanes"]:
                    barcode["lanes"][lane_id] = {
                        "id": lane_id,
                        "ct": [0, 0, 0]
                    }
                lane = barcode["lanes"][lane_id]
                lane["ct"][0] += int(row["# Perfect Index Reads"])
                lane["ct"][1] += int(row["# One Mismatch Index Reads"])
                lane["ct"][2] += int(row["# Two Mismatch Index Reads"])
        # Remove dict only necessary for parsing
        self.samples = list()
        for spl_id, spl in tmp_data.items():
            spl["barcodes"] = spl["barcodes"].values()
            for barcode in spl["barcodes"]:
                barcode["lanes"] = barcode["lanes"].values()
            self.samples.append(spl)

    @property
    def undetermined(self):
        """
        Return undetermined barcodes.

        :return: Undetermined barcodes. Format: [{"seq": "ATGG+TTC", "lanes": [{"id": 1, "ct": 154782}, {"id": 2, "ct": 255567]}].
        :rtype: list
        """
        if self._undet_path is None:
            raise Exception("Top_Unknown_Barcodes.csv must be provide for access to undetermined.")
        barcode_by_seq = dict()
        with open(self._undet_path) as handle:
            reader = csv.DictReader(handle, delimiter=',')
            for row in reader:
                idx_seq = row["index"]
                if "index2" in row:
                    idx_seq = "{}+{}".format(row["index"], row["index2"])
                if idx_seq not in barcode_by_seq:
                    barcode_by_seq[idx_seq] = {
                        "seq": idx_seq,
                        "lanes": list()
                    }
                barcode_by_seq[idx_seq]["lanes"].append({
                    "id": row["Lane"],
                    "ct": int(row["# Reads"])
                })
        return barcode_by_seq.values()
