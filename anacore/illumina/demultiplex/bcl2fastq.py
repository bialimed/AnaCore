# -*- coding: utf-8 -*-
"""Classes to read bcl2fastq logs and stats files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.illumina.demultiplex.base import AbstractDemultStat
import datetime
import json
import re


class DemultLog(object):
    """Reader for bcl2fastq log."""

    def __init__(self, path):
        """
        Build and return an instance of DemultLog.

        :param path: Path to BCL2fastq log file (format: txt).
        :type path: str
        :return: The new instance.
        :rtype: DemultLog
        """
        self.command = None
        self.end_time = None
        self.filepath = path
        self.start_time = None
        self.status = None  # COMPLETED or FAILED or RUNNING
        self.version = None
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        self.command = None
        self.end_time = None
        self.start_time = None
        self.status = "RUNNING"
        self.version = None
        last_info = ""
        with open(self.filepath) as reader:
            for line in reader:
                line = line.strip()
                if " ERROR" in line:
                    self.end_time = datetime.datetime.strptime(
                        " ".join(line.split(" ", 3)[:2]),
                        '%Y-%m-%d %H:%M:%S'
                    )
                    self.status = "FAILED"
                elif self.version is None and line.startswith("bcl2fastq v"):
                    self.version = line.split(" ")[1]
                elif self.command is None and "Command-line invocation:" in line:
                    trace_time, command = line.split(" Command-line invocation: ")
                    self.command = command
                    self.start_time = datetime.datetime.strptime(
                        trace_time.rsplit(" ", 1)[0],
                        '%Y-%m-%d %H:%M:%S'
                    )
                elif line != "":
                    last_info = line
        if "Processing completed" in last_info:
            self.end_time = datetime.datetime.strptime(
                " ".join(last_info.split(" ", 3)[:2]),
                '%Y-%m-%d %H:%M:%S'
            )
            nb_errors = int(
                re.search(r"(\d+) errors", last_info).group(1)
            )
            if nb_errors != 0:
                self.status = "FAILED"
            else:
                self.status = "COMPLETED"


class DemultStat(AbstractDemultStat):
    """Reader for demultipexing statistics file (Data/Intensities/BaseCalls/Stats/Stats.json) from bcl2fastq."""

    def _parse(self):
        """Read self._path content and store information in self.samples."""
        with open(self._path) as reader:
            in_data = json.load(reader)["ConversionResults"]
        tmp_data = dict()
        for in_lane in in_data:
            # Sample
            for in_spl in in_lane["DemuxResults"]:
                spl_id = in_spl["SampleId"]
                if spl_id not in tmp_data:
                    tmp_data[spl_id] = {
                        "id": spl_id,
                        "barcodes": dict()
                    }
                spl = tmp_data[spl_id]
                # Barcode
                for in_barcode in in_spl["IndexMetrics"]:
                    idx_seq = in_barcode["IndexSequence"]
                    if idx_seq not in spl["barcodes"]:
                        spl["barcodes"][idx_seq] = {
                            "seq": idx_seq,
                            "lanes": dict()
                        }
                    barcode = spl["barcodes"][idx_seq]
                    # Lane
                    lane_id = in_lane["LaneNumber"]
                    if lane_id not in barcode["lanes"]:
                        barcode["lanes"][lane_id] = {
                            "id": lane_id,
                            "ct": [0, 0, 0]
                        }
                    lane = barcode["lanes"][lane_id]
                    for nb_mis, ct in in_barcode["MismatchCounts"].items():
                        lane["ct"][int(nb_mis)] += ct
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
        with open(self._path) as reader:
            in_data = json.load(reader)["UnknownBarcodes"]
        barcode_by_seq = dict()
        for in_lane in in_data:
            for idx_seq, count in in_lane["Barcodes"].items():
                if idx_seq not in barcode_by_seq:
                    barcode_by_seq[idx_seq] = {
                        "seq": idx_seq,
                        "lanes": list()
                    }
                barcode_by_seq[idx_seq]["lanes"].append({
                    "id": in_lane["Lane"],
                    "ct": count
                })
        return barcode_by_seq.values()
