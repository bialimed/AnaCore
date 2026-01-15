# -*- coding: utf-8 -*-
"""Classes to read bases2fastq stats files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

from anacore.instrument.demultiplex import AbstractDemultStat
import csv
import datetime
import re


class DemultLog(object):
    """Reader for bcl-convert log."""

    def __init__(self, path):
        """
        Build and return an instance of DemultLog.

        :param path: Path to bases2fastq log file (format: txt).
        :type path: str
        :return: The new instance.
        :rtype: DemultLog
        """
        self.filepath = path
        self.command = None
        self.end_time = None
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
        with open(self.filepath) as reader:
            for line in reader:
                line = line.strip()
                matches = re.match(r"^\<(.+)\> \[(.+)\]: (.+)$", line)
                if matches:
                    log_time, log_lvl, msg, = matches.groups()
                    if msg.startswith("bases2fastq version:"):
                        self.start_time = datetime.datetime.strptime(log_time, '%Y-%m-%d %H:%M:%S')
                        self.version = re.match(r"^bases2fastq version: ([^, ]+)", msg).groups()[0]
                    elif msg.startswith("Run command:"):
                        self.command = msg.replace("Run command: ", "")
                    elif log_lvl == "error":
                        self.status = "FAILED"
                        self.end_time = datetime.datetime.strptime(log_time, '%Y-%m-%d %H:%M:%S')
                    elif msg.startswith("Output stored"):
                        self.end_time = datetime.datetime.strptime(log_time, '%Y-%m-%d %H:%M:%S')
                        self.status = "COMPLETED"


class DemultStat(AbstractDemultStat):
    """Reader for demultipexing statistics file (IndexAssignment.csv and UnassignedSequences.csv) from bases2fastq."""

    def __init__(self, demult_path, undet_path=None):
        """
        Build and return an instance of DemultStat from bases2fastq.

        :param demult_path: Path to the samples demultiplexing statistics file.
        :type demult_path: str
        :param undeter_path: Path to the unassigned demultiplexing statistics file.
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
                if "+" not in row["Lane"]:
                    # Sample
                    spl_id = row["SampleName"]
                    if spl_id not in tmp_data:
                        tmp_data[spl_id] = {
                            "id": spl_id,
                            "barcodes": dict()
                        }
                    spl = tmp_data[spl_id]
                    # Barcode
                    idx_seq = row["I1"]
                    if "I2" in row:
                        idx_seq = "{}+{}".format(row["I1"], row["I2"])
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
                            "ct": [0]
                        }
                    lane = barcode["lanes"][lane_id]
                    lane["ct"][0] += int(row["NumPoloniesAssigned"])
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
        Return unassigned barcodes.

        :return: Undetermined barcodes. Format: [{"seq": "ATGG+TTC", "lanes": [{"id": 1, "ct": 154782}, {"id": 2, "ct": 255567]}].
        :rtype: list
        """
        if self._undet_path is None:
            raise Exception("UnassignedSequences.csv must be provide for access to undetermined.")
        barcode_by_seq = dict()
        with open(self._undet_path) as handle:
            reader = csv.DictReader(handle, delimiter=',')
            for row in reader:
                if "+" not in row["Lane"]:
                    idx_seq = row["I1"]
                    if "I2" in row:
                        idx_seq = "{}+{}".format(row["I1"], row["I2"])
                    if idx_seq not in barcode_by_seq:
                        barcode_by_seq[idx_seq] = {
                            "seq": idx_seq,
                            "lanes": list()
                        }
                    barcode_by_seq[idx_seq]["lanes"].append({
                        "id": row["Lane"],
                        "ct": int(row["Count"])
                    })
        return barcode_by_seq.values()
