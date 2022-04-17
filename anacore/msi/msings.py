# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing mSINGS's outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
from anacore.abstractFile import isGzip
from anacore.sv import HashedSVIO
from anacore.msi import MSILocus, LocusRes, MSISample, MSISplRes, Status


class MSINGSAnalysis(HashedSVIO):
    """Manage output file produced by the command "msi analyzer" of mSINGS (https://bitbucket.org/uwlabmed/msings)."""

    def __init__(self, filepath, mode="r"):
        """
        Return the new instance of MSINGSAnalysis.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        """
        super().__init__(filepath, mode, "\t")
        if mode == "w":
            self.titles = ["Position", "Name", "Average_Depth", "Number_of_Peaks", "Standard_Deviation", "IndelLength:AlleleFraction:SupportingCalls"]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line.
        :rtype: MSILocus
        """
        record = super()._parseLine()
        peaks = record["IndelLength:AlleleFraction:SupportingCalls"].split(" ")
        nb_by_length = {}
        if len(peaks) == 1 and peaks[0] == "0:0.0:0":
            peaks = []
        else:
            for idx, curr_peak in enumerate(peaks):
                indel_length, AF, DP = curr_peak.split(":")
                peaks[idx] = {
                    "indel_length": int(indel_length),
                    "AF": float(AF),
                    "DP": int(DP)
                }
                nb_by_length[int(indel_length)] = int(DP)

        return MSILocus.fromDict({
            "position": record["Position"],
            "name": record["Name"],
            "results": {
                "MSINGS": {
                    "_class": "LocusResDistrib",
                    "status": Status.undetermined,
                    "data": {
                        "avg_depth": record["Average_Depth"],
                        "nb_peaks": record["Number_of_Peaks"],
                        "peaks": peaks,
                        "std_dev": record["Standard_Deviation"],
                        "nb_by_length": nb_by_length
                    }
                }
            }
        })

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: The record containing a result coming from mSINGS.
        :type record: MSILocus
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        formatted_record = {
            "Position": record.position,
            "Name": record.name,
            "Average_Depth": record.results["MSINGS"].data["avg_depth"],
            "Number_of_Peaks": record.results["MSINGS"].data["nb_peaks"],
            "Standard_Deviation": record.results["MSINGS"].data["std_dev"]
        }
        formatted_record["IndelLength:AlleleFraction:SupportingCalls"] = " ".join(
            ["{}:{}:{}".format(curr_peak["indel_length"], curr_peak["AF"], curr_peak["DP"]) for curr_peak in record.results["MSINGS"].data["peaks"]]
        )
        return super().recordToLine(formatted_record)


class MSINGSReport(object):
    """Manage output file produced by the command "msi count_msi_samples" of mSINGS (https://bitbucket.org/uwlabmed/msings)."""

    def __init__(self, filepath):
        """
        Return the new instance of MSINGSReport.

        :param filepath: The filepath.
        :type filepath: str
        """
        self.filepath = filepath
        self.samples = dict()
        self.loci = list()
        self.method_name = "MSINGS"
        self.parse()

    def _parseFileHandle(self, FH):
        """
        Parse file referenced by the file handle FH and store information in self.

        :param FH: The file handle for the filepath.
        :type FH: TextIOWrapper
        """
        # Parse general information
        samples = [elt.strip() for elt in FH.readline().split("\t")[1:]]
        nb_unstable = [int(elt.strip()) for elt in FH.readline().split("\t")[1:]]
        nb_evaluated = [int(elt.strip()) for elt in FH.readline().split("\t")[1:]]
        line = FH.readline()
        if not line.startswith("msing_score"):
            scores = [None for curr_spl in samples]
        else:
            scores = [elt.strip() for elt in line.split("\t")][1:]
            for idx, elt in enumerate(scores):
                curr_score = None
                if elt != "":
                    curr_score = float(elt)
                scores[idx] = curr_score
            line = FH.readline()  # To next line
        status = [elt.strip() for elt in line.split("\t")][1:]
        for idx, elt in enumerate(status):
            curr_status = Status.undetermined
            if nb_evaluated[idx] > 0:
                if elt == "NEG":
                    curr_status = Status.stable
                elif elt == "POS":
                    curr_status = Status.unstable
            status[idx] = curr_status
        for spl_idx, curr_spl in enumerate(samples):
            if status[spl_idx] == Status.stable:
                """
                Convert mSINGS score (ratio unstable / determined) to a confidence score.
                Nb loci unstables:    0    1    2    3    4
                mSINGS score:         0 0.25 0.50 0.75    1
                confi unstable:       0 0.25 0.50 0.75    1
                confi stable:         1 0.75 0.50 0.25    0
                """
                scores[spl_idx] = 1 - scores[spl_idx]
            spl_res = MSISplRes(status[spl_idx], scores[spl_idx], self.method_name)
            self.samples[curr_spl] = MSISample(curr_spl, None, {self.method_name: spl_res})
        # Parse loci information
        for curr_line in FH:
            fields = [elt.strip() for elt in curr_line.split("\t")]
            curr_locus = fields[0]
            self.loci.append(curr_locus)
            for idx, curr_val in enumerate(fields[1:]):
                curr_spl = samples[idx]
                loci_res = None
                if curr_val == "":
                    loci_res = LocusRes(Status.undetermined)
                elif curr_val == "1":
                    loci_res = LocusRes(Status.unstable)
                else:
                    loci_res = LocusRes(Status.stable)
                self.samples[curr_spl].addLocus(
                    MSILocus(curr_locus, None, {self.method_name: loci_res})
                )

    def parse(self):
        """Parse file and store information in self."""
        if isGzip(self.filepath):
            with gzip.open(self.filepath, "rt") as FH:
                self._parseFileHandle(FH)
        else:
            with open(self.filepath, "r") as FH:
                self._parseFileHandle(FH)
