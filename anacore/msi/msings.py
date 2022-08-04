# -*- coding: utf-8 -*-
"""Classes and functions to predict with mSINGS's algorithm and to read/write mSINGS binaries outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.5.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
from anacore.abstractFile import isGzip
from anacore.msi.base import Status
from anacore.msi.locus import Locus, LocusRes
from anacore.msi.sample import MSISample, MSISplRes
from anacore.sv import HashedSVIO
from numpy import average, std


class MSINGSAnalysisIO(HashedSVIO):
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
        :rtype: anacore.msi.base.Locus
        """
        record = super()._parseLine()
        locus_len = int(record["Position"].rsplit("-", 1)[1]) - int(record["Position"].rsplit(":", 1)[1].rsplit("-", 1)[0])
        peaks = record["IndelLength:AlleleFraction:SupportingCalls"].split(" ")
        nb_by_length = {}
        for idx, curr_peak in enumerate(peaks):
            indel_length, AF, DP = curr_peak.split(":")
            if DP != "0":
                nb_by_length[int(indel_length) + locus_len] = int(DP)
        return Locus.fromDict({
            "position": record["Position"],
            "name": record["Name"],
            "results": {
                "mSINGS": {
                    "status": Status.undetermined,
                    "data": {
                        "avg_depth": int(record["Average_Depth"]),
                        "nb_peaks": int(record["Number_of_Peaks"]),
                        "std_dev": float(record["Standard_Deviation"]),
                        "lengths": {"ct_by_len": nb_by_length}
                    }
                }
            }
        })

    def recordToLine(self, record):
        """
        Return the record in SV format.

        :param record: The record containing a result coming from mSINGS.
        :type record: anacore.msi.base.Locus
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        formatted_record = {
            "Position": record.position,
            "Name": record.name,
            "Average_Depth": record.results["mSINGS"].data["avg_depth"],
            "Number_of_Peaks": record.results["mSINGS"].data["nb_peaks"],
            "Standard_Deviation": record.results["mSINGS"].data["std_dev"]
        }
        ct_by_len = record.results["mSINGS"].data["lengths"]
        max_ct = ct_by_len.getMostRepresented()["count"]
        if max_ct == 0:
            formatted_record["IndelLength:AlleleFraction:SupportingCalls"] = "0:0.0000000:0"
        else:
            locus_ref_len = record.length
            locus_min_len = ct_by_len.getMinLength()
            formatted_record["IndelLength:AlleleFraction:SupportingCalls"] = " ".join(
                [
                    "{}:{:.7f}:{}".format((locus_min_len + idx) - locus_ref_len, ct / max_ct, ct)
                    for idx, ct in enumerate(ct_by_len.getDenseCount())
                ]
            )
        return super().recordToLine(formatted_record)


class MSINGSEval:
    """Provides utils to predict with mSINGS: calculate feature and determine feature threshold."""

    @staticmethod
    def getNbPeaks(locus_distrib, highest_peak_ratio=0.05):
        """
        Return feature used to predict stability status by mSINGS: the number of peaks with height >= peak_height_cutoff * highest_peak_height.

        :param locus_distrib: Lengths distribution of microsatellite.
        :type locus_distrib: anacore.msi.locus.LocusDataDistrib
        :param highest_peak_ratio: Minimum rate of highest peak height to consider a peak.
        :type highest_peak_ratio: float
        :return: Feature used to predict stability status by mSINGS: the number of peaks with height >= peak_height_cutoff * highest_peak_height.
        :rtype: int
        """
        highest_peak = locus_distrib.getMostRepresented()
        min_peak_height = highest_peak_ratio * highest_peak["count"]
        return locus_distrib.getNbPeaks(min_peak_height)

    @staticmethod
    def getThresholdFromNbPeaks(models_nb_peaks, std_dev_rate=2.0):
        """
        Return minimum number of peaks to tag distribution as unstable.

        :param models_nb_peaks: Number of peaks in stable population (see MSINGSEval.getNbPeaks).
        :type models_nb_peaks: list
        :param std_dev_rate: Number of standard deviation authorized around average number of peaks in stable population.
        :type std_dev_rate: float
        :return: Minimum number of peaks to tag distribution as unstable.
        :rtype: float
        """
        return average(models_nb_peaks) + std(models_nb_peaks) * std_dev_rate


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
        self.method_name = "mSINGS"
        self._parse()

    def _parse(self):
        """Parse file and store information in self."""
        if isGzip(self.filepath):
            with gzip.open(self.filepath, "rt") as FH:
                self._parseFileHandle(FH)
        else:
            with open(self.filepath, "r") as FH:
                self._parseFileHandle(FH)

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
        if not line.startswith("msings_score"):
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
                elif curr_val == "1" or curr_val == "Unstable":
                    loci_res = LocusRes(Status.unstable)
                elif curr_val == "0" or curr_val == "Stable":
                    loci_res = LocusRes(Status.stable)
                else:
                    raise Exception("{} cannot be parsed by {}.".format(self.filepath, self.__class__.__name__,))
                self.samples[curr_spl].addLocus(
                    Locus(curr_locus, None, {self.method_name: loci_res})
                )
