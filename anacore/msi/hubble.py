# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing Hubble's outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sequenceIO import IdxFastaIO
from anacore.sv import HashedSVIO
from anacore.msi.base import Status
from anacore.msi.locus import Locus
from anacore.msi.sample import MSISample, MSISplRes
import os
import json
from textwrap import wrap


class HubbleDiff(HashedSVIO):
    """Manage differential analysis results file coming from Hubble."""

    def __init__(self, filepath, mode="r", ref_path=None):
        """
        Return the new instance of HubbleDiff.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param ref_path: Path to reference sequences file used in analysis (format: fasta with index).
        :type ref_path: str
        :return: New instance of HubbleDiff.
        :rtype: HubbleDiff
        """
        super().__init__(filepath, mode, "\t", "#")
        if mode == "w":
            self.titles = ["Chromosome", "Start", "RepeatUnit", "Assessed", "Distance", "PValue"]
        elif mode == "r":
            if ref_path is None:
                raise Exception("In read mode {} require reference sequences path.".format(self.__class__.__name__))
            self.ref = IdxFastaIO(ref_path, "r")

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line.
        :rtype: anacore.msi.base.Locus
        """
        record = super()._parseLine()
        microsat_length = nbRepeatFromStart(
            self.ref, record["Chromosome"], int(record["Start"]), record["RepeatUnit"]
        ) * len(record["RepeatUnit"])
        return Locus.fromDict({
            "position": "{}:{}-{}".format(
                record["Chromosome"],
                record["Start"],
                int(record["Start"]) + microsat_length
            ),
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
            "Chromosome": record.position.rsplit(":", 1)[0],
            "Start": int(record.position.rsplit(":", 1)[1].split("-")[0]),
            "RepeatUnit": res.data["nt"],
            "Assessed": "False" if res.data["distance"] is None else "True",
            "Distance": "" if res.data["distance"] is None else res.data["distance"],
            "PValue": "" if res.data["p_value"] is None else res.data["p_value"]
        }
        return super().recordToLine(formatted_record)


class HubbleDist(HashedSVIO):
    """Manage distributions file coming from hubble."""

    def __init__(self, filepath, mode="r", ref_path=None, len_limit=100):
        """
        Return the new instance of HubbleDist.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param ref_path: Path to reference sequences file used in analysis (format: fasta with index).
        :type ref_path: str
        :param len_limit: Max length recorded in length_distribution.
        :type len_limit: int
        :return: New instance of HubbleDist.
        :rtype: HubbleDist
        """
        super().__init__(filepath, mode, "\t", "#")
        self.len_limit = len_limit
        if mode == "w":
            self.titles = ["chromosome", "location", "repeat_unit_bases", "covered", "length_distribution"]
        elif mode == "r":
            if ref_path is None:
                raise Exception("In read mode {} require reference sequences path.".format(self.__class__.__name__))
            self.ref = IdxFastaIO(ref_path, "r")

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
        # Get locus reference length
        microsat_length = nbRepeatFromStart(
            self.ref, record["chromosome"], int(record["location"]), record["repeat_unit_bases"]
        ) * len(record["repeat_unit_bases"])
        # Locus
        return Locus.fromDict({
            "position": "{}:{}-{}".format(
                record["chromosome"],
                record["location"],
                int(record["location"]) + microsat_length
            ),
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
            "chromosome": record.position.rsplit(":", 1)[0],
            "location": int(record.position.rsplit(":", 1)[1].split("-")[0]),
            "repeat_unit_bases": res.data["nt"],
            "covered": "False" if res.data["lengths"].getCount() == 0 else "True",
            "length_distribution": ",".join(
                map(str, res.data["lengths"].getDenseCount(1, self.len_limit))  # Lengths count start from 1
            )
        }
        return super().recordToLine(formatted_record)


def nbRepeatFromStart(ref_reader, chrom, start, motif):
    """
    Return number of occurrences of the pattern at the beginning of the sequence.

    :param ref_reader: File handle to reference sequences file.
    :type ref_reader: anacore.sequenceIO.IdxFastaIO
    :param chrom: Name of reference region.
    :type chrom: str
    :param start: Start position on reference (0-based).
    :type start: int
    :param motif: Repeat motif base.
    :type motif: str
    :return: Number of occurrences of the pattern at the beginning of the sequence.
    :rtype: int
    """
    nb = 0
    motif = motif.upper()
    read_chunk_size = len(motif) * 30
    end_search = False
    while not end_search:
        sub_seq = ref_reader.getSub(
            chrom, start + 1, start + read_chunk_size
        )
        if sub_seq == "":  # Chromosome is ended
            end_search = True
        else:
            sub_seq = sub_seq.upper()
            for sub_seq_win in wrap(sub_seq, len(motif)):
                if sub_seq_win != motif:
                    end_search = True
                    break
                nb += 1
            start += read_chunk_size
    return nb


def parseHubbleResults(sequences_path, summary_path, differential_path, distributions_path=None, unstable_rate_threshold=0.2):
    """
    Return information on loci and analyses about microsatellites for sample.

    :param sequences_path: Path to reference sequences file used in analysis (format: fasta with index).
    :type sequences_path: str
    :param summary_path: Path to microsatellites analysis summary file (format: JSON).
    :type summary_path: str
    :param differential_path: Path to microsatellites diffrential analysis results file (format: JSON).
    :type differential_path: str
    :param distributions_path: Path to microsatellites sizes distributions file (format: JSON).
    :type distributions_path: str
    :param unstable_rate_threshold: Over this rate of unstable loci the sample is tagged as unstable.
    :type unstable_rate_threshold: float
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
        if unstable_rate_threshold:
            if sum_data["PercentageUnstableSites"] > unstable_rate_threshold * 100:
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
    for curr_locus in HubbleDiff(differential_path, ref_path=sequences_path):
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
        for curr_locus in HubbleDist(distributions_path, ref_path=sequences_path):
            old_res = msi_spl.loci[curr_locus.position].results["Hubble"]
            curr_locus.results["Hubble"].status = old_res.status
            curr_locus.results["Hubble"].score = old_res.score
            curr_locus.results["Hubble"].data["distance"] = old_res.data["distance"]
            curr_locus.results["Hubble"].data["p_value"] = old_res.data["p_value"]
            msi_spl.loci[curr_locus.position] = curr_locus
    return msi_spl
