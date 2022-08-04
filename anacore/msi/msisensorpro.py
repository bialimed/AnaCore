# -*- coding: utf-8 -*-
"""Classes and functions to predict with MSISensor-pro's algorithm and to read/write MSIsensor-pro binaries outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
from anacore.abstractFile import AbstractFile, isGzip
from anacore.msi.base import Status
from anacore.msi.locus import getRefSeqInfo, Locus
from anacore.msi.sample import MSISample, MSISplRes
from anacore.sv import HashedSVIO
from numpy import mean, std


class BaselineIO(HashedSVIO):
    """Manage baseline file."""

    def __init__(self, filepath, mode="r"):
        """
        Return the new instance of BaselineIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: New instance of BaselineIO.
        :rtype: BaselineIO
        """
        super().__init__(filepath, mode, "\t")
        if mode == "w":
            self.titles = [
                "chromosome",
                "location",
                "repeat_unit_length",
                "repeat_unit_binary",
                "repeat_times",
                "left_flank_binary",
                "right_flank_binary",
                "repeat_unit_bases",
                "left_flank_bases",
                "right_flank_bases",
                "threshold",
                "supportSamples"
            ]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line.
        :rtype: anacore.msi.msisensor.BaselineRecord
        """
        hash = super()._parseLine()
        for field in ["location", "repeat_times", "supportSamples"]:
            hash[field] = int(hash[field])
        hash["threshold"] = float(hash["threshold"])
        hash["support_samples"] = hash["supportSamples"]
        for field in ["repeat_unit_length", "repeat_unit_binary", "left_flank_binary", "right_flank_binary", "supportSamples"]:
            del(hash[field])
        record = BaselineRecord(**hash)
        return record

    def recordToLine(self, record):
        """
        Return the record in baseline format.

        :param record: The record to process.
        :type record: BaselineRecord
        :return: The baseline line corresponding to the record.
        :rtype: str
        """
        hash = {field: getattr(record, field) for field in self.titles if field != "supportSamples"}
        hash["supportSamples"] = record.support_samples
        return super().recordToLine(hash)


class BaselineRecord:
    """Class to manage one baseline record."""

    def __init__(self, chromosome, location, repeat_times, repeat_unit_bases, left_flank_bases, right_flank_bases, threshold, support_samples):
        """
        Return the new instance of BaselineRecord.

        :param chromosome: Chromosome name.
        :type chromosome: str
        :param location: Start position of locus (0-based).
        :type location: int
        :param repeat_times: Number of repeat of repeat_unit_bases.
        :type repeat_times: int
        :param repeat_unit_bases: Sequence of repeat unit.
        :type repeat_unit_bases: str
        :param left_flank_bases: Sequence before locus.
        :type left_flank_bases: str
        :param right_flank_bases: Sequence after locus.
        :type right_flank_bases: str
        :param threshold: Minimum score to classify a locus as unstable.
        :type threshold: float
        :param support_samples: Number of samples used in threshold calculation.
        :type support_samples: int
        :return: New instance of BaselineRecord.
        :rtype: BaselineRecord
        """
        self.chromosome = chromosome
        self.location = location
        self.repeat_times = repeat_times
        self.repeat_unit_bases = repeat_unit_bases
        self.left_flank_bases = left_flank_bases
        self.right_flank_bases = right_flank_bases
        self.threshold = threshold
        self.support_samples = support_samples

    @property
    def left_flank_binary(self):
        """
        Return MSIsensor binary representation of sequence before locus.

        :return: MSIsensor binary representation of sequence before locus.
        :rtype: int
        """
        return stringToBinary(self.left_flank_bases)

    @property
    def repeat_unit_binary(self):
        """
        Return MSIsensor binary representation of repeat unit sequence.

        :return: MSIsensor binary representation of repeat unit sequence.
        :rtype: int
        """
        return stringToBinary(self.repeat_unit_bases)

    @property
    def repeat_unit_length(self):
        """
        Return length of repeat unit sequence.

        :return: Length of repeat unit sequence.
        :rtype: int
        """
        return len(self.repeat_unit_bases)

    @property
    def right_flank_binary(self):
        """
        Return MSIsensor binary representation of sequence after locus.

        :return: MSIsensor binary representation of sequence after locus.
        :rtype: int
        """
        return stringToBinary(self.right_flank_bases)

    @staticmethod
    def fromModel(ref_fh, target, models, model_method="model", flank_size=5):
        """
        Return BaselineRecord from model.

        :param ref_fh: File handle to reference sequences.
        :type ref_fh: anacore.sequenceIO.IdxFastaIO
        :param target: Microsatellite region.
        :type target: anacore.region.Region
        :param models: Samples in model.
        :type models: list of MSISample
        :param model_method: Name of method used to store model information.
        :type model_method: str
        :param flank_size: Size of microsatellite flanking sequence used in baseline.
        :type flank_size: int
        :return: BaselineRecord from model.
        :rtype: BaselineRecord
        """
        locus_id = "{}:{}-{}".format(target.reference.name, target.start - 1, target.end)
        locus_scan = getRefSeqInfo(ref_fh, target, flank_size)
        locus_scan["threshold"] = ProEval.getThreshold(
            models,
            locus_id,
            method_name=model_method
        )
        locus_scan["support_samples"] = len([
            spl for spl in models if spl.loci[locus_id].results[model_method].status == Status.stable
        ])
        return BaselineRecord(**locus_scan)


def binaryToString(binary, len):
    """
    Return sequence string from MSIsensor binary representation of sequence.

    :param binary: MSIsensor binary representation of sequence.
    :type binary: int
    :param len: Length of the original sequence.
    :type len: int
    :return: Sequence string from MSIsensor binary representation of sequence.
    :rtype: str
    """
    char_by_idx = {0: "A", 1: "C", 2: "G", 3: "T", 4: "N"}
    curr_binary = binary
    str = ""
    for idx in range(len):
        char_idx = curr_binary & 3
        str += char_by_idx[char_idx]
        curr_binary = curr_binary >> 2
    return str[::-1]


class DistIO(AbstractFile):
    """Manage distributions file."""

    def __init__(self, filepath, mode="r", len_limit=100):
        """
        Return the new instance of DistIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param len_limit: Max length recorded in length_distribution.
        :type len_limit: int
        :return: New instance of DistIO.
        :rtype: DistIO
        """
        super().__init__(filepath, mode)
        self.len_limit = len_limit

    def __iter__(self):
        is_end = False
        while not is_end:
            locus, meta = self.nextLocus()
            if locus is not None:
                yield(locus, meta)
            else:
                is_end = True

    def nextLocus(self):
        """
        Return the next locus.

        :return: The next locus object and locus metadata or None if it is the end of file.
        :rtype: anacore.msi.locus.Locus, dict
        """
        locus = None
        meta = None
        try:
            prev_file_pos = self.file_handle.tell()
            header = self.file_handle.readline().strip()
            new_file_pos = self.file_handle.tell()
            if prev_file_pos != new_file_pos:
                # header
                chrom, pos, left_flank, nb_repeat_and_nt, right_flank = header.split()
                nb_repeat, nt = nb_repeat_and_nt.split("[")
                nt = nt[:-1]
                self.current_line_nb += 1
                # lengths
                len_str = self.file_handle.readline().split(":")[1].strip()
                nb_by_length = {idx + 1: int(nb) for idx, nb in enumerate(len_str.split()) if nb != "0"}
                self.current_line_nb += 1
                # record
                meta = {
                    "chromosome": chrom,
                    "location": int(pos),
                    "left_flank_bases": left_flank,
                    "repeat_times": int(nb_repeat),
                    "repeat_unit_bases": nt,
                    "right_flank_bases": right_flank
                }
                locus = Locus.fromDict({
                    "position": "{}:{}-{}".format(chrom, pos, int(pos) + int(nb_repeat) * len(nt)),
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": Status.none,
                            "data": {
                                "lengths": {"ct_by_len": nb_by_length}
                            }
                        }
                    }
                })
        except Exception:
            raise IOError(
                'The line {} in "{}" cannot be parsed by {}.'.format(
                    self.current_line_nb,
                    self.filepath,
                    self.__class__.__name__
                )
            )
        return locus, meta

    def recordToLines(self, record, meta):
        """
        Return the 2 lines in MSIsensor-pro distribution format.

        :param record: Record containing distributions coming from MSIsensor-pro.
        :type record: anacore.msi.base.Locus
        :param meta: Metadata of locus: chromosome, location (0-based), left_flank_bases, repeat_times, repeat_unit_bases and right_flank_bases.
        :type meta: dict
        :return: The 2 lines corresponding to the record.
        :rtype: str
        """
        res = record.results["MSIsensor-pro_pro"]
        first_line = "{} {} {} {}[{}] {}".format(
            meta["chromosome"],
            meta["location"],
            meta["left_flank_bases"],
            meta["repeat_times"],
            meta["repeat_unit_bases"],
            meta["right_flank_bases"]
        )
        second_line = "T: {}".format(
            " ".join(
                map(str, res.data["lengths"].getDenseCount(1, self.len_limit))
            )
        )
        return "{}\n{}".format(first_line, second_line)

    def write(self, record, meta):
        """
        Write record lines in file.

        :param record: Record containing distributions coming from MSIsensor-pro.
        :type record: anacore.msi.base.Locus
        :param meta: Metadata of the locus: chromosome, location (0-based), left_flank_bases, repeat_times, repeat_unit_bases and right_flank_bases.
        :type meta: dict
        """
        self.file_handle.write(self.recordToLines(record, meta) + "\n")
        self.current_line_nb += 2


class ProEval:
    """Provides utils to predict with MSIsensor-pro pro: calculate feature and determine feature threshold."""

    @staticmethod
    def getSlippageScores(locus_distrib, ref_len):
        """
        Return features used to predict stability status by MSIsensor-pro: deletion score and insertion score.

        :param locus_distrib: Lengths distribution of microsatellite.
        :type locus_distrib: anacore.msi.locus.LocusDataDistrib
        :param ref_len: Length of the microsatellite in reference sequence.
        :type ref_len: int
        :return: features used to predict stability status by MSIsensor-pro: deletion score and insertion score.
        :rtype: tuple
        """
        distribution = locus_distrib.getDenseCount(
            1,
            max(ref_len, locus_distrib.getMaxLength())
        )
        del_score = 0
        normal_score = 0
        insert_score = 0
        for idx in range(0, ref_len - 1):
            del_score += distribution[idx] * (ref_len - (idx + 1))
            normal_score += distribution[idx] * (idx + 1)
        for idx in range(ref_len - 1, len(distribution)):
            normal_score += distribution[idx] * ref_len
            insert_score += distribution[idx] * (idx + 1 - ref_len)
        pro_p = del_score / (normal_score + del_score + insert_score)  # Deletion score
        pro_q = insert_score / (normal_score + del_score + insert_score)  # Insertion score
        return pro_p, pro_q

    @staticmethod
    def getThreshold(models, locus_id, method_name="model"):
        """
        Return the minimum score to classify locus as unstable.

        :param models: List of sample in model.
        :type models: list of MSISample
        :param locus_id: Selected ID.
        :type locus_id: str
        :param method_name: Name of the method used to store model information.
        :type method_name: str
        :return: Minimum score to classify locus as unstable.
        :rtype: float
        """
        scores = []
        for spl in models:
            if spl.loci[locus_id].results[method_name].status == Status.stable:
                scores.append(spl.loci[locus_id].results[method_name].data["MSIsensor-pro"]["pro_p"])
        return ProEval.getThresholdFromScores(scores)

    @staticmethod
    def getThresholdFromScores(models_scores, ddof=1):
        """
        Return the minimum score to classify locus as unstable. Due to floating point the result can differ from MSIsensor-pro.

        :param models_scores: List of scores for stable samples.
        :type models_scores: list
        :param ddof: Means Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents the number of elements.
        :type ddof: int
        :return: Minimum score to classify locus as unstable.
        :rtype: float
        """
        rounded_scores = [round(elt, 6) * 1000000 for elt in models_scores]
        return (mean(rounded_scores) + 3 * std(rounded_scores, ddof=ddof)) / 1000000


class ProIO(HashedSVIO):
    """Manage differential analysis results file coming from MSIsensor-pro pro."""

    def __init__(self, filepath, mode="r"):
        """
        Return the new instance of ProIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: New instance of ProIO.
        :rtype: ProIO
        """
        super().__init__(filepath, mode, "\t")
        if mode == "w":
            self.titles = [
                "chromosome", "location", "left_flank_bases", "repeat_times",
                "repeat_unit_bases", "right_flank_bases", "pro_p", "pro_q",
                "CovReads", "threshold"
            ]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record described by the current line (locus result and locus metadata).
        :rtype: (anacore.msi.base.Locus, dict)
        """
        record = super()._parseLine()
        locus_meta = {
            "chromosome": record["chromosome"],
            "location": int(record["location"]),
            "left_flank_bases": record["left_flank_bases"],
            "MSIsensor-pro": {"pro_p_threshold": float(record["threshold"])},
            "repeat_times": int(record["repeat_times"]),
            "repeat_unit_bases": record["repeat_unit_bases"],
            "right_flank_bases": record["right_flank_bases"]
        }
        locus = Locus.fromDict({
            "position": "{}:{}-{}".format(
                record["chromosome"],
                record["location"],
                int(record["location"]) + int(record["repeat_times"]) * len(record["repeat_unit_bases"])
            ),
            "results": {
                "MSIsensor-pro_pro": {
                    "data": {
                        "pro_p": float(record["pro_p"]),
                        "pro_q": float(record["pro_q"]),
                        "depth": int(record["CovReads"])
                    },
                    "score": None,
                    "status": Status.none
                }
            }
        })
        return locus, locus_meta

    @staticmethod
    def loci(filepath):
        """
        Return list of loci present in file.

        :param filepath: Path to the file (format: TSV).
        :type filepath: str
        :return: List of loci present in file.
        :rtype: list
        """
        loci = list()
        handler = open
        handler_options = "r"
        if isGzip(filepath):
            handler = gzip.open
            handler_options = "rt"
        with handler(filepath, handler_options) as reader:
            reader.readline()
            for line in reader:
                chrom, pos, left_flank_bases, repeat_times, repeat_unit_bases, trash = line.split(None, 5)
                loci.append(
                    "{}:{}-{}".format(chrom, pos, int(pos) + int(repeat_times) * len(repeat_unit_bases))
                )
        return loci

    def recordToLine(self, record, meta):
        """
        Return the record in SV format.

        :param record: Record containing all loci results coming from MSIsensor-pro pro.
        :type record: anacore.msi.locus.Locus
        :param meta: Metadata of the locus: chromosome, location (0-based), left_flank_bases, repeat_times, repeat_unit_bases, right_flank_bases and MSIsensor-pro.pro_p_threshold.
        :param meta: dict
        :return: The SV line corresponding to the record.
        :rtype: str
        """
        res = record.results["MSIsensor-pro_pro"]
        formatted_record = {
            "chromosome": record.position.split(":")[0],
            "location": record.position.split(":")[1].split("-")[0],
            "left_flank_bases": meta["left_flank_bases"],
            "repeat_times": meta["repeat_times"],
            "repeat_unit_bases": meta["repeat_unit_bases"],
            "right_flank_bases": meta["right_flank_bases"],
            "pro_p": "{:.6f}".format(res.data["pro_p"]),
            "pro_q": "{:.6f}".format(res.data["pro_q"]),
            "CovReads": res.data["depth"],
            "threshold": "{:.6f}".format(meta["MSIsensor-pro"]["pro_p_threshold"])
        }
        return super().recordToLine(formatted_record)

    def write(self, record, meta):
        """
        Write record lines in file.

        :param record: Record containing distributions coming from MSIsensor-pro.
        :type record: anacore.msi.base.Locus
        :param meta: Metadata of the locus: chromosome, location (0-based), left_flank_bases, repeat_times, repeat_unit_bases, right_flank_bases and MSIsensor-pro.pro_p_threshold.
        :type meta: dict
        """
        # Write header if the file is empty
        if self.current_line_nb == 0:
            self.writeHeader()
        # Write record
        self.file_handle.write(self.recordToLine(record, meta) + "\n")
        self.current_line_nb += 1


def parseProResults(spl_name, all_path, distributions_path=None, unstable_rate_threshold=0.1, determined_rate_threshold=0.5, min_support=20):
    """
    Return information on loci and analyses about microsatellites for sample.

    :param spl_name: Sample name.
    :type spl_name: str
    :param all_path: Path to sample analysis rsult for all loci (format: MSIsensor-pro pro all).
    :type all_path: str
    :param distributions_path: Path to microsatellites sizes distributions file (format: MSIsensor-pro pro dis).
    :type distributions_path: str
    :param unstable_rate_threshold: Over this rate of unstable loci the sample is tagged as unstable.
    :type unstable_rate_threshold: float
    :param determined_rate_threshold: Under this rate of loci (determined / all) the sample status cannot be determined.
    :type determined_rate_threshold: float
    :param min_support: Minimum number of reads on a locus to determine ist status.
    :type min_support: int
    :return: Information on loci and analyses about microsatellites for sample.
    :rtype: anacore.msi.MSISpl
    """
    method = "MSIsensor-pro_pro"
    spl = MSISample(spl_name)
    # Parse analysis results
    for locus, meta in ProIO(all_path):
        res_locus = locus.results[method]
        res_locus.status = Status.undetermined
        if res_locus.data["depth"] >= min_support:
            res_locus.status = Status.stable
            if res_locus.data["pro_p"] > meta["MSIsensor-pro"]["pro_p_threshold"]:
                res_locus.status = Status.unstable
        spl.addLocus(locus)
    # Process sample status
    spl_status = Status.undetermined
    nb_determined = spl.getNbDetermined(method)
    if nb_determined / spl.getNbProcessed(method) > determined_rate_threshold:
        unstable_rate = spl.getNbUnstable(method) / nb_determined
        if unstable_rate > unstable_rate_threshold:
            spl_status = Status.unstable
        else:
            spl_status = Status.stable
    spl.results[method] = MSISplRes(
        spl_status,
        unstable_rate,
        method,
        {
            "aggregation_method": "instability ratio",
            "instability_threshold": unstable_rate_threshold,
            "min_voting_loci": determined_rate_threshold
        }
    )
    # Parse distributions
    if distributions_path is not None:
        for locus, meta in DistIO(distributions_path):
            prev_res = spl.loci[locus.position].results[method]
            prev_res.data["lengths"] = locus.results[method].data["lengths"]
    return spl


def stringToBinary(str):
    """
    Return MSIsensor binary representation of the sequence.

    :param str: Sequence.
    :type str: str
    :return: MSIsensor binary representation of the sequence.
    :rtype: int
    """
    idx_by_char = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
    binary = 0
    for curr_char in str:
        binary = binary << 2  # Switch to next char (decimal: 0-3 -> binary: 0-11)
        char_idx = idx_by_char[curr_char]
        binary = binary | char_idx
    return binary
