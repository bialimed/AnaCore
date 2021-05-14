# -*- coding: utf-8 -*-
"""Utils to works with human genome accessions."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


class ChrAccession:
    """
    Manage Human chromosome accessions and names.

    :Code example: Compare two positions from different sources

        .. highlight:: python
        .. code-block:: python

            from anacore.db.homo_sapiens.accession import ChrAccession

            first_GRCh38_pos = "NC_000001.11:154545-154545"  # from RefSeq
            first_chr, first_pos = first_GRCh38_pos.split(":")
            first_chr_std = ChrAccession.toHumanName(first_chr)  # "1"

            second_GRCh38_pos = "1:154545-154545"  # from Ensembl
            second_chr, second_pos = first_GRCh38_pos.split(":")
            second_chr_std = ChrAccession.toHumanName(second_chr)  # "1"

            if first_chr_std == second_chr_std and first_pos == second_pos:
                print("Same position")
            else:
                print("Different position")

            # Result>
            # Same position
    """

    HUMAN_BY_GENBANK = {
        "CM000663": "1",
        "CM000664": "2",
        "CM000665": "3",
        "CM000666": "4",
        "CM000667": "5",
        "CM000668": "6",
        "CM000669": "7",
        "CM000670": "8",
        "CM000671": "9",
        "CM000672": "10",
        "CM000673": "11",
        "CM000674": "12",
        "CM000675": "13",
        "CM000676": "14",
        "CM000677": "15",
        "CM000678": "16",
        "CM000679": "17",
        "CM000680": "18",
        "CM000681": "19",
        "CM000682": "20",
        "CM000683": "21",
        "CM000684": "22",
        "CM000685": "X",
        "CM000686": "Y",
        "J01415": "MT"  # GRCh38
    }

    HUMAN_BY_REFSEQ = {
        "NC_000001": "1",
        "NC_000002": "2",
        "NC_000003": "3",
        "NC_000004": "4",
        "NC_000005": "5",
        "NC_000006": "6",
        "NC_000007": "7",
        "NC_000008": "8",
        "NC_000009": "9",
        "NC_000010": "10",
        "NC_000011": "11",
        "NC_000012": "12",
        "NC_000013": "13",
        "NC_000014": "14",
        "NC_000015": "15",
        "NC_000016": "16",
        "NC_000017": "17",
        "NC_000018": "18",
        "NC_000019": "19",
        "NC_000020": "20",
        "NC_000021": "21",
        "NC_000022": "22",
        "NC_000023": "X",
        "NC_000024": "Y",
        "NC_0018007": "MT",  # GRCh37
        "NC_012920": "MT"  # GRCh38
    }

    HUMAN_BY_HUMAN_CHR = {
        "chr1": "1",
        "chr2": "2",
        "chr3": "3",
        "chr4": "4",
        "chr5": "5",
        "chr6": "6",
        "chr7": "7",
        "chr8": "8",
        "chr9": "9",
        "chr10": "10",
        "chr11": "11",
        "chr12": "12",
        "chr13": "13",
        "chr14": "14",
        "chr15": "15",
        "chr16": "16",
        "chr17": "17",
        "chr18": "18",
        "chr19": "19",
        "chr20": "20",
        "chr21": "21",
        "chr22": "22",
        "chrX": "X",
        "chrY": "Y",
        "chrM": "MT",  # GRCh37
        "chrMT": "MT"  # GRCh38
    }

    HUMAN_BY_HUMAN_WO_CHR = {
        "1": "1",
        "2": "2",
        "3": "3",
        "4": "4",
        "5": "5",
        "6": "6",
        "7": "7",
        "8": "8",
        "9": "9",
        "10": "10",
        "11": "11",
        "12": "12",
        "13": "13",
        "14": "14",
        "15": "15",
        "16": "16",
        "17": "17",
        "18": "18",
        "19": "19",
        "20": "20",
        "21": "21",
        "22": "22",
        "X": "X",
        "Y": "Y",
        "M": "MT",  # GRCh37
        "MT": "MT"  # GRCh38
    }

    HUMAN_BY_ACC = {
        **HUMAN_BY_GENBANK,
        **HUMAN_BY_REFSEQ,
        **HUMAN_BY_HUMAN_CHR,
        **HUMAN_BY_HUMAN_WO_CHR
    }

    @staticmethod
    def toHumanName(chr_name):
        """
        Return standardized chromosome name from chromosome accession or name.

        :param chr_name: Chromosome name or accession.
        :type chr_name: str
        :return: Standardized chromosome name from chromosome accession.
        :rtype: str
        """
        if chr_name.lower().startswith("chr"):
            chr_name = chr_name[3:]
        chr_name = chr_name.upper()
        chr_name = chr_name.split(".")[0]
        return ChrAccession.HUMAN_BY_ACC[chr_name]
