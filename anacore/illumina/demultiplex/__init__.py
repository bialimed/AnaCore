# -*- coding: utf-8 -*-
"""
Classes and functions to manage Illumina's demultiplex logs and stats files.

:Code example:

    List undetermined barcodes with a count greater than smallest sample

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.demultiplex import DemultStatFactory

        demult_stat = DemultStatFactory.get("my_demultiplex_folder")
        print("count\tbarcode")
        for barcode in demult_stat.unexpectedBarcodes():
            print(
                "{}\t{}".format(barcode["ct"], barcode["seq"])
            )

        # Result (with smallest sample count equal to 10000)>
        # count barcode
        # 14879 AATGC+TTTA
        # 10457 AGCGC+TTGA
"""

from anacore.illumina.demultiplex.bcl2fastq import DemultStat as DemultStatBcl2fastq
from anacore.illumina.demultiplex.bclconvert import DemultStat as DemultStatBclConvert
import os


class DemultStatFactory:
    """Factory to identify and return version compliant handler to DemultStat."""

    @staticmethod
    def get(folder_path):
        """
        Return instance of DemultStat from the demultiplexing folder.

        :param folder_path: Path to the demultiplexing folder.
        :type folder_path: str
        :return: Instance of DemultStat from the demultiplexing folder.
        :rtype: anacore.illumina.AbstractDemultStat
        """
        if os.path.exists(os.path.join(folder_path, "Stats", "Stats.json")):
            return DemultStatBcl2fastq(
                os.path.join(folder_path, "Stats", "Stats.json")
            )
        elif os.path.exists(os.path.join(folder_path, "Reports", "Demultiplex_Stats.csv")):
            return DemultStatBclConvert(
                os.path.join(folder_path, "Reports", "Demultiplex_Stats.csv"),
                os.path.join(folder_path, "Reports", "Top_Unknown_Barcodes.csv"),
            )
        else:
            raise IOError("The folder {} is not valid for DemultStat.".format(folder_path))
