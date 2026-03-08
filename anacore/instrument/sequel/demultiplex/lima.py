# -*- coding: utf-8 -*-
"""Classes and functions to manage lima inputs and outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

from anacore.instrument.sequel.samplesheet import SampleSheet
import os


def filterSamples(in_ss, out_ss, in_demult_folder, samples, file_name_pattern="lima_bcless.{}.bam"):
    """
    Filter samplesheet to keep only samples in provided list and found in demultiplex step.

    :param in_ss: Path to original samplesheet (format: CSV).
    :type in_ss: str
    :param out_ss: Path to filtered samplesheet (format: CSV).
    :type out_ss: str
    :param in_demult_folder: Path to lima output folder.
    :type in_demult_folder: str
    :param samples: List of selectected samples.
    :type samples: list
    :param file_name_pattern: Pattern used to barcodes file name after demultiplex.
    :type file_name_pattern: str
    """
    filtered_samples = list()
    for spl in SampleSheet(in_ss).samples:
        spl_path = os.path.join(in_demult_folder, file_name_pattern.format(spl.barcode_name))
        if spl.id in samples and os.path.exists(spl_path):
            filtered_samples.append(spl)
    SampleSheet.write(out_ss, filtered_samples)


def samplesheetToLima(in_ss, out_ss):
    """
    Write lima SampleSheet (for parameter biosample-csv) from standard SampleSheet.

    :param in_ss: Path to original SampleSheet (format: CSV).
    :type in_ss: str
    :param out_ss: Path to lima biosample file (format: CSV).
    :type out_ss: str
    """
    with open(out_ss, "w") as writer:
        writer.write("Barcodes,Bio Sample\n")
        for spl in SampleSheet(in_ss).samples:
            writer.write("{},{}\n".format(spl.barcode_name, spl.basename))
