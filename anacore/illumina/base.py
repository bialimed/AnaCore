# -*- coding: utf-8 -*-
"""
Functions to manage reads headers and filenames.

:Code example:

    Get info from read ID

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.base import getInfFromSeqID

        read_id = "EAS139:136:FC706VJ:2:2104:15343:197393:ATGCATA+CTAGC"
        print(getInfFromSeqID(read_id))

        # Result>
        # {
        #    "sequencer_id": "EAS139",
        #    "run_id": "136",
        #    "flowcell_id": "FC706VJ",
        #    "lane_id": "2",
        #    "tile_id": "2104",
        #    "x_pos": "15343",
        #    "y_pos": "197393",
        #    "umi": "ATGCATA+CTAGC"
        # }

    Get platform type from instrument serial number

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.base import getPlatformFromSerialNumber

        print(getPlatformFromSerialNumber("EAS139"))

        # Result>
        # HiSeq
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import re


def getIlluminaName(name):
    """
    Return sample name used by Illumina in filename.

    :param name: The name provided to Illumina process (for example in samplesheet).
    :type name: str
    :return: The sample name used by Illumina as part of filename.
    :rtype: str
    """
    return name.replace("_", "-").replace(" ", "-").replace(".", "-").replace("+", "")


def getLibNameFromReadsPath(seq_path):
    """
    Return library name from the path of the sequences file.

    :param seq_path: The path of the sequences file.
    :type seq_path: str
    :return: The library name.
    :rtype: str
    """
    library_name, extensions = os.path.basename(seq_path).split(".", 1)
    for curr_ext in extensions.split("."):
        if curr_ext not in ["fq", "fastq", "fasta", "fa", "gz", "bz", "bz2", "lz", "zip"]:
            raise Exception('The file "{}" cannot be processed by getLibNameFromReadsPath because the extension "{}" is not managed.'.format(seq_path, curr_ext))
    if re.search(r'_[rR][1-2]$', library_name):
        library_name = library_name[:-3]
    elif re.search(r'_[rR][1-2]_\d\d\d$', library_name):
        library_name = library_name[:-7]
    return library_name


def getInfFromSeqID(seq_id):
    """
    Return sequencer id, run id, flowcell id and position of the sequence on the sequencer flowcell.

    :param sequence: The ID of the sequence provided by the sequencer.
    :type sequence: str
    :return: The sequencer id, run id, flowcell id and position of the sequence on the sequencer flowcell.
    :rtype: dict
    """
    fields = seq_id.split(":")
    if len(fields) == 7:  # Illumina's ID: EAS139:136:FC706VJ:2:2104:15343:197393
        fields.append(None)
    # else  Illumina's ID: EAS139:136:FC706VJ:2:2104:15343:197393:ATGCATA+CTAGC
    return {
        "sequencer_id": fields[0],
        "run_id": fields[1],
        "flowcell_id": fields[2],
        "lane_id": int(fields[3]),
        "tile_id": int(fields[4]),
        "x_pos": int(fields[5]),
        "y_pos": int(fields[6]),
        "umi": fields[7]
    }


def getInfFromSeqDesc(seq_desc):
    """
    Return sequencer id, run id, flowcell id and position of the sequence on the sequencer flowcell.

    :param sequence: The ID of the sequence provided by the sequencer.
    :type sequence: str
    :return: The sequencer id, run id, flowcell id and position of the sequence on the sequencer flowcell.
    :rtype: dict
    """
    # Illumina's description: 1:Y:18:ATCACG
    reads_phases, kept_status, control_bits, barcode = seq_desc.split(":")
    return {
        "reads_phases": int(reads_phases),
        "is_kept": kept_status == "N",
        "control_bits": None if control_bits == "0" else int(control_bits),
        "barcode": None if barcode == "" else barcode
    }


def getPlatformFromSerialNumber(instrument_id):
    """
    Return platform name from instrument ID.

    :param instrument_id: The instrument serial number.
    :type instrument_id: str
    :return: The platform name (Hiseq or NextSeq or ...).
    :rtype: str
    """
    platform_by_re = {
        # "?": "iSeq",
        "^MN[0-9]{5}$": "MiniSeq",
        "^ML-..-[0-9]{2}$": "MiniSeq",
        "^M[0-9]{5}$": "MiSeq",
        "^N[SB][0-9]{6}$": "NextSeq",
        "^NDX[0-9]{6}": "NextSeq",
        "^[CDJKE][0-9]{5}$": "HiSeq",
        "^A[0-9]{5}$": "NovaSeq"
    }
    if instrument_id.startswith("HWI"):
        instrument_id = instrument_id[3:]
    platform = None
    for curr_re, curr_instru in platform_by_re.items():
        if platform is None:
            if re.search(curr_re, instrument_id):
                platform = curr_instru
    return platform
