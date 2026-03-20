# -*- coding: utf-8 -*-
"""
Functions to manage reads headers and filenames.

:Code example:

    Get info from read ID

    .. highlight:: python
    .. code-block:: python

        from anacore.instrument.elementbio.base import getInfFromSeqID

        read_id = "AVITI1:SIDEB_ADA_2510541191:2510541191:1:10102:0144:0037:ATGCATA+CTAGC"
        print(getInfFromSeqID(read_id))

        # Result>
        # {
        #    "sequencer_name": "AVITI1",
        #    "run_name": "SIDEB_ADA_2510541191",
        #    "flowcell_id": "2510541191",
        #    "lane_id": "1",
        #    "tile_id": "10102",
        #    "x_pos": "0144",
        #    "y_pos": "0037",
        #    "umi": "ATGCATA+CTAGC"
        # }
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'


def getInfFromSeqID(seq_id):
    """
    Return sequencer name, run name, flowcell id and position of the sequence on the sequencer flowcell.

    :param sequence: The ID of the sequence provided by the sequencer.
    :type sequence: str
    :return: The sequencer id, run id, flowcell id and position of the sequence on the sequencer flowcell.
    :rtype: dict
    """
    # https://docs.elembio.io/docs/bases2fastq/outputs/#sequence-identifiers-for-sequencing
    fields = seq_id.split(":")
    if len(fields) == 7:  # AVITI's ID: @AVITI1:SIDEB_ADA_2510541191:2510541191:1:10102:0144:0037
        fields.append(None)
    # else AVITI's ID: @AVITI1:SIDEB_ADA_2510541191:2510541191:1:10102:0144:0037:ATGCATA+CTAGC
    return {
        "sequencer_name": fields[0],
        "run_name": fields[1],
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
    # AVITI's description: 1:Y:0:ATCACG
    reads_phases, kept_status, control_bits, barcode = seq_desc.split(":")
    return {
        "reads_phases": int(reads_phases),
        "is_kept": kept_status == "N",
        "control_bits": None if control_bits == "0" else int(control_bits),
        "barcode": None if barcode == "" or barcode.isnumeric() else barcode  # No indexing: The sample number
    }
