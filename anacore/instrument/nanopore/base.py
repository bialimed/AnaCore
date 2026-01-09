# -*- coding: utf-8 -*-
"""
Functions to manage reads headers and filenames.

:Code example:

    Get info from read description

    .. highlight:: python
    .. code-block:: python

        from anacore.nanopore.base import getInfFromSeqDesc

        read_header = "@945da76d-b7ef-4c65-91a1-e69c85832185 runid=2cf9e7c4dd1888d9a1090ebc394dd5b7cb2fe3f0 read=126 ch=181 start_time=2017-10-17T06:43:19Z barcode=barcode01"
        print(getInfFromSeqDesc(read_header.split(" ", 1)[1]))
        # Result>
        # {
        #    "run_id": "2cf9e7c4dd1888d9a1090ebc394dd5b7cb2fe3f0",
        #    "read_number": 126,
        #    "channel_id": 181,
        #    "start_time": "2017-10-17T06:43:19Z",
        #    "barcode_id": "barcode01",
        #    "parent_read_id": None
        # }
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'


def getInfFromSeqDesc(seq_desc):
    """
    Return information (run id, channel, barcode, ..) from read header's description.

    :param seq_desc: The description part of read header.
    :type seq_desc: str
    :return: Information (run id, channel, barcode, ..) from read header's description.
    :rtype: dict
    """
    # NanoPore's description: runid={run_id} read={read_number} ch={channel_id} start_time={start_time_utc} [barcode={barcode_id} parent_read_id={parent_read_id}]
    # https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/input-and-output-files
    info = {
        "run_id": None,
        "read_number": None,
        "channel_id": None,
        "start_time": None,
        "barcode_id": None,
        "parent_read_id": None
    }
    for field in seq_desc.split(" "):
        tag, value = field.split("=")
        if tag == "runid":
            tag = "run_id"
        elif tag == "read":
            tag = "read_number"
            value = int(value)
        elif tag == "ch":
            tag = "channel_id"
            value = int(value)
        if tag in {"barcode", "barcodeid"}:
            tag = "barcode_id"
        if value == "":
            value = None
        info[tag] = value
    return info
