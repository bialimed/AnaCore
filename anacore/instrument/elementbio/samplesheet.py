# -*- coding: utf-8 -*-
"""
Classes and functions for reading AVITI's run manifest.

:Code example:

    List samples and index from the samplesheet

    .. highlight:: python
    .. code-block:: python

        from anacore.instrument.aviti.samplesheet import SampleSheet

        samplesheet = SampleSheet("my_run/RunManifest.csv"):
        print("Samples:")
        for spl in samplesheet.samples:
            print("", spl.basename, spl.barcodes)

        # Result>
        # Samples:
        #  sapmleA {"index": "TATCTTCAGC", "index2": "CGAATATTGG"}
        #  sampleB {"index": "TGCACGGATA", "index2": "GTTGTCTGCG"}
        #  sampleC {"index": "GGTTGATAGA", "index2": "TCTGCTCCGA"}
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

from anacore.instrument.samplesheet import cleanedEnd, Sample
import csv
import re


class SampleSheet:
    """Reader for SampleSheet (see: https://docs.elembio.io/docs/run-manifest/prepare-manifest)."""

    def __init__(self, path):
        """
        Build and return an instance of SampleSheet.

        :param path: Path to samplesheet (format: CSV).
        :type path: str
        :return: The new instance.
        :rtype: SampleSheet
        """
        self.filepath = path
        self.run_values = dict()
        self.samples = list()
        self.settings = dict()
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        with open(self.filepath) as fh:
            reader = csv.reader(fh)
            spl_obj_fields = {"SampleName", "Index1", "Index2", "Description"}
            section = None
            fields_by_section = {}
            for rec in reader:
                if len(rec) == 0 or rec[0].startswith("#") or "".join(rec) == "":
                    pass
                elif rec[0].startswith("["):
                    section = rec[0][1:-1].lower()
                    fields_by_section[section] = [elt.strip() for elt in next(reader) if elt != ""]
                else:
                    rec_info = {key: val.strip() for key, val in zip(fields_by_section[section], rec)}
                    if section == "runvalues":
                        self.run_values[rec_info["KeyName"]] = rec_info["Value"]
                    elif section == "settings":
                        lanes = rec_info["Lane"].split("+") if "Lane" in rec_info and rec_info["Lane"] != "" else ["1", "2"]
                        self.settings[rec_info["SettingName"]] = {curr: rec_info["Value"] for curr in lanes}
                    elif section == "samples":
                        if "SampleName" not in rec_info or rec_info["SampleName"] == "":
                            raise Exception("SampleName is required for each sample.")
                        barcodes = None
                        if "Index1" in rec_info:
                            barcodes = {"index": rec_info["Index1"]}
                        if "Index2" in rec_info:
                            barcodes["index2"] = rec_info["Index2"]
                        description = rec_info["Description"] if "Description" in rec_info and rec_info["Description"] != "" else None
                        metadata = {k: v for k, v in rec_info.items() if k not in spl_obj_fields}
                        if "Lane" in metadata:
                            metadata["Lane"] = metadata["Lane"].split("+")
                        self.samples.append(
                            Sample(rec_info["SampleName"], barcodes, description, metadata)
                        )

    @property
    def header(self):
        return self.run_values

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in SamplSheet format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in SamplSheet format.
        :rtype: bool
        """
        is_valid = False
        required_sections = {"Samples"}
        authorized_sections = {"RunValues", "Settings", "Samples"}
        try:
            sections = set()
            with open(filepath) as reader:
                curr_title = None
                for line in reader:
                    cleaned_line = cleanedEnd(line)
                    if re.fullmatch(r"\[[^]]+\]", cleaned_line):
                        curr_title = cleaned_line[1:-1]
                        sections.add(curr_title)
                    elif not cleaned_line.startswith("#") and cleaned_line != "":
                        if curr_title == "RunValues":
                            if cleaned_line.count(",") != 2:
                                raise Exception("Rows in RunValues section must contain 2 fields: key,value.")
                        elif curr_title == "Settings":
                            if cleaned_line.count(",") != 2 and cleaned_line.count(",") != 3:
                                raise Exception("Rows in Settings section must contain 2 or 3 fields: key,value[,lane].")
            if len(required_sections - sections) != 0:
                raise Exception('Missing core sections: {}.'.format(required_sections - sections))
            if len(sections - authorized_sections) != 0:
                raise Exception('Unauthorized sections: {}.'.format(sections - authorized_sections))
            SampleSheet(filepath)
            is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid
