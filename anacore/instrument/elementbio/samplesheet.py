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
        self.settings = {"all": dict(), "1": dict(), "2": dict()}  # Settings by lane
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        with open(self.filepath) as fh:
            reader = csv.Reader(fh)
            spl_obj_fields = {"SampleName", "Index1", "Index2", "Description"}
            section = None
            for rec in reader:
                if rec[0].startswith("#") or "".join(rec) == "":
                    pass
                elif rec[0].startswith("["):
                    section = rec[0][1:-1].lower()
                    if section == "samples":
                        spl_fields = [elt.strip() for elt in next(reader) if elt != ""]
                else:
                    if section == "run_values":
                        self.run_values[rec[0].strip()] = rec[1].strip()
                    elif section == "settings":
                        lane = "all" if len(rec) >= 3 and rec[2] != "" else rec[2]
                        self.settings[lane][rec[0].strip()] = rec[1].strip()
                    elif section != "samples":
                        spl_rec = {key: val.strip() for key, val in zip(spl_fields, rec)}
                        barcodes = None
                        if "Index1" in spl_rec:
                            barcodes = {"index": spl_rec["Index1"]}
                        if "Index2" in spl_rec:
                            barcodes["index2"] = spl_rec["Index2"]
                        metadata = {k: v for k, v in spl_rec.items() if k not in spl_obj_fields}
                        if "SampleName" not in spl_rec or spl_rec["SampleName"] == "":
                            raise Exception("SampleName is required for each sample.")
                        self.samples.append(
                            Sample(spl_rec["SampleName"], barcodes, spl_rec.get("Description"), metadata)
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
