# -*- coding: utf-8 -*-
"""
Classes and functions for reading Nanopore's samplesheet.

:Code example:

    List samples and index from the samplesheet

    .. highlight:: python
    .. code-block:: python

        from anacore.nanopore.samplesheet import SampleSheet

        samplesheet = SampleSheet("my_run/SampleSheet.csv"):
        print("Samples:")
        for spl in samplesheet.samples:
            print("", spl.basename, spl.barcodes)

        # Result>
        # Samples:
        #  sapmleA {"index": "barcode01", "index2": "barcode11"}
        #  sampleB {"index": "barcode02", "index2": "barcode12"}
        #  sampleC {"index": "barcode03", "index2": "barcode13"}
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'

from anacore.instrument.samplesheet import Sample
import csv


class SampleSheet:
    """Reader for SampleSheet (see: https://community.nanoporetech.com/docs/prepare/library_prep_protocols/experiment-companion-minknow/v/mke_1013_v1_revdc_11apr2016/sample-sheet-upload)."""

    def __init__(self, path):
        """
        Build and return an instance of SampleSheet.

        :param path: Path to samplesheet (format: CSV).
        :type path: str
        :return: The new instance.
        :rtype: SampleSheet
        """
        self.filepath = path
        self.samples = None
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        self.samples = list()
        with open(self.filepath) as fh:
            reader = csv.DictReader(fh)
            barcode = "no"
            if "barcode" in reader.fieldnames:
                barcode = "single"
            elif "internal_barcode" in reader.fieldnames:
                barcode = "dual"
            for rec in reader:
                id = rec["sample_id"] if "sample_id" in rec else "no_sample"
                barcodes = None
                if barcode == "single":
                    barcodes = {"index": rec["barcode"]}
                    if "alias" in rec and rec["alias"]:
                        barcodes["alias"] = rec["alias"]
                elif barcode == "dual":
                    barcodes = {"index": rec["internal_barcode"], "index2": rec["external_barcode"]}
                    if "alias" in rec and rec["alias"]:
                        barcodes["alias"] = rec["alias"]
                metadata = dict()
                description = None
                for key in reader.fieldnames:
                    if key not in {"alias", "barcode", "description", "external_barcode", "internal_barcode", "sample_id"}:
                        metadata[key] = rec[key]
                    elif key == "description":
                        description = rec[key]
                self.samples.append(
                    Sample(id, barcodes, description, metadata)
                )

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
        try:
            with open(filepath) as fh:
                reader = csv.DictReader(fh)
                if len({"experiment_id", "kit"} - set(reader.fieldnames)):
                    raise Exception('Missing mandatory fields: {}'.format({"experiment_id", "kit"} - set(reader.fieldnames)))
                if len({"flow_cell_id", "position_id"} & set(reader.fieldnames)) == 0:
                    raise Exception('Missing mandatory fields: "flow_cell_id" or "position_id"')
                nb_col = len(reader.fieldnames)
                nb_spl = 0
                for rec in reader:
                    nb_spl += 1
                    if len(rec) != nb_col:
                        raise Exception('Inconsistant number of columns')
                if nb_spl == 0:
                    raise Exception('Empty samplesheet')
            is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid
