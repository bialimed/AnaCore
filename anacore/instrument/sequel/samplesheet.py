# -*- coding: utf-8 -*-
"""
Classes and functions for reading PacificBioscience's sequel samplesheet.

:Code example:

    List samples and index from the samplesheet

    .. highlight:: python
    .. code-block:: python

        from anacore.instrument.sequel.samplesheet import SampleSheet

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

from anacore.instrument.samplesheet import Sample as SampleBase
import csv


class Sample(SampleBase):
    """Class to manage a sample from sample sheet."""

    @property
    def barcode_name(self):
        barcode_name = None
        if "index" in self.barcodes:
            barcode_name = self.barcodes["index"]
            if "index2" in self.barcodes:
                barcode_name += "--" + self.barcodes["index2"]
        return barcode_name

    @property
    def library_basename(self):
        return self.id

    def toDict(self):
        """
        Return dict representation of the instance.

        :return: Dict representation of the instance.
        :rtype: dict
        """
        res = super().toDict()
        res["barcode_name"] = self.barcode_name
        return res


class SampleSheet:
    """Reader for SampleSheet (see: https://www.pacb.com/wp-content/uploads/SMRT_Analysis_Barcoding_Overview_v600.pdf)."""

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
            name_field = None
            barcode_field = None
            if "Barcode" in reader.fieldnames:
                name_field = "Barcode"
                barcode_field = "Barcode"
            elif "Barcodes" in reader.fieldnames:
                name_field = "Barcodes"
                barcode_field = "Barcodes"
            elif "Barcode Name" in reader.fieldnames:
                name_field = "Barcode Name"
                barcode_field = "Barcode Name"
            if "Bio Sample Name" in reader.fieldnames:
                name_field = "Bio Sample Name"
            elif "Bio Sample" in reader.fieldnames:
                name_field = "Bio Sample"
            elif "Biosample" in reader.fieldnames:
                name_field = "Biosample"
            for rec in reader:
                id = rec[name_field]
                if barcode_field is not None:
                    if "--" in rec[barcode_field]:
                        idx1, idx2 = rec[barcode_field].split("--")
                        barcodes = {"index": idx1, "index2": idx2}
                    else:
                        barcodes = {"index": rec[barcode_field]}
                metadata = dict()
                description = None
                for key in reader.fieldnames:
                    if key not in {name_field, barcode_field, "description"}:
                        metadata[key] = rec[key]
                    elif key == "description":
                        description = rec[key]
                self.samples.append(
                    Sample(id, barcodes, description, metadata)
                )

    @staticmethod
    def filterSamples(in_ss, out_ss, samples):
        filtered_samples = [spl for spl in SampleSheet(in_ss).samples if spl.id in samples]
        SampleSheet.write(out_ss, filtered_samples)

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

    @staticmethod
    def write(filepath, samples, fields_names=None):
        if fields_names is None:
            metadata_keys = set()
            for spl in samples:
                for tag in spl.metadata.keys():
                    metadata_keys.add(tag)
            fields_names = ["Barcode Name", "Bio Sample Name"] + sorted(list(metadata_keys))  # https://www.pacb.com/wp-content/uploads/SMRT_Link_User_Guide_v10.2.pdf
        with open(filepath, "w") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields_names)
            writer.writeheader()
            for spl in samples:
                spl_dict = spl.metadata
                spl_dict["Bio Sample Name"] = ""
                if spl.id != spl.barcode_name:
                    spl_dict["Bio Sample Name"] = spl.id
                spl_dict["Barcode Name"] = ""
                if spl.barcodes:
                    spl_dict["Barcode Name"] = spl.barcode_name
                writer.writerow(spl_dict)
