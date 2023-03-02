# -*- coding: utf-8 -*-
"""
Classes and functions for reading Illumina's samplesheet.

:Code example:

    List samples and index from the samplesheet

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.samplesheet import SampleSheetFactory

        samplesheet = SampleSheetFactory.get("my_run/SampleSheet.csv"):
        print("Samples:")
        for spl in samplesheet.samples:
            print("", spl.library_basename, spl.barcodes)

        # Result>
        # Samples:
        #  sapmleA_S1 {"index": "TATCTTCAGC", "index2": "CGAATATTGG"}
        #  sampleB_S2 {"index": "TGCACGGATA", "index2": "GTTGTCTGCG"}
        #  sampleC_S3 {"index": "GGTTGATAGA", "index2": "TCTGCTCCGA"}
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.illumina.base import getIlluminaName
import glob
import os
import re


class AbstractSampleSheet(object):
    """Abstract reader for SampleSheet.csv."""

    # REQUIRED_SECTIONS = {"Header", "Reads", "Settings", "Data"}
    # DEMULTIPLEX_SECTION = "Reads"
    # SAMPLES_SECTION = "Data"

    def __init__(self, path):
        """
        Build and return an instance of AbstractSampleSheet.

        :param path: Path to samplesheet (format: CSV).
        :type path: str
        :return: The new instance.
        :rtype: AbstractSampleSheet
        """
        self.extra = None
        self.filepath = path
        self.header = None
        self.manifests = None
        self.reads = None
        self.samples = None
        self._parse()

    @staticmethod
    def cleanedEnd(value):
        """
        Return line without ending empty CSV fields. This method can be used to removes invalid comma used in CSV to obtain the same number of columns in whole file.

        :param value: Line to clean.
        :type value: str
        :return: Line without ending empty CSV fields.
        :rtype: str
        """
        value = value.strip()
        while value.endswith(","):
            value = value[:-1].strip()
        return value

    def _getInfoFromSection(self, section):
        """
        Return keys and values contained in section.

        :param section: Lines from section. Each row contains one key and corresponding value.
        :type section: list
        :return: Keys and values contained in section.
        :rtype: dict
        """
        info = dict()
        for line in section:
            key, value = [field.strip() for field in line.split(",", 1)]
            value = AbstractSampleSheet.cleanedEnd(value)
            info[key] = value
        return info

    def _getInfoFromTitledSection(self, section):
        """
        Return titles and list of entries (Keys and values for one element).

        :param section: Lines from section with a titles row. Each row contains one entry where columns correspond to key and cells to values.
        :type section: list
        :return: Titles and list of entries (Keys and values for one element).
        :rtype: (list, list)
        """
        titles_line = AbstractSampleSheet.cleanedEnd(section[0])
        titles = titles_line.split(",")
        rows = list()
        for line in section[1:]:
            fields = [elt.strip() for elt in line.split(",")]
            while len(fields) < len(titles):
                fields.append("")
            rows.append(
                {title: fields[idx] for idx, title in enumerate(titles)}
            )
        return titles, rows

    def _getReads(self, reads_section):
        """
        Returns reads phase(s). Each phase is described by its type (is_index) and by its number of cycles (nb_cycles).

        :param reads_section: Lines from desction "Reads".
        :type reads_section: list
        :return: Reads phase(s). Each phase is described by its type (is_index) and by its number of cycles (nb_cycles).
        :rtype: list
        """
        raise NotImplementedError

    def _getSamplesFromData(self, data_section):
        """
        Return samples list containing for each: keys and values from data section as dict.

        :param section: Lines from section with a titles row. Each row contains one entry where columns correspond to key and cells to values.
        :type section: list
        :return: Samples list as anacore.illumina.samplesheet.Sample.
        :rtype: list
        """
        samples_obj = list()
        titles, samples = self._getInfoFromTitledSection(data_section)
        protected_fields = {"sample_id", "sample_description", "description", "index", "index2"}
        spl_idx = 1  # Index come from position in Data (for V1) or BCLConvert_Data (for V3) list afterremove same Sample_ID
        spl_idx_by_id = dict()
        lib_by_barcode = dict()
        for spl in samples:
            spl_id = spl["Sample_ID"]
            barcodes = dict()
            if "index" in spl:
                barcodes["index"] = spl["index"]
            if "index2" in spl:
                if "index" not in spl:
                    barcodes["index"] = ""
                barcodes["index2"] = spl["index2"]
            barcodes_str = "+".join([str(val) for key, val in sorted(barcodes.items())])
            if barcodes_str in lib_by_barcode:  # Same lib on multiple lanes
                if "Lane" in spl:
                    lib_by_barcode[barcodes_str].metadata["Lane"] += "," + spl["Lane"]
                else:
                    raise Exception("Sample {} is repeated with same barcodes: {}".format(spl_id, barcodes_str))
            else:  # New library
                # Get description
                description = None
                if "Sample_Description" in spl:
                    description = spl["Sample_Description"]
                elif "Description" in spl:
                    description = spl["Description"]
                # Get metadata
                meta = dict()
                for key, val in spl.items():
                    lower_key = key.lower()
                    if lower_key not in protected_fields:
                        meta[key] = val
                # Store library
                lib_spl_idx = spl_idx
                if spl_id in spl_idx_by_id:  # Different libraries with same sample ID
                    lib_spl_idx = spl_idx_by_id[spl_id]  # Libraries with same sample ID keep the same S{idx}
                else:  # New sample
                    spl_idx_by_id[spl_id] = lib_spl_idx
                    spl_idx += 1
                lib_obj = Sample(lib_spl_idx, spl_id, barcodes, description, meta)
                samples_obj.append(lib_obj)
                lib_by_barcode[barcodes_str] = lib_obj
        return samples_obj

    def _getSections(self):
        """
        Read file content and return by section its lines.

        :return: By section title its content as list of rows.
        :rtype: dict
        """
        sections_by_title = dict()
        with open(self.filepath) as reader:
            sections_by_title = dict()
            section_title = None
            for line in reader:
                striped_line = line.strip()
                if re.fullmatch(r"\[[^]]+\],*", striped_line):  # New section
                    striped_line = AbstractSampleSheet.cleanedEnd(striped_line)
                    section_title = striped_line[1:-1]
                    sections_by_title[section_title] = list()
                elif re.fullmatch(r",*", striped_line):  # Empty line
                    section_title = None
                else:  # Section content
                    sections_by_title[section_title].append(striped_line)
        return sections_by_title

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        # Retrieve information by section
        sections_by_title = self._getSections()
        # Process information
        cls = self.__class__
        self.header = self._getInfoFromSection(
            sections_by_title["Header"]
        )
        self.manifests = dict()
        if "Manifests" in sections_by_title:
            self.manifests = self._getInfoFromSection(sections_by_title["Manifests"])
        self.samples = self._getSamplesFromData(
            sections_by_title[cls.SAMPLES_SECTION]
        )
        self.reads = self._getReads(
            sections_by_title[cls.DEMULTIPLEX_SECTION]
        )
        self.extra = dict()
        for title, data in sections_by_title.items():
            if title not in cls.REQUIRED_SECTIONS and title != "Manifests":
                if title.endswith("_Data") or title.endswith("_Settings"):
                    title, sub_title = title.rsplit("_", 1)
                    sub_title = sub_title.lower()
                    if title not in self.extra:
                        self.extra[title] = dict()
                    if sub_title == "data":
                        self.extra[title][sub_title] = self._getInfoFromTitledSection(data)[1]
                    else:
                        self.extra[title][sub_title] = self._getInfoFromSection(data)
                else:
                    self.extra[title] = self._getInfoFromSection(data)
        # Post process
        if "Description" not in self.header:
            self.header["Description"] = ""


class Sample:
    """Class to manage a sample from sample sheet. One sample could represent several libraries if they are the same sample ID."""

    def __init__(self, sheet_index, id, barcodes=None, description=None, metadata=None):
        """
        Build and return an instance of Sample.

        :param sheet_index: Index of sample in samplesheet (start from 1). This element is used for library name: <splID>_S<splSheetIndex>.
        :type sheet_index: int
        :param id: Sample ID (example: splA).
        :type id: str
        :param barcodes: Sample molecular index used in demultiplex. Example:
            * For single-end: {"index": "ATGC"}
            * For paired-end: {"index": "ATGC", "index2": "CGCA"}
            * For paired-end with skipped index 1: {"index": "", "index2": "CGCA"}
        :type barcodes: dict
        :param description: Sample description.
        :type description: str
        :param metadata: Additional data on sample (example: {"panel": "Enrich_Soli-Tumor_V3", "project": "Rimo"}). These data come from additional column in section Data.
        :type metadata: dict
        :return: The new instance.
        :rtype: Sample
        """
        self.barcodes = dict() if barcodes is None else barcodes
        self.description = description
        self.id = id
        self.sheet_index = sheet_index
        self.metadata = dict() if metadata is None else metadata

    @property
    def basename(self):
        """
        Return the start of file basename corresponding to the sample (example: splA).

        :return: Sample basename.
        :rtype: str
        """
        return getIlluminaName(self.id)

    @property
    def library_basename(self):
        """
        Return the start of file basename corresponding to the library (example: splA_S3).

        :return: Library basename.
        :rtype: str
        """
        return "{}_S{}".format(self.basename, self.sheet_index)

    def toDict(self):
        """
        Return dict representation of the instance.

        :return: Dict representation of the instance.
        :rtype: dict
        """
        res = self.__dict__
        res["basename"] = self.basename
        res["library_basename"] = self.library_basename
        return res


class SampleADS(Sample):
    """Class to manage an amplicon double strand sample from sample sheet. One sample could represent several libraries if they are the same sample ID."""

    @property
    def basename(self):
        """
        Return the start of file basename corresponding to the sample (example: splA).

        :return: Sample basename.
        :rtype: str
        """
        return getIlluminaName(self.metadata["Sample_Name"])

    def findSplFiles(self, directory, end_pattern, subject="library"):
        """
        Return paths to files in directory corresponding to the subject.

        :param directory: The path to the directory where the files names are evaluated.
        :type directory: str
        :param end_pattern: File pattern added at the end of sample basename or library basename to select files.
        :type end_pattern: str
        :param subject: "library" or "sample" to select files corresponding to libraries or corresponding to samples.
        :type subject: str
        :return: Paths of the corresponding files.
        :rtype: list
        """
        selector = self.libray_basename if subject == "library" else self.basename
        return sorted(
            glob.glob(
                os.path.join(directory, selector + end_pattern)
            )
        )

    def setSplFiles(self, tag, directory, end_pattern, subject="library"):
        """
        Add the files in directory corresponding to the subject. These information are keyed by provided tag for each element of self.samples.

        :param tag: The attribute used in self.sample to store the files paths.
        :type tag: str
        :param directory: The path to the directory where the files names are evaluated.
        :type directory: str
        :param end_pattern: File pattern added at the end of each sample basename or library basename to select files.
        :type end_pattern: str
        :param subject: "library" or "sample" to select files corresponding to libraries or corresponding to samples.
        :type subject: str
        """
        self.metadata[tag] = self.findSplFiles(directory, end_pattern, subject)


class SampleSheetFactory(object):
    """Factory to identify and return version compliant handler to SampleSheet."""

    @staticmethod
    def get(filepath, *args, **kwargs):
        """
        Return instance of SampleSheet corresponding to the file.

        :param filepath: Path to the samplesheet (format: CSV).
        :type filepath: str
        :return: Instance of SampleSheet corresponding to the file.
        :rtype: AbstractSampleSheet
        """
        if SampleSheetV2.isValid(filepath):
            return SampleSheetV2(filepath, *args, **kwargs)
        elif SampleSheetADS.isValid(filepath):
            return SampleSheetADS(filepath, *args, **kwargs)
        elif SampleSheetV1.isValid(filepath):
            return SampleSheetV1(filepath, *args, **kwargs)
        else:
            raise IOError("The file {} does not have a valid format for SampleSheetReader.".format(filepath))


class SampleSheetV1(AbstractSampleSheet):
    """Reader for SampleSheet.csv in V1 format (see: https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/sequencing-sheet-format-specifications-technical-note-970-2017-004.pdf)."""

    REQUIRED_SECTIONS = {"Header", "Reads", "Settings", "Data"}
    DEMULTIPLEX_SECTION = "Reads"
    SAMPLES_SECTION = "Data"

    def _getReads(self, reads_section):
        """
        Returns number of cycles in reads phase(s). Index cycles are ignored.

        :param reads_section: Lines from desction "Reads".
        :type reads_section: list
        :return: Number of cycles in reads phase(s). Index cycles are ignored.
        :rtype: dict
        """
        cycles = dict()
        read_idx = 1
        for line in reads_section:
            striped_line = AbstractSampleSheet.cleanedEnd(line)
            cycles["R{}".format(read_idx)] = int(striped_line)
            read_idx += 1
        return {"nb_cycles": cycles}

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in SamplSheet V1 format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in SamplSheet V1 format.
        :rtype: bool
        """
        is_valid = False
        try:
            sections = set()
            with open(filepath) as reader:
                for line in reader:
                    striped_line = line.strip()
                    if re.fullmatch(r"\[[^]]+\],*", striped_line):
                        title = AbstractSampleSheet.cleanedEnd(striped_line)[1:-1]
                        sections.add(title)
                    elif striped_line.startswith("FileFormatVersion,"):
                        if striped_line.split(",")[1] != "1":
                            raise Exception('Invalid format declaration.')
            if len(SampleSheetV1.REQUIRED_SECTIONS - sections) != 0:
                raise Exception('Missing core sections: {}.'.format(SampleSheetV1.REQUIRED_SECTIONS - sections))
            is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid


class SampleSheetV2(AbstractSampleSheet):
    """Reader for SampleSheet.csv in V2 format (see: https://sapac.support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/NextSeq2000/SampleSheet_v2_template.csv)."""

    REQUIRED_SECTIONS = {"Header", "Reads", "BCLConvert_Data", "BCLConvert_Settings"}
    DEMULTIPLEX_SECTION = "Reads"
    SAMPLES_SECTION = "BCLConvert_Data"

    def _getReads(self, reads_section):
        """
        Returns reads phase(s). Each phase is described by its type (is_index) and by its number of cycles (nb_cycles).

        :param reads_section: Lines from desction "Reads".
        :type reads_section: list
        :return: Reads phase(s). Each phase is described by its type (is_index) and by its number of cycles (nb_cycles).
        :rtype: list
        """
        reads_info = self._getInfoFromSection(reads_section)
        reads_phases = list()
        for phase in ["Read1Cycles", "Index1Cycles", "Index2Cycles", "Read2Cycles"]:
            if phase in reads_info:
                reads_phases.append(
                    {'is_index': phase.startswith("Index"), 'nb_cycles': int(reads_info[phase])}
                )
        return {"phases": reads_phases}

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        super()._parse()
        # Add samples data from software data to self.samples like in V1 format
        spl_by_id = dict()
        for spl in self.samples:
            if spl.id not in spl_by_id:
                spl_by_id[spl.id] = list()
            spl_by_id[spl.id].append(spl)
        for title, curr_extra in self.extra.items():
            if "data" in curr_extra and isinstance(curr_extra["data"], list):
                if len(curr_extra["data"]) != 0:
                    lower_keys = [key.lower() for key in curr_extra["data"][0].keys()]
                    if "sample_id" in lower_keys:  # Section contains data by sample
                        for spl_metadata in curr_extra["data"]:
                            spl_id = [val for key, val in spl_metadata.items() if key.lower() == "sample_id"][0]
                            for key, val in spl_metadata.items():
                                lower_key = key.lower()
                                if lower_key != "sample_id":
                                    if lower_key in {"sample_description", "description"}:
                                        for curr_spl in spl_by_id[spl_id]:
                                            curr_spl.description = val
                                    else:
                                        for curr_spl in spl_by_id[spl_id]:
                                            curr_spl.metadata[key] = val

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in SamplSheet V2 format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in SamplSheet V2 format.
        :rtype: bool
        """
        is_valid = False
        try:
            format_declaration = None
            sections = set()
            with open(filepath) as reader:
                for line in reader:
                    striped_line = line.strip()
                    if re.fullmatch(r"\[[^]]+\],*", striped_line):
                        title = AbstractSampleSheet.cleanedEnd(striped_line)[1:-1]
                        sections.add(title)
                    elif striped_line.startswith("FileFormatVersion,2"):
                        format_declaration = striped_line
            if format_declaration is None:
                raise Exception('Invalid format declaration.')
            if len(SampleSheetV2.REQUIRED_SECTIONS - sections) != 0:
                raise Exception('Missing core sections: {}.'.format(SampleSheetV2.REQUIRED_SECTIONS - sections))
            for extra_section in sections - SampleSheetV2.REQUIRED_SECTIONS:
                if not extra_section.endswith("_Data") and not extra_section.endswith("_Settings"):
                    raise Exception('Section "{}" is invalid.'.fomat(extra_section))
            is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid


class SampleSheetADS(SampleSheetV1):
    """Reader for SampleSheet.csv designed for AmpliconDS analysis."""

    def _getSamplesFromData(self, data_section):
        """
        Return samples list containing for each: keys and values from data section as dict.

        :param section: Lines from section with a titles row. Each row contains one entry where columns correspond to key and cells to values.
        :type section: list
        :return: Samples list as anacore.illumina.samplesheet.SampleADS.
        :rtype: list
        """
        samples = super()._getSamplesFromData(data_section)
        return [
            SampleADS(spl.sheet_index, spl.id, spl.barcodes, spl.description, spl.metadata)
            for spl in samples
        ]

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in SamplSheet Amplicon Double Strand format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in SamplSheet Amplicon Double Strand format.
        :rtype: bool
        """
        is_valid = False
        if SampleSheetV1.isValid(filepath):
            try:
                samplesheet = SampleSheetV1(filepath)
                if "Application" in samplesheet.header and \
                   samplesheet.header["Application"] == "Amplicon - DS" and \
                   len(samplesheet.manifests) != 0:
                    if len(samplesheet.samples) == 0:
                        is_valid = True
                    elif "Manifest" in samplesheet.samples[0].metadata:
                        is_valid = True
            except FileNotFoundError:
                raise
            except Exception:
                pass
        return is_valid

    def filterPanels(self, selected_manifests):
        """
        Filter samples and manifests on their panel.

        :param selected_manifests: The manifests corresponding to the selected panel.
        :type selected_manifests: list
        """
        # Get kept
        new_manifests = {}
        new_by_old_id = {}
        new_idx = 0
        for manifest_id in sorted(self.manifests):
            manifest_file = self.manifests[manifest_id]
            if manifest_file in selected_manifests:
                new_manifest_id = chr(ord("A") + new_idx)
                new_by_old_id[manifest_id] = new_manifest_id
                new_manifests[new_manifest_id] = self.manifests[manifest_id]
                new_idx += 1
        # Update manifests
        self.manifests = new_manifests
        # Update samples
        new_spl = list()
        for spl in self.samples:
            if spl.metadata["Manifest"] in new_by_old_id:
                if spl.id.endswith("_" + spl.metadata["Manifest"]):
                    spl.id = spl.id.rsplit("_" + spl.metadata["Manifest"], 1)[0]
                    spl.id += "_" + new_by_old_id[spl.metadata["Manifest"]]
                spl.metadata["Manifest"] = new_by_old_id[spl.metadata["Manifest"]]
                new_spl.append(spl)
        self.samples = new_spl
