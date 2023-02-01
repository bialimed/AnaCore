# -*- coding: utf-8 -*-
"""
Classes and functions for reading Illumina's samplesheet.

:Code example:

    List samples from the samplesheet

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.samplesheet import SampleSheetFactory

        samplesheet = SampleSheetFactory.get("my_run/SampleSheet.csv"):
        print("Samples:")
        for spl in samplesheet.samples:
            print("", spl["Library_Basename"])

        # Result>
        # Samples:
        #   sapmleA_S1
        #   sampleB_S2
        #   sampleC_S3
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
        :return: Samples list containing for each: keys and values from data section as dict.
        :rtype: list
        """
        titles, samples = self._getInfoFromTitledSection(data_section)
        for idx, spl in enumerate(samples):
            spl["Sample_Basename"] = getIlluminaName(spl["Sample_ID"])
            spl["Library_Basename"] = spl["Sample_Basename"] + "_S" + str(idx + 1)
            if "Sample_Description" in spl:
                spl["Description"] = spl["Sample_Description"]
        return samples

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
        self.extra = dict()
        for title, data in sections_by_title.items():
            if title not in cls.REQUIRED_SECTIONS:
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
        # Post process
        if "Description" not in self.header:
            self.header["Description"] = ""


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
        samples = super()._getSamplesFromData(data_section)
        for spl_idx, spl in enumerate(samples):
            spl["Sample_Basename"] = getIlluminaName(spl["Sample_Name"])
            spl["Library_Basename"] = spl["Sample_Basename"] + "_S" + str(spl_idx + 1)
        return samples

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
            if spl["Manifest"] in new_by_old_id:
                if spl["Sample_ID"].endswith("_" + spl["Manifest"]):
                    spl["Sample_ID"] = spl["Sample_ID"].rsplit("_" + spl["Manifest"], 1)[0]
                    spl["Sample_ID"] += "_" + new_by_old_id[spl["Manifest"]]
                spl["Manifest"] = new_by_old_id[spl["Manifest"]]
                new_spl.append(spl)
        self.samples = new_spl

    def findSplFiles(self, directory, end_pattern, subject="library"):
        """
        Return by sample name or library name all the files in directory corresponding to the subject.

        :param directory: The path to the directory where the files names are evaluated.
        :type directory: str
        :param end_pattern: File pattern added at the end of each sample basename or library basename to select files.
        :type end_pattern: str
        :param subject: "library" or "sample" to select files corresponding to libraries or corresponding to samples.
        :type subject: str
        :return: By subject the paths of the corresponding files.
        :rtype: dict
        """
        subject_start_tag = subject.capitalize()
        files_by_elt = {}
        for spl in self.samples:
            if subject == "library" or spl["Sample_Name"] not in files_by_elt:  # The subject is the library or the subject is the sample and no library has been already processed for this sample
                files_by_elt[spl[subject_start_tag + "_Name"]] = sorted(glob.glob(
                    os.path.join(directory, spl[subject_start_tag + "_Basename"] + end_pattern)
                ))
        return files_by_elt

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
        subject_tag = subject.capitalize() + "_Name"
        files_by_spl = self.findSplFiles(directory, end_pattern, subject)
        for spl in self.samples:
            spl[tag] = files_by_spl[spl[subject_tag]]
