# -*- coding: utf-8 -*-
"""Classes and functions for reading Illumina's files (samplesheet, runInfo, runParameters, ...) and parsing header, filenames and run folder informations."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.21.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import re
import glob
import datetime
import xml.etree.ElementTree as ET


class SampleSheetIO(object):
    """Reader for SampleSheet.csv."""

    SECTIONS = ["Header", "Manifests", "Reads", "Settings", "Data"]

    def __init__(self, path):
        self.filepath = path
        self.run = None
        self.samples = None
        self.header = None
        self.manifests = None
        self._parse()

    def _parse(self):
        # Retrieve information by section
        with open(self.filepath) as FH_sheet:
            sections_by_title = dict()
            section_title = None
            for line in FH_sheet:
                striped_line = line.strip()
                if any([re.fullmatch("\[" + elt + "\],*", striped_line) for elt in SampleSheetIO.SECTIONS]):
                    while striped_line.endswith(","):  # Removes invalid comma used in CSV to obtain the same number of columns in whole file
                        striped_line = striped_line[:-1]
                    section_title = striped_line[1:-1]
                    sections_by_title[section_title] = list()
                elif re.fullmatch(",*", striped_line):
                    section_title = None
                else:
                    sections_by_title[section_title].append(striped_line)
        # Process information
        self.samples = self._getSamplesFromData(sections_by_title["Data"])
        self.header = self._getInfoFromSection(sections_by_title["Header"])
        self.manifests = self._getInfoFromSection(sections_by_title["Manifests"]) if "Manifests" in sections_by_title else dict()
        self.run = self._getRunFromHeaderAndReads(sections_by_title["Header"], sections_by_title["Reads"])
        # Post process
        if "Description" not in self.header:
            self.header["Description"] = ""

    def _getSamplesFromData(self, data_section):
        samples = list()
        data_titles = [field.strip() for field in data_section[0].split(",")]
        for spl_idx, line in enumerate(data_section[1:]):
            spl = {data_titles[idx]: field.strip() for idx, field in enumerate(line.split(","))}
            spl["Sample_Basename"] = getIlluminaName(spl["Sample_ID"])
            spl["Library_Basename"] = spl["Sample_Basename"] + "_S" + str(spl_idx + 1)
            if "Sample_Description" in spl:
                spl["Description"] = spl["Sample_Description"]
            samples.append(spl)
        return samples

    def _getInfoFromSection(self, section):
        info = dict()
        for line in section:
            key, value = [field.strip() for field in line.split(",", 1)]
            while value.endswith(","):  # Removes invalid comma used in CSV to obtain the same number of columns in whole file
                value = value[:-1]
            info[key] = value
        return info

    def _getRunFromHeaderAndReads(self, header_section, reads_section):
        run = dict()
        # Header
        for line in header_section:
            key, value = line.split(",", 1)
            while value.endswith(","):  # Removes invalid comma used in CSV to obtain the same number of columns in whole file
                value = value[:-1]
            run[key] = value
        # Reads
        read_idx = 1
        for line in reads_section:
            striped_line = line.strip()
            while striped_line.endswith(","):  # Removes invalid comma used in CSV to obtain the same number of columns in whole file
                striped_line = striped_line[:-1]
            run["nb_cycles_R" + str(read_idx)] = int(striped_line)
            read_idx += 1
        return run


class ADSSampleSheetIO(SampleSheetIO):
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


class Bcl2fastqLog(object):
    """Reader for bcl2fastq log."""

    def __init__(self, path):
        self.command = None
        self.complete = None
        self.end_time = None
        self.filepath = path
        self.start_time = None
        self.version = None
        self._parse()

    def _parse(self):
        self.command = None
        self.complete = None
        self.end_time = None
        self.start_time = None
        self.version = None
        last_info = ""
        with open(self.filepath) as reader:
            for line in reader:
                line = line.strip()
                if self.version is None and line.startswith("bcl2fastq v"):
                    self.version = line.split(" ")[1]
                elif self.command is None and "Command-line invocation:" in line:
                    trace_time, command = line.split(" Command-line invocation: ")
                    self.command = command
                    self.start_time = datetime.datetime.strptime(
                        trace_time.rsplit(" ", 1)[0],
                        '%Y-%m-%d %H:%M:%S'
                    )
                elif line != "":
                    last_info = line
        if "Processing completed" in last_info:
            self.end_time = datetime.datetime.strptime(
                " ".join(last_info.split(" ", 3)[:2]),
                '%Y-%m-%d %H:%M:%S'
            )
            self.complete = {
                "nb_errors": int(
                    re.search(r"(\d+) errors", last_info).group(1)
                ),
                "nb_warnings": int(
                    re.search(r"(\d+) warnings", last_info).group(1)
                )
            }


class RTAComplete(object):
    """Reader for RTAComplete.txt."""

    def __init__(self, path):
        self.filepath = path
        self.end_date = None
        self.RTA_version = None
        self._parse()

    def _parse(self):
        with open(self.filepath) as FH_in:
            content = FH_in.read().strip()
            match = re.fullmatch("(.+/.+/.+,.+),Illumina RTA (.+)", content)
            if match:  # 11/1/2017,15:11:43.174,Illumina RTA 1.18.54
                date_str, self.RTA_version = match.groups()
                if "." in date_str and ".":
                    date_str = date_str.split(".")[0]
                self.end_date = datetime.datetime.strptime(date_str, '%m/%d/%Y,%H:%M:%S')
            else:
                match = re.fullmatch("RTA (.+) completed on (.+/.+/.+ .+:.+:.+ ..)", content)
                if match:  # RTA 2.4.11 completed on 11/14/2019 4:56:45 AM
                    self.RTA_version, date_str = match.groups()
                    self.end_date = datetime.datetime.strptime(date_str, '%m/%d/%Y %I:%M:%S %p')
                else:
                    raise Exception('"{}" in {} cannot be parsed by {}.'.format(content, self.filepath, self.__class__.__name__))


class RunInfo(object):
    """Reader for RunInfo.xml."""

    def __init__(self, path):
        self.filepath = path
        self.flowcell = None  # {'id': 'HN33VBGX5', 'layout': {'LaneCount': '4', 'SurfaceCount': '2', 'SwathCount': '3', 'TileCount': '12', 'SectionPerLane': '3', 'LanePerSection': '2'}}
        self.instrument = None  # {'id': 'NS500523', 'platform': 'NextSeq'}
        self.reads_phases = None  # [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}]
        self.run = None  # {'number': '133', 'id': '180406_NS500523_0133_AHN33VBGX5', 'start_date': datetime.datetime(2018, 4, 6, 0, 0)}
        self._parse()

    def _parse(self):
        # Retrieve information by section
        tree = ET.parse(self.filepath)
        run = tree.getroot().find("Run")
        # Process information
        self.instrument = self._getInstrumentFromRun(run)
        self.reads_phases = self._getReadsFromRun(run.find("Reads"))
        self.run = self._getRunFromRun(run)
        self.flowcell = self._getFlowcellFromRun(run)

    def _getReadsFromRun(self, reads_tree):
        reads = list()
        for child in reads_tree:
            reads.append({
                "is_index": (child.get("IsIndexedRead") == "Y"),
                "nb_cycles": int(child.get("NumCycles"))
            })
        return reads

    def _getInstrumentFromRun(self, run_tree):
        serial_number = run_tree.find("Instrument").text
        return {
            "id": serial_number,
            "platform": platformFromInstrumentSerialNumber(serial_number)
        }

    def _getRunFromRun(self, run_tree):
        date_str = run_tree.find("Date").text
        start_date = None
        if len(date_str) == 6:
            start_date = datetime.datetime.strptime(date_str, '%y%m%d')
        elif len(date_str) == 8:
            start_date = datetime.datetime.strptime(date_str, '%Y%m%d')
        else:
            start_date = datetime.datetime.strptime(date_str, '%m/%d/%Y %I:%M:%S %p')
        return {
            "number": run_tree.attrib["Number"],
            "id": run_tree.attrib["Id"],
            "start_date": start_date
        }

    def _getFlowcellFromRun(self, run_tree):
        return {
            "id": run_tree.find("Flowcell").text,
            "layout": run_tree.find("FlowcellLayout").attrib
        }


class RunParameters(object):
    """Reader for runParameters.xml."""

    def __init__(self, path):
        self.filepath = path
        self.kit = None  # {"flowcell_id": "H14U5BGXX", "reagent_kit_id": "NS2044178-REAGT"}
        self.instrument = None  # {"id": "NS500413", "platform": "NextSeq"}
        self.reads_phases = None  # [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}]
        self.run = None  # {"number": "0033", "id": "141107_NS500413_0033_H14U5BGXX", "start_date": datetime}
        self.post_process = None  # "GenerateFASTQ"
        self.software = None  # {"RTA": "2.11.3", "CS": "4.0.1.41"}
        self._parse()

    def _parse(self):
        # Retrieve information by section
        tree = ET.parse(self.filepath)
        root = tree.getroot()
        if {"Version", "Setup"} == {child.tag for child in root}:  # For HiSeq information are all in <Setup> under <RunParameters>
            root = [child for child in root if child.tag == "Setup"][0]
        # Process information
        self.instrument = self._getInstrumentFromRoot(root)
        reads_subtree = root.find("Reads")
        if reads_subtree is not None:
            self.reads_phases = self._getReadsFromSection(reads_subtree)
        else:
            self.reads_phases = self._getReadsFromSetup(root.find("Setup"))
        self.run = self._getRunFromRoot(root)
        self.kit = self._getKitFromRoot(root)
        self.post_process = self._getPostProcessFromRoot(root)
        self.software = self._getSoftwareFromRoot(root)

    def _getReadsFromSetup(self, subtree):
        reads = list()
        indices = list()
        for child in subtree:
            if child.tag.startswith("Read"):
                reads.append({
                    "is_index": False,
                    "nb_cycles": int(child.text)
                })
            elif child.tag.startswith("Index"):
                indices.append({
                    "is_index": True,
                    "nb_cycles": int(child.text)
                })
        phases = [reads[0]]
        for curr in indices:
            phases.append(curr)
        for curr in reads[1:]:
            phases.append(curr)
        return phases

    def _getReadsFromSection(self, subtree):
        reads = list()
        for child in subtree:
            reads.append({
                "is_index": (child.get("IsIndexedRead") == "Y"),
                "nb_cycles": int(child.get("NumCycles"))
            })
        return reads

    def _getInstrumentFromRoot(self, root):
        serial_number = root.find("ScannerID")
        if serial_number is None:
            serial_number = root.find("InstrumentID")
        serial_number = serial_number.text
        return {
            "id": serial_number,
            "platform": platformFromInstrumentSerialNumber(serial_number)
        }

    def _getRunFromRoot(self, root):
        run_number = root.find("RunNumber")
        if run_number is None:
            run_number = root.find("ScanNumber")
        run_number = run_number.text.lstrip("0")
        return {
            "number": run_number,
            "id": root.find("RunID").text,
            "start_date": datetime.datetime.strptime(root.find("RunStartDate").text, '%y%m%d')
        }

    def _getPostProcessFromRoot(self, root):
        workflow = root.find("AnalysisWorkflowType")  # NextSeq style
        if workflow is None:
            workflow_markup = root.find("Workflow")  # MiSeq style
            if workflow_markup is not None:
                workflow = workflow_markup.find("Analysis")
            # Else HiSeq does not process any post-processing
        if workflow is not None:
            workflow = None if workflow.text == "" else workflow.text
        return workflow

    def _getKitFromRoot(self, root):
        # Flowcell ID
        flowcell_id = root.find("FlowCellSerial")  # NextSeq style
        if flowcell_id is not None:
            flowcell_id = flowcell_id.text
        else:
            flowcell_markup = root.find("FlowcellRFIDTag")  # MiSeq style
            if flowcell_markup is not None:
                flowcell_id = flowcell_markup.find("SerialNumber").text
            else:
                flowcell_id = root.find("Barcode")  # HiSeq style
                if flowcell_id is not None:
                    flowcell_id = flowcell_id.text
                else:
                    flowcell_id = self.run["id"].rsplit("_", 1)[1]
        # Reagent kit ID
        reagent_kit_id = root.find("ReagentKitBarcode")
        if reagent_kit_id is None:
            reagent_kit_id = root.find("ReagentKitSerial")
        if reagent_kit_id is not None:
            reagent_kit_id = reagent_kit_id.text
        return {
            "flowcell_id": flowcell_id,
            "reagent_kit_id": reagent_kit_id
        }

    def _getSoftwareFromRoot(self, root):
        cs_root_markup = root.find("Setup")
        if cs_root_markup is None:
            cs_root_markup = root
        return {
            "RTA": root.find("RTAVersion").text,
            "CS": cs_root_markup.find("ApplicationVersion").text
        }


class CompletedJobInfo(object):
    """Reader for CompletedJobInfo.xml."""

    def __init__(self, path, date_format='%Y-%m-%dT%H:%M:%S'):
        self.filepath = path
        self.date_format = date_format
        self.start_datetime = None
        self.end_datetime = None
        self.version = None
        self.workflow_name = None
        self.parameters = None
        self._parse()

    def _parse(self):
        # Retrieve information by section
        tree = ET.parse(self.filepath)
        root = tree.getroot()
        # Process information
        self.end_datetime = root.find("CompletionTime").text  # 2017-11-23T17:19:55.0941558+01:00
        self.end_datetime = self.end_datetime.split(".")[0]
        self.end_datetime = datetime.datetime.strptime(self.end_datetime, self.date_format)
        self.start_datetime = root.find("StartTime").text  # 2017-11-23T17:19:55.0941558+01:00
        self.start_datetime = self.start_datetime.split(".")[0]
        self.start_datetime = datetime.datetime.strptime(self.start_datetime, self.date_format)
        self.version = root.find("Workflow").find("WorkflowVersion").text
        self.workflow_name = root.find("Workflow").find("Analysis").text
        workflow_param = etreeToDict(root.find("Workflow").find("WorkflowSettings"))
        self.parameters = {"WorkflowSettings": workflow_param}
        if self.workflow_name == "Amplicon - DS":
            self.parameters["AmpliconSettings"] = etreeToDict(root.find("Workflow").find("AmpliconSettings"))


def etreeToDict(node):
    data = {"key": node.tag}
    # Attributes
    if len(node.attrib):
        data["attributes"] = node.attrib
    # Children
    children = list()
    for child in node:
        children.append(etreeToDict(child))
    if len(children):
        data["children"] = children
    else:
        data["val"] = node.text
    # Return
    return data


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
    if re.search('_[rR][1-2]$', library_name):
        library_name = library_name[:-3]
    elif re.search('_[rR][1-2]_\d\d\d$', library_name):
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


def platformFromInstrumentSerialNumber(instrument_id):
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


def getRunFolderInfo(run_folder):
    """
    Return run information (instrument, run configuration, samples, ...) coming from several file in run folder.

    :param run_folder: Path to the run folder.
    :type run_folder: str
    :return: Run information.
    :rtype: dict
    """
    # RunInfo
    run_info_path = os.path.join(run_folder, "RunInfo.xml")
    run_info = RunInfo(run_info_path)
    run = {
        "reads_phases": run_info.reads_phases,
        "flowcell": run_info.flowcell,
        "instrument": run_info.instrument,
        "info": run_info.run,
        "samples": None
    }
    run["info"]["end_date"] = None
    # RTAComplete
    rta_complete_path = os.path.join(run_folder, "RTAComplete.txt")
    if os.path.exists(rta_complete_path):
        rta_complete = RTAComplete(rta_complete_path)
        run["info"]["end_date"] = rta_complete.end_date
    # SampleSheet
    samplesheet_path = os.path.join(run_folder, "SampleSheet.csv")
    if os.path.exists(samplesheet_path):
        run["samples"] = SampleSheetIO(samplesheet_path).samples
    return run
