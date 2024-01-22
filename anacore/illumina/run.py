# -*- coding: utf-8 -*-
"""
Classes and functions for reading Illumina's files on run: runInfo, runParameters, ... .

:Code example:

    Check if run is ended

    .. highlight:: python
    .. code-block:: python

        from anacore.illumina.run import Run

        run = Run("my_run/"):
        print("Run steps:")
        print("", "1- sequencing ended:", run.isSequenced())
        print("", "2- copy ended:", run.isCopied())
        print("", "3- pos-processing ended:", run.isEnded())

        # Result>
        # Run steps:
        #   1- sequencing ended: True
        #   2- copy ended: False
        #   3- pos-processing ended: False
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'

from anacore.illumina.base import getPlatformFromSerialNumber
from anacore.illumina.samplesheet import SampleSheetFactory
import datetime
import os
import re
import xml.etree.ElementTree as ET


class CompletedJobInfo(object):
    """Reader for CompletedJobInfo.xml."""

    def __init__(self, path, date_format='%Y-%m-%dT%H:%M:%S'):
        """
        Build and return an instance of CompletedJobInfo.

        :param path: Path to CompletedJobInfo.xml (format: XML).
        :type path: str
        :param date_format: Date format on sequencer.
        :type date_format: str
        :return: The new instance.
        :rtype: CompletedJobInfo
        """
        self.filepath = path
        self.date_format = date_format
        self.start_datetime = None
        self.end_datetime = None
        self.version = None
        self.workflow_name = None
        self.parameters = None
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
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
    """
    Return dict from XML node.

    :param node: Node to convert.
    :type node: xml.etree.ElementTree.Element
    :return: Dict corresponding to XML node.
    :rtype: dict
    """
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


def findAny(node, searches):
    """
    Return the first subelement matching any term of searches.

    :param node: Node to convert.
    :type node: xml.etree.ElementTree.Element
    :param searches: List accepted element tag or XPath. Order determines the first match.
    :type searche: str
    :return: Subelement matching any term of searches.
    :rtype: xml.etree.ElementTree.Element
    """
    for curr in searches:
        match = node.find(curr)
        if match is not None:
            return match
    return None


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
        run["samples"] = SampleSheetFactory.get(samplesheet_path).samples
    return run


class RTAComplete(object):
    """Reader for RTAComplete.txt."""

    def __init__(self, path):
        """
        Build and return an instance of RTAComplete.

        :param path: Path to RTAComplete.txt (format: txt).
        :type path: str
        :return: The new instance.
        :rtype: RTAComplete
        """
        self.filepath = path
        self.end_date = None
        self.RTA_version = None
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes. If file come from NovaSeq end_date is retrieve from file last modification date."""
        with open(self.filepath) as FH_in:
            content = FH_in.read().strip()
        if content == "":  # is empty (NovaSeq)
            self.RTA_version = None
            self.end_date = datetime.datetime.fromtimestamp(os.path.getmtime(self.filepath))
        else:
            match = re.fullmatch(r"(.+/.+/.+,.+),Illumina RTA (.+)", content)
            if match:  # 11/1/2017,15:11:43.174,Illumina RTA 1.18.54
                date_str, self.RTA_version = match.groups()
                if "." in date_str and ".":
                    date_str = date_str.split(".")[0]
                self.end_date = datetime.datetime.strptime(date_str, '%m/%d/%Y,%H:%M:%S')
            else:
                match = re.fullmatch(r"RTA (.+) completed on (.+/.+/.+ .+:.+:.+ ..)", content)
                if match:  # RTA 2.4.11 completed on 11/14/2019 4:56:45 AM
                    self.RTA_version, date_str = match.groups()
                    self.end_date = datetime.datetime.strptime(date_str, '%m/%d/%Y %I:%M:%S %p')
                else:
                    raise Exception('"{}" in {} cannot be parsed by {}.'.format(content, self.filepath, self.__class__.__name__))


class Run:
    """Read and provide getters on Run from run folder."""

    def __init__(self, path):
        """
        Build and return an instance of Run.

        :param path: Path to the run folder.
        :type path: str
        :return: The new instance.
        :rtype: Run
        """
        self.path = path

    def isCopied(self):
        """
        Return True if copy from temp to output folder is ended by instrument. Copy step is followed post processing step on instrument.

        :return: True if copy from temp to output folder is ended by instrument.
        :rtype: boolean
        """
        is_completed = False
        if self.isSequenced():
            run_params = self.parameters
            completed_copy_file = os.path.join(self.path, "CopyComplete.txt")
            if run_params.instrument["platform"] == "HiSeq":
                is_completed = True
            elif run_params.instrument["platform"] == "MiSeq":
                mcs_version = run_params.software["CS"]
                if int(mcs_version.split(".")[0]) < 4:
                    is_completed = True
                else:  # The universal copy service replace the previous system since MCS v4
                    if os.path.exists(completed_copy_file):
                        is_completed = True
            elif run_params.instrument["platform"] == "NextSeq":
                ncs_version = run_params.software["CS"]
                if int(ncs_version.split(".")[0]) < 4:
                    is_completed = True
                else:  # The universal copy service replace the previous system since NCS v4
                    if os.path.exists(completed_copy_file):
                        is_completed = True
            elif os.path.exists(completed_copy_file):
                is_completed = True
        return is_completed

    def isEnded(self):
        """
        Return True if all instrument steps are ended (sequencing, copy and post-process).

        :return: True if all instrument steps are ended (sequencing, copy and post-process).
        :rtype: boolean
        """
        is_completed = False
        if self.isCopied():  # Copy is ended
            completed_illumina_analysis_file = os.path.join(self.path, "CompletedJobInfo.xml")
            if self.parameters.post_process is None or os.path.exists(completed_illumina_analysis_file):  # The sequencer does not produce post-processing or they are ended
                is_completed = True
        return is_completed

    def isSequenced(self):
        """
        Return True if sequencing step is ended on instrument. Sequencing step is followed by copy and post processing step on instrument.

        :return: True if sequencing step is ended on instrument.
        :rtype: boolean
        """
        rta_complete_file = os.path.join(self.path, "RTAComplete.txt")
        return os.path.exists(rta_complete_file)

    @property
    def parameters(self):
        """
        Return run parameters from [rR]unParameters.xml.

        :return: Run parameters from RunParameters.xml.
        :rtype: anacore.illumina.RunParameters
        """
        run_parameters = os.path.join(self.path, "RunParameters.xml")
        if not os.path.exists(run_parameters):
            run_parameters = os.path.join(self.path, "runParameters.xml")
            if not os.path.exists(run_parameters):
                raise Exception("RunParameters cannot be found in {}.".format(self.path))
        return RunParameters(run_parameters)


class RunInfo(object):
    """Reader for RunInfo.xml (example: see my_run/RunInfo.xml)."""

    def __init__(self, path):
        """
        Build and return an instance of RunInfo.

        :param path: Path to RunInfo.xml (format: XML).
        :type path: str
        :return: The new instance.
        :rtype: RunInfo
        """
        self.filepath = path
        self.flowcell = None  # {'id': 'HN33VBGX5', 'layout': {'LaneCount': '4', 'SurfaceCount': '2', 'SwathCount': '3', 'TileCount': '12', 'SectionPerLane': '3', 'LanePerSection': '2'}}
        self.instrument = None  # {'id': 'NS500523', 'platform': 'NextSeq'}
        self.reads_phases = None  # [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}]
        self.run = None  # {'number': '133', 'id': '180406_NS500523_0133_AHN33VBGX5', 'start_date': datetime.datetime(2018, 4, 6, 0, 0)}
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
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
            "platform": getPlatformFromSerialNumber(serial_number)
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
    """Reader for runParameters.xml (example: see my_run/RunParameters.xml)."""

    def __init__(self, path):
        """
        Build and return an instance of RunParameters.

        :param path: Path to runParameters.xml (format: XML).
        :type path: str
        :return: The new instance.
        :rtype: RunParameters
        """
        self.filepath = path
        self.kit = None  # {"flowcell_id": "H14U5BGXX", "reagent_kit_id": "NS2044178-REAGT"}
        self.instrument = None  # {"id": "NS500413", "platform": "NextSeq"}
        self.reads_phases = None  # [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}]
        self.run = None  # {"number": "0033", "id": "141107_NS500413_0033_H14U5BGXX", "start_date": datetime}
        self.post_process = None  # "GenerateFASTQ"
        self.software = None  # {"RTA": "2.11.3", "CS": "4.0.1.41"}
        self._parse()

    def _parse(self):
        """Read file content and store information on the instance's attributes."""
        # Retrieve information by section
        tree = ET.parse(self.filepath)
        root = tree.getroot()
        if {"Version", "Setup"} == {child.tag for child in root}:  # For HiSeq information are all in <Setup> under <RunParameters>
            root = [child for child in root if child.tag == "Setup"][0]
        # Process information
        self.instrument = self._getInstrumentFromRoot(root)
        self.reads_phases = self._getReadsFromRoot(root)
        self.run = self._getRunFromRoot(root)
        self.kit = self._getKitFromRoot(root)
        self.post_process = self._getPostProcessFromRoot(root)
        self.software = self._getSoftwareFromRoot(root)

    def _getReadsFrom(self, subtree, prefix=""):
        reads = list()
        indices = list()
        for child in subtree:
            if child.tag.startswith(prefix + "Read"):
                reads.append({
                    "is_index": False,
                    "nb_cycles": int(child.text)
                })
            elif child.tag.startswith(prefix + "Index"):
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

    def _getReadsFromRoot(self, root):
        reads_phases = None
        reads_subtree = root.find("Reads")
        if reads_subtree is not None:
            reads_phases = self._getReadsFromReads(reads_subtree)
        else:
            setup_subtree = root.find("Setup")
            if setup_subtree is not None:
                reads_phases = self._getReadsFrom(root.find("Setup"))
            else:
                reads_phases = self._getReadsFrom(root, "Planned")  # Nova style
        return reads_phases

    def _getReadsFromReads(self, subtree):
        reads = list()
        for child in subtree:
            reads.append({
                "is_index": (child.get("IsIndexedRead") == "Y"),
                "nb_cycles": int(child.get("NumCycles"))
            })
        return reads

    def _getInstrumentFromRoot(self, root):
        serial_number = findAny(root, ("ScannerID", "InstrumentID", "InstrumentName")).text
        return {
            "id": serial_number,
            "platform": getPlatformFromSerialNumber(serial_number)
        }

    def _getRunFromRoot(self, root):
        run_number = findAny(root, ("RunNumber", "ScanNumber"))
        run_number = run_number.text.lstrip("0")
        run_id = findAny(root, ("RunID", "RunId")).text
        return {
            "number": run_number,
            "id": run_id,
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
            if flowcell_markup is not None and flowcell_markup.find("SerialNumber") is not None:
                flowcell_id = flowcell_markup.find("SerialNumber").text
            else:
                flowcell_id = root.find("Barcode")  # HiSeq style or MiSeq with invalid RFID markup
                if flowcell_id is not None:
                    flowcell_id = flowcell_id.text
                else:
                    flowcell_id = self.run["id"].rsplit("_", 1)[1]
        # Reagent kit ID
        reagent_kit_id = findAny(root, ("ReagentKitBarcode", "ReagentKitSerial"))
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
        rta_version = findAny(root, ("RTAVersion", "RtaVersion")).text
        if rta_version.startswith("v"):
            rta_version = rta_version[1:]
        return {
            "RTA": rta_version,
            "CS": cs_root_markup.find("ApplicationVersion").text
        }
