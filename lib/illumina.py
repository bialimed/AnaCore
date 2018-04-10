#
# Copyright (C) 2017 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.8.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import re
import glob
import datetime
import xml.etree.ElementTree as ET


class SampleSheetIO(object):
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
                if line.strip() in ["[Header]", "[Manifests]", "[Reads]", "[Settings]", "[Data]"]:
                    section_title = line.strip()[1:-1]
                    sections_by_title[section_title] = list()
                elif line.strip() == "":
                    section_title = None
                else:
                    sections_by_title[section_title].append(line.strip())
        # Process information
        self.samples = self._getSamplesFromData(sections_by_title["Data"])
        self.header = self._getInfoFromSection(sections_by_title["Header"])
        self.manifests = self._getInfoFromSection(sections_by_title["Manifests"]) if "Manifests" in sections_by_title else dict()
        self.run = self._getRunFromHeaderAndReads(sections_by_title["Header"], sections_by_title["Reads"])

    def _getSamplesFromData(self, data_section):
        samples = list()
        data_titles = [field.strip() for field in data_section[0].split(",")]
        for line in data_section[1:]:
            samples.append(
                {data_titles[idx]: field.strip() for idx, field in enumerate(line.split(","))}
            )
        return samples

    def _getInfoFromSection(self, section):
        info = dict()
        for line in section:
            key, value = [field.strip() for field in line.split(",", 1)]
            info[key] = value
        return info

    def _getRunFromHeaderAndReads(self, header_section, reads_section):
        run = dict()
        # Header
        for line in header_section:
            key, value = line.split(",")
            run[key] = value
        # Reads
        read_idx = 1
        for line in reads_section:
            run["nb_cycles_R" + str(read_idx)] = int(line.strip())
            read_idx += 1
        return run


class ADSSampleSheetIO(SampleSheetIO):
    """
    @summary: Manage SampleSheet designed for AmpliconDS analysis.
    """

    def _getSamplesFromData(self, data_section):
        samples = super()._getSamplesFromData(data_section)
        for spl_idx, spl in enumerate(samples):
            spl["Sample_Basename"] = getIlluminaName(spl["Sample_Name"])
            spl["Library_Basename"] = spl["Sample_Basename"] + "_S" + str(spl_idx + 1)
            spl["Library_Name"] = spl["Library_Basename"]
        return samples


    def filterPanels(self, selected_manifests):
        """
        @summary: Filters samples and manifests on their panel.
        @param selected_manifests: [list] The manifests corresponding to the selected panel.
        """
        # Get kept
        new_manifests = {}
        new_by_old_id = {}
        new_idx = 0
        for manifest_id in sorted(self.manifests):
            manifest_file = self.manifests[manifest_id]
            if manifest_file in selected_manifests:
                print(manifest_file, selected_manifests)
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
        @summary: Returns by sample name or library name all the files in directory corresponding to the subject.
        @param directory: [str] The path to the directory where the files names are evaluated.
        @param end_pattern: [str] File pattern added at the end of each sample basename or library basename to select files.
        @param subject: [str] "library" or "sample" to select files corresponding to libraries or corresponding to samples.
        @return: [dict] By subject the pathes of the corresponding files.
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
        @summary: Adds the files in directory corresponding to the subject. These information are keyed by provided tag for each element of self.samples.
        @param tag: [str] The attribute used in self.sample to store the files pathes.
        @param directory: [str] The path to the directory where the files names are evaluated.
        @param end_pattern: [str] File pattern added at the end of each sample basename or library basename to select files.
        @param subject: [str] "library" or "sample" to select files corresponding to libraries or corresponding to samples.
        """
        subject_tag = subject.capitalize() + "_Name"
        files_by_spl = self.findSplFiles(directory, end_pattern, subject)
        for spl in self.samples:
            spl[tag] = files_by_spl[spl[subject_tag]]


class RTAComplete(object):
    def __init__(self, path, date_format='%m/%d/%Y,%H:%M:%S'):
        self.filepath = path
        self.end_date = None
        self.RTA_version = None
        self.date_format = date_format
        self._parse()

    def _parse(self):
        with open(self.filepath) as FH_in:
            fields = FH_in.read().rsplit(",", 1)  # 11/1/2017,15:11:43.174,Illumina RTA 1.18.54
            if "." in fields[0] and "." not in self.date_format:
                fields[0] = fields[0].split(".")[0]
            self.end_date = datetime.datetime.strptime(fields[0], self.date_format)
            self.RTA_version = fields[1].rsplit(" ", 1)[0]


class RunParameters(object):
    PLATFORMS = ["MiniSeq", "MiSeq", "NextSeq", "HiSeq", "NovaSeq"]

    def __init__(self, path):
        self.filepath = path
        self.kit = None
        self.instrument = None
        self.reads = None
        self.run = None
        self._parse()

    def _parse(self):
        # Retrieve information by section
        tree = ET.parse(self.filepath)
        root = tree.getroot()
        # Process information
        self.instrument = self._getInstrumentFromRoot(root)
        reads_subtree = root.find("Reads")
        if reads_subtree is not None:
            self.reads = self._getReadsFromSection(reads_subtree)
        else:
            self.reads = self._getReadsFromSetup(root.find("Setup"))
        self.run = self._getRunFromRoot(root)
        self.kit = self._getKitFromRoot(root)

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
        # Get instrument platform
        platform = None
        max_nb_occur = 0
        raw_file = str(ET.tostring(root))
        for curr_platform in RunParameters.PLATFORMS:
            nb_occur = raw_file.count(curr_platform)
            if nb_occur > max_nb_occur:
                platform = curr_platform
        # Get instrument serial
        serial_number = root.find("ScannerID")
        if serial_number is None:
            serial_number = root.find("InstrumentID")
        serial_number = serial_number.text
        # Return
        return {
            "id": serial_number,
            "platform": platform
        }

    def _getRunFromRoot(self, root):
        run_number = root.find("RunNumber")
        if run_number is None:
            run_number = root.find("ScanNumber")
        run_number = run_number.text
        return {
            "number": run_number,
            "id": root.find("RunID").text,
            "start_date": datetime.datetime.strptime(root.find("RunStartDate").text, '%y%m%d')
        }

    def _getKitFromRoot(self, root):
        # Flowcell ID
        flowcell_id = self.run["id"].rsplit("_", 1)[1]  # HiSeq: <RunID>141110_D00267_0193_BHAMVRADXX</RunID> ; NextSeq: <FlowCellSerial> ; MiSeq: <FlowcellRFIDTag> > <SerialNumber>
        # Reagent kit ID
        reagent_kit_id = root.find("ReagentKitBarcode")
        if reagent_kit_id is None:
            reagent_kit_id = root.find("ReagentKitSerial")
        reagent_kit_id = reagent_kit_id.text
        return {
            "flowcell_id": flowcell_id,
            "reagent_kit_id": reagent_kit_id
        }


class CompletedJobInfo(object):
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
    @sumary: Returns sample name used by Illumina in filename.
    @param name: [str] The name provided to Illumina process (for example in samplesheet).
    @return: [str] The sample name used by Illumina as part of filename.
    """
    return name.replace("_", "-").replace(" ", "-").replace(".", "-").replace("+", "")


def getLibNameFromReadsPath(seq_path):
    """
    @sumary: Returns library name from the path of the sequences file.
    @param seq_path: [str] The path of the sequences file.
    @return: [str] The library name.
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
