#!/usr/bin/env python3
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
__version__ = '2.2.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import time
import glob
import warnings
import argparse
import subprocess
import importlib.util

# Load Illumina
bin_dir = os.path.abspath(os.path.dirname(__file__))
illumina_lib = os.path.join(os.path.dirname(bin_dir), "lib", "illumina.py")
spec = importlib.util.spec_from_file_location("illumina", illumina_lib)
illumina = importlib.util.module_from_spec(spec)
spec.loader.exec_module(illumina)


########################################################################
#
# FUNCTIONS
#
########################################################################
def exec_cmd(cmd):
    """
    @summary: Display and submit the command.
    @param cmd: [list] The command to execute.
    """
    print(" ".join(cmd))
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL)


def getProtocol(in_spl_folder):
    """
    @summary: Returns the analysis protocol (design and Illumina workflow) declarated in samplesheet.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @return: [dict] The analysis protocol.
    @warning: For Amplicon - DS the manifests names must be <DESIGN>_A.txt and <DESIGN>_B.txt.
    """
    protocol = {"workflow": None, "design": None}
    samplesheet = illumina.SampleSheetIO(os.path.join(in_spl_folder, "SampleSheet.csv"))
    protocol["workflow"] = samplesheet.header["Workflow"]
    if samplesheet.header["Application"] == "Amplicon - DS" or samplesheet.header["Workflow"] == "Amplicon - DS":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        protocol["design"] = manifest_A.split("_A.txt")[0]
    elif samplesheet.header["Workflow"] == "Enrichment":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        protocol["design"] = manifest_A.split(".txt")[0]
    return protocol


def getWorkflows(protocol):
    """
    @summary: Returns the list of workflows to execute on data.
    @param protocol: [str] The analysis protocol (design and Illumina workflow).
    @return: [list] The names of the workflows to execute.
    """
    wf = list()
    if protocol["workflow"] == "Amplicon - DS":
        wf = ["ADIVaR", "AmpliconDS"]
    return wf


def getRunCmd(workflow, in_spl_folder, out_run_folder):
    """
    @summary: Returns the command to launch the workflow.
    @param workflow: [str] The name of the workflow.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param out_run_folder: [str] Path to the workflow output folder.
    @return: [None/list] The command if the corresponding workflow exists and None otherwise.
    """
    cmd = None
    if workflow == "ADIVaR":
        cmd = getADIVaRCmd(in_spl_folder, out_run_folder, protocol["design"])
    elif workflow == "AmpliconDS":
        cmd = getADSACmd(in_spl_folder, out_run_folder, protocol["design"])
    return cmd


def getADSACmd(in_spl_folder, out_run_folder, design):
    """
    @summary: Returns the command to launch the variant annotation on MiSeq reporter outputs.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param design: [str] Path to the folder containing the files describing the amplicons.
    @return: [list] The command.
    """
    ressources_folder = args.design_path_pattern.replace("WF_NAME", "AmpliconDSAnnot")
    cmd = [
        args.workflow_path_pattern.replace("WF_NAME", "AmpliconDSAnnot"), "amplicondsannot",
        "--RNA-selection", os.path.join(ressources_folder, design, "reference_RNA.tsv"),
        "--assembly-version", genome["assembly"],
        "--filters", os.path.join(ressources_folder, "ampliDS_filters_wfAmpliconDS.json"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = os.path.join(ressources_folder, design, genome["assembly"] + "_chr", "pos_ctrl_expected.vcf")
    if os.path.exists(positive_ctrl_ref):
        if design == "INCa_V1":
            cmd.extend(["--pos-ctrl-names", ".*[Hh][Oo][Rr][Ii].*"])
        else:
            cmd.extend(["--pos-ctrl-names", ".*[Tt][Ee][Mm][Oo][Ii][Nn].*"])
        cmd.extend(["--pos-ctrl-expected", positive_ctrl_ref])
    return cmd


def getADIVaRCmd(in_spl_folder, out_run_folder, design):
    """
    @summary: Returns the command to launch the amplicon double strand workflow.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param out_run_folder: [str] Path to the workflow output folder.
    @param design: [str] The name of the folder containing the files describing the amplicons.
    @return: [list] The command.
    """
    ressources_folder = args.design_path_pattern.replace("WF_NAME", "AmpliconDSAnnot")
    cmd = [
        args.workflow_path_pattern.replace("WF_NAME", "ADIVaR"), "adivar",
        "--R1-end-adapter", os.path.join(ressources_folder, "adapters", "Illumina_3prim_adapter.fasta"),
        "--R2-end-adapter", os.path.join(ressources_folder, "adapters", "Illumina_5prim_adapter_rvc.fasta"),
        "--libA", "folderA=" + os.path.join(ressources_folder, design, genome["assembly"], "libA"),
        "--libB", "folderB=" + os.path.join(ressources_folder, design, genome["assembly"], "libB"),
        "--RNA-selection", os.path.join(ressources_folder, design, "reference_RNA.tsv"),
        "--assembly-version", genome["assembly"],
        "--genome-seq", genome["sequences"],
        "--filters", os.path.join(ressources_folder, "ampliDS_filters.json"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = os.path.join(ressources_folder, design, genome["assembly"], "pos_ctrl_expected.vcf")
    if os.path.exists(positive_ctrl_ref):
        if design == "INCa_V1":
            cmd.extend(["--pos-ctrl-names", ".*[Hh][Oo][Rr][Ii].*"])
        else:
            cmd.extend(["--pos-ctrl-names", ".*[Tt][Ee][Mm][Oo][Ii][Nn].*"])
        cmd.extend(["--pos-ctrl-expected", positive_ctrl_ref])
    return cmd


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="This script uses an infinite loop to listen the Illumina's sequencer output folder and launch the appropriate workflow when a run is ended.")
    parser.add_argument('-r', '--roll-time', type=int, default=(60*20), help="The time between each sequencer output folder evaluation (in seconds). [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-l', '--listened-folder', required=True, help="The sequencer output folder.")
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-s', '--storage-folder', help='Path to the storage folder. The run folder will be moved to this folder after process.')
    group_output.add_argument('-a', '--analysis-folder', required=True, help='Path to the workflows output folder. Each workflow will create a sub-folder in this parent folder (example: for the analysis foler "170101_run01" the result can be two sub-folders 170101_run01/ADIVaR and 170101_run01/AmpliconDS).')
    group_output = parser.add_argument_group('Ressources')  # Ressources
    group_output.add_argument('-w', '--workflow-path-pattern', required=True, help='The path to the jflow client of the workflows. The tag WF_NAME is used as placeholder to adapt the command at the workflow. Example: /usr/local/bioinfo/WF_NAME/current/bin/jflow_cli.py.')
    group_output.add_argument('-d', '--design-path-pattern', required=True, help='The path to the folder containing the designs definitions folders. The tag WF_NAME is used as placeholder to adapt the command at the workflow. Example: /usr/local/bioinfo/WF_NAME/current/ressources.')
    group_output.add_argument('-g', '--genome', required=True, help='The assembly version and the path to the reference genome used. The two information are separated by ":" like in the following example: GRCh37:/bank/Homo_sapiens/GRCH37/dna.fasta.')
    args = parser.parse_args()
    genome = {
        'assembly': args.genome.split(":")[0],
        'sequences': args.genome.split(":")[1]
    }

    # Process
    while True:
        for filename in os.listdir(args.listened_folder):
            in_run_folder = os.path.join(args.listened_folder, filename)
            if os.path.isdir(in_run_folder):
                out_run_folder = os.path.join(args.analysis_folder, filename)
                if os.path.exists(os.path.join(in_run_folder, "CompletedJobInfo.xml")):  # The run is ended
                    in_basecalls_folder = os.path.join(in_run_folder, "Data", "Intensities", "BaseCalls")
                    protocol = getProtocol(in_basecalls_folder)
                    # Process Anapath's workflows
                    if not os.path.exists(out_run_folder):  #################################################### To remove
                        if not os.path.exists(out_run_folder):
                            os.mkdir(out_run_folder)
                        # Analysis workflows
                        for curr_wf in getWorkflows(protocol):
                            out_wf_folder = os.path.join(out_run_folder, curr_wf)
                            if os.path.exists(out_wf_folder):
                                warnings.warn('The workflow "{}" has already be processed for run "{}".'.format(curr_wf, in_run_folder))
                            else:
                                os.mkdir(out_wf_folder)
                                cmd_analysis = getRunCmd(curr_wf, in_basecalls_folder, out_wf_folder)
                                exec_cmd(cmd_analysis)
                    # Copy raw data if it not already exist
                    if args.storage_folder is not None:
                        out_run_storage = os.path.join(args.storage_folder, filename)
                        if not os.path.exists(out_run_storage):  # The run folder has not been copied in storage folder
                            cmd_copy_raw = ["rsync", "--recursive", "--perms", "--times", in_run_folder + os.sep, out_run_storage]
                            exec_cmd(cmd_copy_raw)
                            # Remove copy of the Illumina's workflows
                            storage_basecalls_folder = os.path.join(out_run_storage, "Data", "Intensities", "BaseCalls")
                            alignments_files = glob.glob(os.path.join(storage_basecalls_folder, "Alignment*"))
                            for curr_file in alignments_files:
                                if os.path.isdir(curr_file):
                                    cmd_rm_analysis = ["rm", "-r", curr_file]
                                    exec_cmd(cmd_rm_analysis)
                    # Copy Illumina's workflows if they not already exist
                    if protocol["workflow"] != "GenerateFASTQ":  # Illumina reporter has processed an analysis
                        # Get workflow name
                        analysis_basename = protocol["workflow"].replace(" ", "_")
                        if analysis_basename == "Amplicon_-_DS":
                            analysis_basename = "AmpliconDS"
                        # Process workflow results
                        in_basecalls_folder = os.path.join(in_run_folder, "Data", "Intensities", "BaseCalls")
                        in_wf_folder_prefix = os.path.join(in_basecalls_folder, "Alignment")
                        alignments_files = glob.glob(in_wf_folder_prefix + "*")
                        for in_wf_folder in alignments_files:
                            if os.path.isdir(in_wf_folder):
                                out_wf_folder = os.path.join(out_run_folder, analysis_basename)
                                out_wf_folder += in_wf_folder.replace(in_wf_folder_prefix, "")  # Add workflow suffix (the index of the algnment folder)
                                log_file = os.path.join(out_wf_folder, "CompletedJobInfo.xml")
                                if not os.path.exists(log_file):
                                    cmd_copy_analysis = ["rsync", "--recursive", "--perms", "--times", in_wf_folder + os.sep, out_wf_folder]
                                    exec_cmd(cmd_copy_analysis)
        time.sleep(args.roll_time)
