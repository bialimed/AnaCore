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
__version__ = '2.13.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import yaml
import time
import glob
import smtplib
import logging
import argparse
import datetime
import traceback
import subprocess
from email.mime.text import MIMEText

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)  # Load libraries for current script but not in bprocess

from anacore.illumina import SampleSheetIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class Design:
    """Class to describe design and his resources."""

    def __init__(self, name, pos_ctrl_re=None, neg_ctrl_re=None, resources_folder=None):
        """
        Build and return an instance of Design.

        :param name: Name/ID of pannel.
        :type name: str
        :param pos_ctrl_re: Regexp used to indetify positive control samples from a list of files pathes.
        :type pos_ctrl_re: str
        :param neg_ctrl_re: Regexp used to indetify negative control samples from a list of files pathes.
        :type neg_ctrl_re: str
        :param resources_folder: Path to the folder containing one folder by design. The name of the subfolders is the name of designs.
        :type resources_folder: str
        :return: The new instance.
        :rtype: Design
        """
        self.name = name
        self.pos_ctrl_re = pos_ctrl_re
        self.neg_ctrl_re = neg_ctrl_re
        self.resources_folder = resources_folder

    def getAssemblyFolder(self, assembly):
        """
        Return the resource folder for the design in the selected assembly.

        :param assembly: The selected assembly.
        :type assembly: str
        :return: Path to the resource folder.
        :rtype: str
        """
        return os.path.join(self.resources_folder, self.name, assembly)

    def getPosCtrlRef(self, assembly):
        """
        Return the file containing the reference values for the positive control (for example the VCF containing the expected variants).

        :param assembly: The selected assembly.
        :type assembly: str
        :return: Path to the reference values file.
        :rtype: str
        """
        return os.path.join(self.getAssemblyFolder(assembly), "pos_ctrl_expected.vcf")

    def getVarNoise(self, assembly):
        """
        Return the file containing the list of constitutive variants with their maximum AF.

        :param assembly: The selected assembly.
        :type assembly: str
        :return: Path to the reference values file.
        :rtype: str
        """
        return os.path.join(self.getAssemblyFolder(assembly), "variants_noise.tsv")

    def getSelectedRef(self, ref_type):
        """
        Return the file containing the reference selected (for example the list of RNA where the variants are kept).

        :param ref_type: The type of reference.
        :type ref_type: str
        :return: Path to the file.
        :rtype: str
        """
        return os.path.join(self.resources_folder, self.name, "reference_" + ref_type + ".tsv")

    def getFilters(self, wf_name):
        """
        Return the file containing the filters to apply on the workflow results.

        :param wf_name: The workflow name.
        :type wf_name: str
        :return: Path to the filters description file.
        :rtype: str
        """
        return os.path.join(self.resources_folder, self.name, wf_name + "_filters.json")


class Cmd:
    """Class to manage command line(s)."""

    def __init__(self, cmd, is_shell=False):
        """
        Build and return an instance of Cmd.

        :param cmd: The command to execute.
        :type cmd: list
        :param is_shell: If shell is True, the specified command will be executed through the shell.
        :type is_shell: boolean
        :return: The new instance.
        :rtype: Cmd
        """
        self.cmd = cmd
        self.is_shell = is_shell

    def run(self, log):
        """
        @summary: Log and submit the command.
        @param log: [Logger] The logger of the script.
        """
        log.info("Start sub-command: " + " ".join([str(elt) for elt in self.cmd]))
        if self.is_shell:
            subprocess.check_call(" ".join([str(elt) for elt in self.cmd]), stdout=subprocess.DEVNULL, shell=True)
        else:
            subprocess.check_call(self.cmd, stdout=subprocess.DEVNULL)
        log.info("End sub-command")


def getProtocol(in_spl_folder):
    """
    Return the analysis protocol (design and Illumina workflow) declarated in samplesheet.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :return: The analysis protocol. This dictionary contains 2 keys: workflow and design.
    :rtype: dict
    :warning: For Amplicon - DS the manifests names must be <DESIGN>_A.txt and <DESIGN>_B.txt.
    """
    protocol = {"workflow": None, "design": None}
    samplesheet = SampleSheetIO(os.path.join(in_spl_folder, "SampleSheet.csv"))
    protocol["workflow"] = samplesheet.header["Workflow"]
    if samplesheet.header["Application"] == "Amplicon - DS" or samplesheet.header["Workflow"] == "Amplicon - DS":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        design_name = manifest_A.split("_A.txt")[0]
        pos_ctrl_re = "[Tt][Ee][Mm].*"
        protocol["design"] = Design(design_name, pos_ctrl_re, "[Nn][Tt][Cc].*")
    elif samplesheet.header["Workflow"] == "Enrichment":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        protocol["design"] = Design(
            manifest_A.split(".txt")[0], None, "[Nn][Tt][Cc].*"
        )
    return protocol


def getWfConfig(cfg_path):
    """
    Return the workflows configuration file (format: YAML).

    :param cfg_path: Path to the configuration file (format: YAML).
    :type cfg_path: str
    :return: The analysis protocol. This dictionary contains 2 keys: workflow and design.
    :rtype: dict

    Example of configuration file:
        default:
            reference:
                sequences: /work/fescudie/bank/Homo_sapiens/DNA/GRCh37_Ensembl75_std/without_contig/Homo_sapiens.GRCh37.75.dna.woutContigs.fa
                assembly: GRCh37
                with_chr: False
            app_dir: /save/fescudie/IUCT/WF_NAME/current
            adapters_folder: /save/fescudie/IUCT/WF_NAME/current/resources/adapters
            designs_folder: /save/fescudie/IUCT/WF_NAME/current/resources/designs
        workflows:
            MIAmS:
                apply_on:
                    protocol: Amplicon - DS
                    designs: [TSCA_INCa_V2]  # None: without design ; []: all designs ; [A, B]: with design A or B
                app_dir:
                adapters_folder:
                designs_folder:
                reference:
                    sequences:
                    assembly:
                    with_chr: False  # Mandatory for mSINGS
    """
    cfg = None
    # Get parameters
    with open(cfg_path, 'r') as FH_cfg:
        cfg = yaml.load(FH_cfg)
    # Replace empty by default values
    for wf_name, wf_param in cfg["workflows"].items():
        # Application
        if wf_param["app_dir"] is None:
            wf_param["app_dir"] = cfg["default"]["app_dir"].replace("WF_NAME", wf_name)
        if "adapters_folder" in wf_param and wf_param["adapters_folder"] is None:
            wf_param["adapters_folder"] = cfg["default"]["adapters_folder"].replace("WF_NAME", wf_name)
        if "designs_folder" in wf_param and wf_param["designs_folder"] is None:
            wf_param["designs_folder"] = cfg["default"]["designs_folder"].replace("WF_NAME", wf_name)
        # Reference
        if "reference" in wf_param:
            if "sequences" in wf_param["reference"] and wf_param["reference"]["sequences"] is None:
                wf_param["reference"]["sequences"] = cfg["default"]["reference"]["sequences"]
            if "assembly" in wf_param["reference"] and wf_param["reference"]["assembly"] is None:
                wf_param["reference"]["assembly"] = cfg["default"]["reference"]["assembly"]
            if "with_chr" in wf_param["reference"] and wf_param["reference"]["with_chr"] is None:
                wf_param["reference"]["with_chr"] = cfg["default"]["reference"]["with_chr"]
    return cfg


def getWorkflows(protocol, cfg):
    """
    Return the list of workflows to execute on data.

    :param protocol: The analysis protocol (design and Illumina workflow).
    :type protocol: str
    :param cfg: Configuration for workflows.
    :type cfg: dict
    :return: The names of the workflows to execute.
    :rtype: list
    """
    expected_wf = list()
    for wf_name, wf_param in cfg["workflows"].items():
        if wf_param["apply_on"]["protocol"] == "." or protocol["workflow"] == wf_param["apply_on"]["protocol"]:
            if protocol["design"] is None and wf_param["apply_on"]["designs"] is None:
                expected_wf.append(wf_name)
            elif len(wf_param["apply_on"]["designs"]) == 0 or protocol["design"].name in wf_param["apply_on"]["designs"]:
                expected_wf.append(wf_name)
    return expected_wf


def getRunCmd(wf_name, in_spl_folder, out_run_folder, cfg):
    """
    Return the command to launch the workflow.

    :param wf_name: The name of the workflow.
    :type wf_name: str
    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param cfg: Configuration for workflows.
    :type cfg: dict
    :return: The command if the corresponding workflow exists and None otherwise.
    :rtype: None/Cmd
    """
    wf_cfg = cfg["workflows"][wf_name]
    fct_name = "get{}Cmd".format(wf_name)
    if fct_name not in globals():
        raise Exception("No {} has been implemented for the workflow {}.".format(fct_name, wf_name))
    cmd = globals()[fct_name](in_spl_folder, out_run_folder, protocol["design"], wf_cfg)
    return cmd


def getLoadCmd(design, raw_folder, out_by_wf):
    """
    Return the command to load run and workflows data in laboratory database.

    :param design: The object describing the design.
    :type design: Design
    :param raw_folder: Path to the run folder.
    :type raw_folder: str
    :param out_by_wf: The workflow output path by workflow name.
    :type out_by_wf: dict
    :return: The command to load data.
    :rtype: Cmd
    """
    cmd = [
        "lisImportRun.py",
        "--results-folder", "/Anapath/Illumina_Run_datas/Routine/database_files",  ################################
        "--input-raw", raw_folder
    ]
    # Control samples
    if design.neg_ctrl_re is not None:
        cmd.extend(["--negative-ctrl-regexp", design.neg_ctrl_re])
    if design.pos_ctrl_re is not None:
        cmd.extend(["--positive-ctrl-regexp", design.pos_ctrl_re])
    # Workflows data
    if len(out_by_wf) != 0:
        cmd.append("--input-workflows")
        for wf, wf_folder in out_by_wf.items():
            cmd.append("{}:{}".format(wf, wf_folder))
    return Cmd(cmd)


def getADSACmd(in_spl_folder, out_run_folder, design, wf_cfg):
    """
    Return the command to launch the variant annotation on MiSeq reporter outputs.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param design: The object describing the design.
    :type design: Design
    :param wf_cfg: Configuration for the workflow (subset of autoRunWf workflows's configutration).
    :type wf_cfg: dict
    :return: The command to run ADSA.
    :rtype: Cmd
    """
    design.resources_folder = wf_cfg["designs_folder"]
    assembly_design = wf_cfg["reference"]["assembly"] + "_chr" if wf_cfg["reference"]["with_chr"] else wf_cfg["reference"]["assembly"]
    app_exec = os.path.join(wf_cfg["app_dir"], "app", "bin", "jflow_cli.py")
    cmd = [
        app_exec, "amplicondsannot",
        "--RNA-selection", design.getSelectedRef("RNA"),
        "--assembly-version", wf_cfg["reference"]["assembly"],
        "--filters", design.getFilters("ADSA"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = design.getPosCtrlRef(assembly_design)
    if os.path.exists(positive_ctrl_ref):
        cmd.extend([
            "--pos-ctrl-names", design.pos_ctrl_re,
            "--pos-ctrl-expected", positive_ctrl_ref
        ])
    design.resources_folder = None
    return Cmd(cmd)


def getADIVaRCmd(in_spl_folder, out_run_folder, design, wf_cfg):
    """
    Return the command to launch the amplicon double strand workflow.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param design: The object describing the design.
    :type design: Design
    :param wf_cfg: Configuration for the workflow (subset of autoRunWf workflows's configutration).
    :type wf_cfg: dict
    :return: The command to run ADIVaR.
    :rtype: Cmd
    """
    design.resources_folder = wf_cfg["designs_folder"]
    assembly_design = wf_cfg["reference"]["assembly"] + "_chr" if wf_cfg["reference"]["with_chr"] else wf_cfg["reference"]["assembly"]
    app_exec = os.path.join(wf_cfg["app_dir"], "app", "bin", "jflow_cli.py")
    assembly_folder = design.getAssemblyFolder(assembly_design)
    cmd = [
        app_exec, "adivar",
        "--R1-end-adapter", os.path.join(wf_cfg["adapters_folder"], "Illumina_3prim_adapter.fasta"),
        "--R2-end-adapter", os.path.join(wf_cfg["adapters_folder"], "Illumina_5prim_adapter_rvc.fasta"),
        "--libA", "folderA=" + os.path.join(assembly_folder, "libA"),
        "--libB", "folderB=" + os.path.join(assembly_folder, "libB"),
        "--RNA-selection", design.getSelectedRef("RNA"),
        "--assembly-version", wf_cfg["reference"]["assembly"],
        "--genome-seq", wf_cfg["reference"]["sequences"],
        "--filters", design.getFilters("ADIVaR"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = design.getPosCtrlRef(assembly_design)
    if os.path.exists(positive_ctrl_ref):
        cmd.extend([
            "--pos-ctrl-names", design.pos_ctrl_re,
            "--pos-ctrl-expected", positive_ctrl_ref
        ])
    constit_variants = design.getVarNoise(assembly_design)
    if os.path.exists(constit_variants):
        cmd.extend([
            "--constit-variants", constit_variants
        ])
    design.resources_folder = None
    return Cmd(cmd)


def getMIAmSCmd(in_spl_folder, out_run_folder, design, wf_cfg):
    """
    Return the command to launch microsatellite instability detection by next-generation sequencing on amplicons.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param design: The object describing the design.
    :type design: Design
    :param wf_cfg: Configuration for the workflow (subset of autoRunWf workflows's configutration).
    :type wf_cfg: dict
    :return: The command to run MIAmS.
    :rtype: Cmd
    """
    app_env_bin = os.path.join(wf_cfg["app_dir"], "envs", "miniconda3", "bin")
    app_exec = os.path.join(wf_cfg["app_dir"], "jflow", "bin", "jflow_cli.py")
    design.resources_folder = wf_cfg["designs_folder"]
    assembly_folder = design.getAssemblyFolder(wf_cfg["reference"]["assembly"])
    cmd = [
        "unset", "PYTHONPATH", ";",
        "source", os.path.join(app_env_bin, "activate"), "MIAmS", "&&",
        app_exec, "miamstag",
        "--max-mismatch-ratio", 0.25,
        "--min-pair-overlap", 40,
        "--min-support-reads", 300,
        "--random-seed", 42,
        "--R1-end-adapter", os.path.join(wf_cfg["adapters_folder"], "Illumina_3prim_adapter.fasta"),
        "--R2-end-adapter", os.path.join(wf_cfg["adapters_folder"], "Illumina_5prim_adapter_rvc.fasta"),
        "--models", os.path.join(assembly_folder, "MSI", "models.json"),
        "--targets", os.path.join(assembly_folder, "MSI", "targets.bed"),
        "--intervals", os.path.join(assembly_folder, "MSI", "intervals.tsv"),
        "--baseline", os.path.join(assembly_folder, "MSI", "baseline.tsv"),
        "--genome-seq", wf_cfg["reference"]["sequences"],
        "--R1-pattern", "'" + os.path.join(in_spl_folder, "*R1_???.fastq.gz") + "'",
        "--R2-pattern", "'" + os.path.join(in_spl_folder, "*R2_???.fastq.gz") + "'",
        "--exclusion-pattern", "'" + "*Undetermined_S0_L???_R?_???.fastq.gz" + "'",
        "--output-dir", out_run_folder, "&&",
        "source", os.path.join(app_env_bin, "deactivate")
    ]
    design.resources_folder = None
    return Cmd(cmd, True)


def getArchiveCmd(profiles_folder, run, workflows):
    """
    Return the command to archivate run and workflows data.

    :param profiles_folder: Path to the folder containing profiles describing the data to archivate.
    :type profiles_folder: str
    :param run: Name of the run folder.
    :type run: str
    :param workflows: The names of the workflows to archivate.
    :type workflows: list
    :return: The command to archivate data.
    :rtype: Cmd
    """
    cmd = [
        "archive.py",
        "--placeholders", "'@YEAR={}' '@RUN={}'".format("20" + run[0:2], run),
        "--resources-profiles", os.path.join(profiles_folder, "Illumina", "*.yml")
    ]
    # Workflows data
    if len(workflows) != 0:
        for curr_wf in workflows:
            wf_profiles_folder = os.path.join(profiles_folder, "workflows", curr_wf)
            if os.path.exists(wf_profiles_folder):
                cmd.append(os.path.join(wf_profiles_folder, "*.yml"))
            else:
                print("The data for the workflow {} of the run {} cannot be archivated.".format(curr_wf, run))
    return Cmd(cmd, True)


def sendMail(smtp_adress, sender, recipients, subject, content):
    """
    Send execution log mail.

    :param smtp_adress: The SMTP server used to send mail.
    :type smtp_adress: str
    :param sender: Mail of the sender (it may be a not  existing adrsesses).
    :type sender: str
    :param recipients: Mails of the recipients.
    :type recipients: list
    :param subject: The mail subject.
    :type subject: str
    :param content: The mail content.
    :type content: str
    """
    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = ", ".join(recipients)
    smtp = smtplib.SMTP(smtp_adress)
    smtp.send_message(msg)
    smtp.quit()


def sendStatusMail(smtp_adress, sender, recipients, status, run_folder):
    """
    Send execution log mail.

    :param smtp_adress: The SMTP server used to send mail.
    :type smtp_adress: str
    :param sender: Mail of the sender (it may be a not  existing adrsesses).
    :type sender: str
    :param recipients: Mails of the recipients.
    :type recipients: list
    :param status: The status of the message: start or end or fail.
    :type status: str
    :param run_folder: Path to the run folder currently processed.
    :type run_folder: str
    """
    subject = "[NGS] {} sequencer post-processing".format(status.capitalize())
    content = "{}: The sequencer post-processing {} {}ed for {}.".format(
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        ("has" if status == "fail" else "is"),
        status.lower(),
        run_folder
    )
    sendMail(smtp_adress, sender, recipients, subject, content)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="This script uses an infinite loop to listen the Illumina's sequencer output folder and launch the appropriate workflow when a run is ended.")
    parser.add_argument('-r', '--roll-time', type=int, default=(60*20), help="The time between each sequencer output folder evaluation (in seconds). [Default: %(default)s]")
    parser.add_argument('-p', '--archivage-profiles-folder', help='Path to the folder containing profiles describing the data to archivate.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-l', '--listened-folder', required=True, help="The sequencer output folder.")
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-s', '--storage-folder', help='Path to the storage folder. The run folder will be moved to this folder after process.')
    group_output.add_argument('-a', '--analysis-folder', required=True, help='Path to the workflows output folder. Each workflow will create a sub-folder in this parent folder (example: for the analysis foler "170101_run01" the result can be two sub-folders 170101_run01/ADIVaR and 170101_run01/ADSA).')
    group_output = parser.add_argument_group('Workflows')  # Workflows configuration
    group_output.add_argument('-w', '--workflow-config', required=True, help='The path to the file containing configuration for workflows.')
    group_mail = parser.add_argument_group('Mail')  # Mail
    group_mail.add_argument('--mail-recipients', nargs='+', help='Mails of the recipients.')
    group_mail.add_argument('--mail-sender', help='Mail of the sender (it may be a not  existing adrsesses).')
    group_mail.add_argument('--mail-smtp', help='The SMTP server used to send mail.')
    args = parser.parse_args()
    analyses_cfg = getWfConfig(args.workflow_config)
    send_mail = False
    if args.mail_recipients or args.mail_sender or args.mail_smtp:  # At least one mail parameter is setted
        send_mail = True
        if not(args.mail_recipients and args.mail_sender and args.mail_smtp):  # At least one mail parameter is not setted
            raise Exception("To send mail all the mail parameters must be set.")

    # Logger
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger("autoRunWf")
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))
    log.info("Version: " + str(__version__))

    # Process
    while True:
        for run_id in os.listdir(args.listened_folder):
            in_run_folder = os.path.join(args.listened_folder, run_id)
            if os.path.isdir(in_run_folder):
                out_run_folder = os.path.join(args.analysis_folder, run_id)
                completed_illumina_file = os.path.join(in_run_folder, "CompletedJobInfo.xml")
                completed_analyses_file = os.path.join(out_run_folder, "processCompleted.txt")
                failed_analyses_file = os.path.join(out_run_folder, "processFailed.txt")
                if os.path.exists(completed_illumina_file) and not os.path.exists(completed_analyses_file) and not os.path.exists(failed_analyses_file):  # The run is ended
                    log.info("Start post-process on run {}.".format(run_id))
                    step = "start"
                    if send_mail:
                        sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "start", in_run_folder)
                    try:
                        if not os.path.exists(out_run_folder):
                            os.mkdir(out_run_folder)
                        in_basecalls_folder = os.path.join(in_run_folder, "Data", "Intensities", "BaseCalls")
                        protocol = getProtocol(in_basecalls_folder)
                        out_folder_by_wf = dict()
                        # Launch workflows
                        step = "Launch analysis workflows"
                        for curr_wf in getWorkflows(protocol, analyses_cfg):
                            out_wf_folder = os.path.join(out_run_folder, curr_wf)
                            out_folder_by_wf[curr_wf] = out_wf_folder
                            if os.path.exists(out_wf_folder):
                                log.warning('The workflow "{}" has already be processed for run "{}".'.format(curr_wf, in_run_folder))
                            else:
                                os.mkdir(out_wf_folder)
                                cmd_analysis = getRunCmd(curr_wf, in_basecalls_folder, out_wf_folder, analyses_cfg)
                                cmd_analysis.run(log)
                        # Copy raw data
                        step = "Copy raw"
                        if args.storage_folder is not None:
                            out_run_storage = os.path.join(args.storage_folder, run_id)
                            if not os.path.exists(out_run_storage):  # The run folder has not been copied in storage folder
                                cmd_copy_raw = ["rsync", "--recursive", "--perms", "--times", in_run_folder + os.sep, out_run_storage]
                                Cmd(cmd_copy_raw).run(log)
                                # Remove copy of the Illumina's workflows
                                storage_basecalls_folder = os.path.join(out_run_storage, "Data", "Intensities", "BaseCalls")
                                alignments_files = glob.glob(os.path.join(storage_basecalls_folder, "Alignment*"))
                                for curr_file in alignments_files:
                                    if os.path.isdir(curr_file):
                                        cmd_rm_analysis = ["rm", "-r", curr_file]
                                        Cmd(cmd_rm_analysis).run(log)
                        # Copy Illumina's workflows
                        step = "Copy Illumina's workflows"
                        if protocol["workflow"] != "GenerateFASTQ":  # Illumina reporter has processed an analysis
                            # Get workflow name
                            analysis_basename = protocol["workflow"].replace(" ", "_")
                            if analysis_basename == "Amplicon_-_DS":
                                analysis_basename = "ADSA"
                            # Process workflow results
                            in_basecalls_folder = os.path.join(in_run_folder, "Data", "Intensities", "BaseCalls")
                            in_wf_folder_prefix = os.path.join(in_basecalls_folder, "Alignment")
                            alignments_files = glob.glob(in_wf_folder_prefix + "*")
                            for in_wf_folder in alignments_files:
                                if os.path.isdir(in_wf_folder):
                                    out_wf_folder = os.path.join(out_run_folder, analysis_basename)
                                    out_wf_folder += in_wf_folder.replace(in_wf_folder_prefix, "")  # Add workflow suffix (the index of the algnment folder)
                                    if analysis_basename not in ["ADSA"]: out_folder_by_wf[analysis_basename] = out_wf_folder
                                    log_file = os.path.join(out_wf_folder, "CompletedJobInfo.xml")
                                    if not os.path.exists(log_file):  # The results have not been copied in workflow folder
                                        cmd_copy_analysis = ["rsync", "--recursive", "--perms", "--times", in_wf_folder + os.sep, out_wf_folder]
                                        Cmd(cmd_copy_analysis).run(log)
                        # Save in database
                        step = "Database storage"
                        out_run_storage = in_run_folder
                        if args.storage_folder is not None:
                            out_run_storage = os.path.join(args.storage_folder, run_id)
                        cmd_load = getLoadCmd(protocol["design"], out_run_storage, out_folder_by_wf)
                        Cmd(cmd_load).run(log)
                        # Log end process
                        with open(completed_analyses_file, "w") as FH:
                            FH.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                        if send_mail:
                            sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "end", in_run_folder)
                        # Send to archive
                        step = "Send to archive"
                        if args.archivage_profiles_folder is not None:
                            cmd_archivage = getArchiveCmd(args.archivage_profiles_folder, run_id, list(out_folder_by_wf.keys()))
                            cmd_archivage.run(log)
                            if send_mail:
                                sendMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "[NGS] Sequencer data archivage", "The run {} has been archived.".format(in_run_folder))
                    except Exception:
                        log.error("The post-process on run {} has failed.".format(run_id))
                        if send_mail:
                            sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "fail", in_run_folder)
                        with open(failed_analyses_file, "w") as FH:
                            FH.write('{}: error in step "{}"\n\nStacktrace:\n{}'.format(
                                datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                step,
                                traceback.format_exc()
                            ))
                    log.info("End post-process on run {}.".format(run_id))
        time.sleep(args.roll_time)
