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
__version__ = '2.15.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import yaml
import time
import glob
import hashlib
import smtplib
import logging
import argparse
import datetime
import traceback
import subprocess
from threading import Thread
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
class ThreadWrapper(Thread):
    """Thread wrapper to manage exit code like multiprocessig.Process."""

    def run(self):
        """
        Method representing the thread’s activity.

        You may override this method in a subclass. The standard run() method invokes the callable object passed to the object’s constructor as the target argument, if any, with sequential and keyword arguments taken from the args and kwargs arguments, respectively.
        """
        self.errcode = None
        try:
            Thread.run(self)
        except Exception:
            self.errcode = -1
            raise
        else:
            self.errcode = 0

    def exitcode(self):
        """
        Return the thread’s exit code.

        :return: This will be None if the process has not yet terminated. A negative value -N indicates that the child was terminated by signal N.
        :rtype: int
        """
        return self.errcode


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

    def run(self, log, cmd_id=None):
        """
        Log and submit the command.

        :param log: The logger of the script.
        :type log: Logger
        :param cmd_id: The ID used to identify command in log.
        :type cmd_id: str
        """
        cmd_str = " ".join([str(elt) for elt in self.cmd])
        if cmd_id is None:
            cmd_id = hashlib.md5(cmd_str.encode()).hexdigest()
        log.info("Start sub-command {}: {}".format(cmd_id, cmd_str))
        if self.is_shell:
            subprocess.check_call(" ".join([str(elt) for elt in self.cmd]), stdout=subprocess.DEVNULL, shell=True)
        else:
            subprocess.check_call(self.cmd, stdout=subprocess.DEVNULL)
        log.info("End sub-command {}".format(cmd_id))


def getIlluminaWorkflow(samplesheet_path):
    """
    Return the workflow declarated in the samplesheet.

    :param samplesheet_path: Path to the run samplesheet.
    :type samplesheet_path: str
    :return: The workflow declarated in the samplesheet.
    :rtype: str
    """
    workflow = None
    samplesheet = SampleSheetIO(samplesheet_path)
    if "Workflow" in samplesheet.header and samplesheet.header["Workflow"].strip() != "":
        workflow = samplesheet.header["Workflow"].strip()
    return workflow


def getDesigns(samplesheet_path):
    """
    Return the designs declarated in samplesheet.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :return: By design ID (<protocol>_<panel>_<version>) the design and the samples coming from this.
    :rtype: dict
    :warning: For Amplicon - DS the manifests names must be <DESIGN>_A.txt and <DESIGN>_B.txt.
    """
    info_by_design = {}
    samplesheet = SampleSheetIO(samplesheet_path)
    # Get design tag from samplesheet
    samples_are_tagged = False
    for curr_spl in samplesheet.samples:
        for curr_annot in curr_spl["Description"].split(","):
            if curr_annot.startswith("DESIGN#"):
                samples_are_tagged = True
    # By samples
    if samples_are_tagged:
        for curr_spl in samplesheet.samples:
            for curr_annot in curr_spl["Description"].split(","):
                if curr_annot.startswith("DESIGN#"):
                    design_id = curr_annot[7:]
                    if design_id in info_by_design:
                        info_by_design[design_id]["samples"].append(curr_spl)
                    else:
                        info_by_design[design_id] = {
                            "design": {
                                "id": design_id,
                                "protocol": design_id.split("_")[0],
                                "name": design_id.split("_")[1],
                                "version": design_id.split("_")[2]
                            },
                            "samples": [curr_spl]
                        }
    # By header
    if "DESIGN#" in samplesheet.header["Description"]:
        for curr_annot in samplesheet.header["Description"].split(","):
            if curr_annot.startswith("DESIGN#"):
                design_id = curr_annot[7:]
                info_by_design[design_id] = {
                    "design": {
                        "id": design_id,
                        "protocol": design_id.split("_")[0],
                        "name": design_id.split("_")[1],
                        "version": design_id.split("_")[2]
                    },
                    "samples": []
                }
                for curr_spl in samplesheet.samples:
                    if design_id in info_by_design:
                        info_by_design[design_id]["samples"].append(curr_spl)
    elif samplesheet.header["Application"] == "Amplicon - DS" or samplesheet.header["Workflow"] == "Amplicon - DS":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        design_id = manifest_A.split("_A.txt")[0]
        info_by_design[design_id] = {
            "design": {
                "id": design_id,
                "protocol": design_id.split("_")[0],
                "name": design_id.split("_")[1],
                "version": design_id.split("_")[2]
            },
            "samples": []
        }
        for curr_spl in samplesheet.samples:
            if design_id in info_by_design:
                info_by_design[design_id]["samples"].append(curr_spl)
    return info_by_design.values()


def getRenderedDefault(instance, default, replace_none=True):
    """
    Return configuration with missing values replaced by default values.

    :param instance: The configuration node to take as subject.
    :type instance: *
    :param default: The default configuration node. It must be at the same level as instance.
    :type default: *
    :param replace_none: If False the parameters declarated in config (key in config file) but without defined values does not take default values.
    :type replace_none: bool
    :return: The copy of instance with missing values replaced by default values.
    :rtype: *
    """
    new_val = None
    if not isinstance(default, dict):
        new_val = instance
    else:  # Value is dict
        new_val = {}
        for default_param, default_value in default.items():
            if default_param not in instance or (replace_none and instance[default_param] is None):
                new_val[default_param] = default_value
            elif isinstance(instance[default_param], dict):  # Value is dict
                new_val[default_param] = getRenderedDefault(default_value, instance[default_param])
        for param, value in instance.items():
            if param not in new_val:
                new_val[param] = value
    return new_val


def getRenderedPlaceholders(instance, placeholders):
    """
    Return configuration with placeholders replaced by values.

    :param instance: The configuration node to take as subject.
    :type instance: *
    :param placeholders: The values to use by placeholder.
    :type placeholders: dict
    :return: The copy of field_value with placeholders replaced by values.
    :rtype: *
    """
    new_val = None
    if isinstance(instance, dict):  # Value is dict
        new_val = {}
        for key, val in instance.items():
            new_val[key] = getRenderedPlaceholders(val, placeholders)
    elif isinstance(instance, (list, tuple)):  # Value is iterable
        new_val = []
        for val in instance:
            new_val.append(getRenderedPlaceholders(val, placeholders))
    elif isinstance(instance, set):  # Value is set
        new_val = set()
        for val in instance:
            new_val.add(getRenderedPlaceholders(val, placeholders))
    else:  # Value is number or string
        new_val = instance
        if isinstance(instance, str):  # Value is str
            for placeholder_name, placeholder_value in placeholders.items():
                new_val = instance.replace("$" + placeholder_name, placeholder_value)
    return new_val


def getWfConfig(cfg_path):
    """
    Return the workflows configuration.

    :param cfg_path: Path to the configuration file (format: YAML).
    :type cfg_path: str
    :return: workflows configuration by workflow name.
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
    # load parameters from file
    with open(cfg_path, 'r') as FH_cfg:
        cfg = yaml.load(FH_cfg)
    # Render value
    new_cfg = {}
    for wf_name, wf_param in cfg["workflows"].items():
        # Replace empty by default values
        new_cfg[wf_name] = getRenderedDefault(wf_param, cfg["default"])
        # Replace placeholders
        new_cfg[wf_name] = getRenderedPlaceholders(new_cfg[wf_name], {"WF_NAME": wf_name})
    return new_cfg


def getWorkflows(designs, cfg):
    """
    Return the workflows to execute on data.

    :param designs: The designs declarated for the run.
    :type designs: list
    :param cfg: Configuration for workflows.
    :type cfg: dict
    :return: By ID the workflows information.
    :rtype: dict
    """
    workflows = {}
    for curr_design in designs:
        for wf_name, wf_param in cfg.items():
            if wf_param["apply_on"]["protocol"] == "." or wf_param["apply_on"]["protocol"] == curr_design["design"]["protocol"]:
                if curr_design["design"] is None and wf_param["apply_on"]["designs"] is None:
                    wf_id = wf_name + "_" + curr_design["design"]["id"]
                    workflows[wf_id] = {
                        "id": wf_id,
                        "name": wf_name,
                        "design": curr_design["design"],
                        "samples": curr_design["samples"],
                        "resources": getRenderedPlaceholders(wf_param, {"DESIGN_ID": curr_design["design"]["id"]}),
                        "out_folder": None
                    }
                elif len(wf_param["apply_on"]["designs"]) == 0 or curr_design["design"]["id"] in wf_param["apply_on"]["designs"]:
                    wf_id = wf_name + "_" + curr_design["design"]["id"]
                    workflows[wf_id] = {
                        "id": wf_id,
                        "name": wf_name,
                        "design": curr_design["design"],
                        "samples": curr_design["samples"],
                        "resources": getRenderedPlaceholders(wf_param, {"DESIGN_ID": curr_design["design"]["id"]}),
                        "out_folder": None
                    }
    return workflows


def getRunCmd(wf_info, in_spl_folder, out_run_folder):
    """
    Return the command to launch the workflow.

    :param wf_info: The workflow information.
    :type wf_info: dict
    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type param out_run_folder: str
    :return: The command if the corresponding workflow exists and None otherwise.
    :rtype: Cmd
    """
    fct_name = "get{}Cmd".format(wf_info["name"])
    if fct_name not in globals():
        raise Exception("No {} has been implemented for the workflow {}.".format(fct_name, wf_info["name"]))
    cmd = globals()[fct_name](in_spl_folder, out_run_folder, wf_info)
    return cmd


def getLoadCmd(raw_folder, wf_by_id):
    """
    Return the command to load run and workflows data in laboratory database.

    :param raw_folder: Path to the run folder.
    :type raw_folder: str
    :param wf_by_id: The workflow information by design ID.
    :type wf_by_id: dict
    :return: The command to load data.
    :rtype: Cmd
    """
    cmd = [
        "lisImportRun.py",
        "--results-folder", "/Anapath/Illumina_Run_datas/Routine/database_files",  ################################
        "--input-raw", raw_folder
    ]
    # Control samples
    pos_ctrl = set()
    neg_ctrl = set()
    for wf_id, wf_info in wf_by_id.items():
        if "control" in wf_info["resources"]:
            if "neg_ctrl_regexp" in wf_info["resources"]["control"] and wf_info["resources"]["control"]["neg_ctrl_regexp"] is not None:
                neg_ctrl.add(wf_info["resources"]["control"]["neg_ctrl_regexp"])
            if "pos_ctrl_regexp" in wf_info["resources"]["control"] and wf_info["resources"]["control"]["pos_ctrl_regexp"] is not None:
                pos_ctrl.add(wf_info["resources"]["control"]["pos_ctrl_regexp"])
    if len(pos_ctrl) != 0:
        if len(pos_ctrl) > 1:
            raise Exception("The run cannot be load with different positive-ctrl-regexp: {}".format(pos_ctrl))
        cmd.extend(["--positive-ctrl-regexp", list(pos_ctrl)[0]])
    if len(neg_ctrl) != 0:
        if len(neg_ctrl) > 1:
            raise Exception("The run cannot be load with different negative-ctrl-regexp: {}".format(neg_ctrl))
        cmd.extend(["--negative-ctrl-regexp", list(neg_ctrl)[0]])
    # Workflows data
    if len(wf_by_id) != 0:
        cmd.append("--input-workflows")
        for wf_id, wf_info in wf_by_id.items():
            cmd.append("{}:{}".format(wf_info["name"], wf_info["out_folder"]))
    return Cmd(cmd)


def getADSACmd(in_spl_folder, out_run_folder, wf_param):
    """
    Return the command to launch the variant annotation on MiSeq reporter outputs.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param wf_param: Workflow's configuration.
    :type wf_param: dict
    :return: The command to run ADSA.
    :rtype: Cmd
    """
    wf_rsc = getRenderedPlaceholders(
        wf_param["resources"],
        {"ASSEMBLY": wf_param["resources"]["reference"]["assembly"]}
    )
    app_exec = os.path.join(wf_rsc["app_dir"], "app", "bin", "jflow_cli.py")
    cmd = [
        app_exec, "amplicondsannot",
        "--RNA-selection", wf_rsc["design"]["selected_RNA"],
        "--assembly-version", wf_rsc["reference"]["assembly"],
        "--filters", wf_rsc["design"]["filters"],
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    ctrl = wf_rsc["control"]
    if os.path.exists(ctrl["pos_ctrl_path"]):
        cmd.extend([
            "--pos-ctrl-names", ctrl["pos_ctrl_regexp"],
            "--pos-ctrl-expected", ctrl["pos_ctrl_path"]
        ])
    return Cmd(cmd)


def getADIVaRCmd(in_spl_folder, out_run_folder, wf_param):
    """
    Return the command to launch the amplicon double strand workflow.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param wf_param: Workflow's configuration.
    :type wf_param: dict
    :return: The command to run ADIVaR.
    :rtype: Cmd
    """
    wf_rsc = getRenderedPlaceholders(
        wf_param["resources"],
        {"ASSEMBLY": wf_param["resources"]["reference"]["assembly"]}
    )
    app_exec = os.path.join(wf_rsc["app_dir"], "app", "bin", "jflow_cli.py")
    cmd = [
        app_exec, "adivar",
        "--R1-end-adapter", wf_rsc["design"]["3prim_adapter"],
        "--R2-end-adapter", wf_rsc["design"]["5prim_adapter"],
        "--libA", "folderA=" + wf_rsc["design"]["libA_folder"],
        "--libB", "folderB=" + wf_rsc["design"]["libB_folder"],
        "--RNA-selection", wf_rsc["design"]["selected_RNA"],
        "--assembly-version", wf_rsc["reference"]["assembly"],
        "--genome-seq", wf_rsc["reference"]["sequences"],
        "--filters", wf_rsc["design"]["filters"],
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    ctrl = wf_rsc["control"]
    if os.path.exists(ctrl["pos_ctrl_path"]):
        cmd.extend([
            "--pos-ctrl-names", ctrl["pos_ctrl_regexp"],
            "--pos-ctrl-expected", ctrl["pos_ctrl_path"]
        ])
    constit_variants = wf_rsc["design"]["variants_noise"]
    if os.path.exists(constit_variants):
        cmd.extend([
            "--constit-variants", constit_variants
        ])
    return Cmd(cmd)


def getMIAmSCmd(in_spl_folder, out_run_folder, wf_param):
    """
    Return the command to launch microsatellite instability detection by next-generation sequencing on amplicons.

    :param in_spl_folder: Path to the sample folder (containing fastq and the samplesheet).
    :type in_spl_folder: str
    :param out_run_folder: Path to the workflow output folder.
    :type out_run_folder: str
    :param wf_param: Workflow's configuration.
    :type wf_param: dict
    :return: The command to run MIAmS.
    :rtype: Cmd
    """
    wf_rsc = getRenderedPlaceholders(
        wf_param["resources"],
        {"ASSEMBLY": wf_param["resources"]["reference"]["assembly"]}
    )
    app_env_bin = os.path.join(wf_rsc["app_dir"], "envs", "miniconda3", "bin")
    app_exec = os.path.join(wf_rsc["app_dir"], "jflow", "bin", "jflow_cli.py")
    cmd = [
        "unset", "PYTHONPATH", ";",
        "source", os.path.join(app_env_bin, "activate"), "MIAmS", "&&",
        app_exec, "miamstag",
        "--max-mismatch-ratio", 0.25,
        "--min-pair-overlap", 40,
        "--min-support-reads", 300,
        "--random-seed", 42,
        "--R1-end-adapter", wf_rsc["design"]["3prim_adapter"],
        "--R2-end-adapter", wf_rsc["design"]["5prim_adapter"],
        "--models", os.path.join(wf_rsc["design"]["MSI_folder"], "models.json"),
        "--targets", os.path.join(wf_rsc["design"]["MSI_folder"], "targets.bed"),
        "--intervals", os.path.join(wf_rsc["design"]["MSI_folder"], "intervals.tsv"),
        "--baseline", os.path.join(wf_rsc["design"]["MSI_folder"], "baseline.tsv"),
        "--genome-seq", wf_rsc["reference"]["sequences"],
        "--R1-pattern", "'" + os.path.join(in_spl_folder, "*R1_???.fastq.gz") + "'",
        "--R2-pattern", "'" + os.path.join(in_spl_folder, "*R2_???.fastq.gz") + "'",
        "--exclusion-pattern", "'" + "*Undetermined_S0_L???_R?_???.fastq.gz" + "'",
        "--output-dir", out_run_folder, "&&",
        "source", os.path.join(app_env_bin, "deactivate")
    ]
    return Cmd(cmd, True)


def archivate(profiles_folder, run, workflows, log):
    """
    Archivate run and workflows data.

    :param profiles_folder: Path to the folder containing profiles describing the data to archivate.
    :type profiles_folder: str
    :param run: Name of the run folder.
    :type run: str
    :param workflows: The workflows to archivate.
    :type workflows: list
    :param log: The logger instance.
    :type log: logging.Logger
    """
    # Raw data
    Cmd(
        [
            "archive.py",
            "--placeholders", "'@YEAR={}'".format("20" + run[0:2]), "'@RUN={}'".format(run),
            "--resources-profiles", os.path.join(profiles_folder, "Illumina", "*.yml")
        ],
        True
    ).run(log)
    # Workflows data
    for curr_wf in workflows:
        wf_profiles_folder = os.path.join(profiles_folder, "workflows", curr_wf["name"])
        if os.path.exists(wf_profiles_folder):
            Cmd(
                [
                    "archive.py",
                    "--placeholders", "'@YEAR={}'".format("20" + run[0:2]), "'@RUN={}'".format(run), "'@WF={}'".format(curr_wf["id"]),
                    "--resources-profiles", os.path.join(wf_profiles_folder, "*.yml")
                ],
                True
            ).run(log)
        else:
            log.warning("The data for the workflow {} of the run {} cannot be archivated.".format(curr_wf, run))


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
    parser.add_argument('-r', '--roll-time', type=int, default=(60 * 20), help="The time between each sequencer output folder evaluation (in seconds). [Default: %(default)s]")
    parser.add_argument('-p', '--archive-profiles-folder', help='Path to the folder containing profiles describing the data to archivate.')
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
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
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
                        samplesheet_path = os.path.join(in_basecalls_folder, "SampleSheet.csv")
                        out_folder_by_wf = dict()
                        # Launch workflows
                        step = "Launch analysis workflows"
                        parallel_processes = []
                        wf_by_id = getWorkflows(getDesigns(samplesheet_path), analyses_cfg)
                        for wf_id, wf_info in wf_by_id.items():
                            out_wf_folder = os.path.join(out_run_folder, wf_id)
                            wf_info["out_folder"] = out_wf_folder
                            if os.path.exists(out_wf_folder):
                                log.warning('The workflow "{}" has already be processed for run "{}".'.format(wf_id, in_run_folder))
                            else:
                                os.mkdir(out_wf_folder)
                                cmd_analysis = getRunCmd(wf_info, in_basecalls_folder, out_wf_folder)
                                wf_process = ThreadWrapper(name=wf_id, target=cmd_analysis.run, args=[log, wf_id])
                                parallel_processes.append(wf_process)
                                wf_process.start()
                        for wf_process in parallel_processes:  # Wait end of all workflows
                            wf_process.join()
                        wf_failed = [wf._name for wf in parallel_processes if wf.exitcode() != 0]  # Check workflows end status
                        if len(wf_failed) > 0:
                            raise Exception("Error in workflow(s): {}.".format(wf_failed))
                        out_folder_by_wf.pop("GenerateFASTQ", None)  # Remove unanalysed workflow
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
                        illumina_workflow = getIlluminaWorkflow(samplesheet_path)
                        if illumina_workflow is not None and illumina_workflow != "GenerateFASTQ":  # Illumina reporter has processed an analysis
                            # Get workflow name
                            analysis_basename = illumina_workflow.replace(" ", "_")
                            if analysis_basename == "Amplicon_-_DS":
                                analysis_basename = "ADSA"
                            else:
                                analysis_basename = "Illumina_" + analysis_basename
                            # Process workflow results
                            in_wf_folder_prefix = os.path.join(in_basecalls_folder, "Alignment")
                            alignments_files = glob.glob(in_wf_folder_prefix + "*")
                            for in_wf_folder in alignments_files:
                                if os.path.isdir(in_wf_folder):
                                    out_wf_folder = os.path.join(out_run_folder, analysis_basename)
                                    out_wf_folder += in_wf_folder.replace(in_wf_folder_prefix, "")  # Add workflow suffix (the index of the algnment folder)
                                    if analysis_basename not in ["ADSA"]:
                                        wf_by_id[analysis_basename] = {
                                            "id": analysis_basename,
                                            "name": analysis_basename,
                                            "out_folder": out_wf_folder,
                                        }
                                    log_file = os.path.join(out_wf_folder, "CompletedJobInfo.xml")
                                    if not os.path.exists(log_file):  # The results have not been copied in workflow folder
                                        cmd_copy_analysis = ["rsync", "--recursive", "--perms", "--times", in_wf_folder + os.sep, out_wf_folder]
                                        Cmd(cmd_copy_analysis).run(log)
                        # Save in database
                        step = "Database storage"
                        out_run_storage = in_run_folder
                        if args.storage_folder is not None:
                            out_run_storage = os.path.join(args.storage_folder, run_id)
                        cmd_load = getLoadCmd(out_run_storage, wf_by_id, log)
                        cmd_load.run(log)
                        # Log end process
                        with open(completed_analyses_file, "w") as FH:
                            FH.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                        if send_mail:
                            sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "end", in_run_folder)
                        # Send to archive
                        step = "Send to archive"
                        if args.archive_profiles_folder is not None:
                            archivate(args.archive_profiles_folder, run_id, list(wf_by_id.values()), log)
                            if send_mail:
                                sendMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "[NGS] Sequencer data archiving", "The run {} has been archived.".format(in_run_folder))
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
