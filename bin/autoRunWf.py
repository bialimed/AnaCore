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
__version__ = '2.11.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
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
    """
    @summary: Class to describe design and his resources.
    """
    def __init__(self, name, pos_ctrl_re=None, neg_ctrl_re=None, resources_folder=None):
        self.name = name
        self.pos_ctrl_re = pos_ctrl_re
        self.neg_ctrl_re = neg_ctrl_re
        self.resources_folder = resources_folder

    def getAssemblyFolder(self, assembly):
        """
        @summary: Returns the resource folder for the design in the selected assembly.
        @param assembly: [str] The selected assembly.
        @return: [str] Path to the resource folder.
        """
        return os.path.join(self.resources_folder, assembly)

    def getPosCtrlRef(self, assembly):
        """
        @summary: Returns the file containing the reference values for the positive control (for example the VCF containing the expected variants).
        @param assembly: [str] The selected assembly.
        @return: [str] Path to the reference values file.
        """
        return os.path.join(self.getAssemblyFolder(assembly), "pos_ctrl_expected.vcf")

    def getVarNoise(self, assembly):
        """
        @summary: Returns the file containing the list of constitutive variants with their maximum AF.
        @param assembly: [str] The selected assembly.
        @return: [str] Path to the reference values file.
        """
        return os.path.join(self.getAssemblyFolder(assembly), "variants_noise.tsv")

    def getSelectedRef(self, ref_type):
        """
        @summary: Returns the file containing the reference selected (for example the list of RNA where the variants are kept).
        @param ref_type: [str] The type of reference.
        @return: [str] Path to the file.
        """
        return os.path.join(self.resources_folder, "reference_" + ref_type + ".tsv")

    def getFilters(self, wf_name):
        """
        @summary: Returns the file containing the filters to apply on the workflow results.
        @param wf_name: [str] The workflow name.
        @return: [str] Path to the filters description file.
        """
        return os.path.join(self.resources_folder, wf_name + "_filters.json")


def exec_cmd(cmd, log, shell=False):
    """
    @summary: Log and submit the command.
    @param cmd: [list] The command to execute.
    @param log: [Logger] The logger of the script.
    @param shell: [boolean] If shell is True, the specified command will be executed through the shell.
    """
    log.info("Start sub-command: " + " ".join(cmd))
    if shell:
        subprocess.check_call(" ".join(cmd), stdout=subprocess.DEVNULL, shell=shell)
    else:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
    log.info("End sub-command")


def getProtocol(in_spl_folder):
    """
    @summary: Returns the analysis protocol (design and Illumina workflow) declarated in samplesheet.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @return: [dict] The analysis protocol.
    @warning: For Amplicon - DS the manifests names must be <DESIGN>_A.txt and <DESIGN>_B.txt.
    """
    protocol = {"workflow": None, "design": None}
    samplesheet = SampleSheetIO(os.path.join(in_spl_folder, "SampleSheet.csv"))
    protocol["workflow"] = samplesheet.header["Workflow"]
    if samplesheet.header["Application"] == "Amplicon - DS" or samplesheet.header["Workflow"] == "Amplicon - DS":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        design_name = manifest_A.split("_A.txt")[0]
        pos_ctrl_re = ".*[Tt][Ee][Mm][Oo][Ii][Nn].*"
        if design_name.startswith("INCa") or design_name.startswith("TSCA_INCa"):
            pos_ctrl_re = ".*[Hh][Oo][Rr][Ii].*"
        protocol["design"] = Design(design_name, pos_ctrl_re, ".*[Nn][Tt][Cc].*")
    elif samplesheet.header["Workflow"] == "Enrichment":
        manifest_A = os.path.basename(samplesheet.manifests["A"])
        protocol["design"] = Design(
            manifest_A.split(".txt")[0], None, ".*[Nn][Tt][Cc].*"
        )
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


def getRunCmd(workflow, in_spl_folder, out_run_folder, reference):
    """
    @summary: Returns the command to launch the workflow.
    @param workflow: [str] The name of the workflow.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param out_run_folder: [str] Path to the workflow output folder.
    @param reference: [dict] Reference description: assembly and sequences.
    @return: [None/list] The command if the corresponding workflow exists and None otherwise.
    """
    cmd = None
    if workflow == "ADIVaR":
        cmd = getADIVaRCmd(in_spl_folder, out_run_folder, protocol["design"], reference)
    elif workflow == "AmpliconDS":
        cmd = getADSACmd(in_spl_folder, out_run_folder, protocol["design"], reference)
    return cmd


def getLoadCmd(design, raw_folder, out_by_wf):
    """
    @summary: Returns the command to load run and workflows data in laboratory database.
    @param design: [Design] The object describing the design.
    @param raw_folder: [str] Path to the run folder.
    @param out_by_wf: [dict] The workflow output path by workflow name.
    @return: [list] The command to load data.
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
    return cmd


def getADSACmd(in_spl_folder, out_run_folder, design, reference):
    """
    @summary: Returns the command to launch the variant annotation on MiSeq reporter outputs.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param design: [Design] The object describing the design.
    @param reference: [dict] Reference description: assembly and sequences.
    @return: [list] The command.
    """
    cmd = [
        args.workflow_path_pattern.replace("WF_NAME", "AmpliconDSAnnot"), "amplicondsannot",
        "--RNA-selection", design.getSelectedRef("RNA").replace("WF_NAME", "AmpliconDSAnnot"),
        "--assembly-version", reference["assembly"],
        "--filters", design.getFilters("ADSA").replace("WF_NAME", "AmpliconDSAnnot"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = design.getPosCtrlRef(reference["assembly"] + "_chr").replace("WF_NAME", "AmpliconDSAnnot")
    if os.path.exists(positive_ctrl_ref):
        cmd.extend([
            "--pos-ctrl-names", design.pos_ctrl_re,
            "--pos-ctrl-expected", positive_ctrl_ref
        ])
    return cmd


def getADIVaRCmd(in_spl_folder, out_run_folder, design, reference):
    """
    @summary: Returns the command to launch the amplicon double strand workflow.
    @param in_spl_folder: [str] Path to the sample folder (containing fastq and the samplesheet).
    @param out_run_folder: [str] Path to the workflow output folder.
    @param design: [Design] The object describing the design.
    @param reference: [dict] Reference description: assembly and sequences.
    @return: [list] The command.
    """
    resources_folder = os.path.dirname(args.design_path_pattern.replace("WF_NAME", "ADIVaR"))
    assembly_folder = design.getAssemblyFolder(reference["assembly"]).replace("WF_NAME", "ADIVaR")
    cmd = [
        args.workflow_path_pattern.replace("WF_NAME", "ADIVaR"), "adivar",
        "--R1-end-adapter", os.path.join(resources_folder, "adapters", "Illumina_3prim_adapter.fasta"),
        "--R2-end-adapter", os.path.join(resources_folder, "adapters", "Illumina_5prim_adapter_rvc.fasta"),
        "--libA", "folderA=" + os.path.join(assembly_folder, "libA"),
        "--libB", "folderB=" + os.path.join(assembly_folder, "libB"),
        "--RNA-selection", design.getSelectedRef("RNA").replace("WF_NAME", "ADIVaR"),
        "--assembly-version", reference["assembly"],
        "--genome-seq", reference["sequences"],
        "--filters", design.getFilters("ADIVaR").replace("WF_NAME", "ADIVaR"),
        "--samplesheet", os.path.join(in_spl_folder, "SampleSheet.csv"),
        "--output-dir", out_run_folder
    ]
    positive_ctrl_ref = design.getPosCtrlRef(reference["assembly"]).replace("WF_NAME", "ADIVaR")
    if os.path.exists(positive_ctrl_ref):
        cmd.extend([
            "--pos-ctrl-names", design.pos_ctrl_re,
            "--pos-ctrl-expected", positive_ctrl_ref
        ])
    constit_variants = design.getVarNoise(reference["assembly"]).replace("WF_NAME", "ADIVaR")
    if os.path.exists(constit_variants):
        cmd.extend([
            "--constit-variants", constit_variants
        ])
    return cmd


def getArchiveCmd(profiles_folder, run, workflows):
    """
    @summary: Returns the command to archivate run and workflows data.
    @param profiles_folder: [str] Path to the folder containing profiles describing the data to archivate.
    @param run: [str] Name of the run folder.
    @param workflows: [list] The names of the workflows to archivate.
    @return: [list] The command to archivate data.
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
    return cmd


def sendMail(smtp_adress, sender, recipients, subject, content):
    """
    @summary: Send execution log mail.
    @param smtp_adress: [str] The SMTP server used to send mail.
    @param sender: [str] Mail of the sender (it may be a not  existing adrsesses).
    @param recipients: [list] Mails of the recipients.
    @param subject: [str] The mail subject.
    @param content: [str] The mail content.
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
    @summary: Send execution log mail.
    @param smtp_adress: [str] The SMTP server used to send mail.
    @param sender: [str] Mail of the sender (it may be a not  existing adrsesses).
    @param recipients: [list] Mails of the recipients.
    @param status: [str] The status of the message: start or end or fail.
    @param run_folder: [str] Path to the run folder currently processed.
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
    group_output.add_argument('-a', '--analysis-folder', required=True, help='Path to the workflows output folder. Each workflow will create a sub-folder in this parent folder (example: for the analysis foler "170101_run01" the result can be two sub-folders 170101_run01/ADIVaR and 170101_run01/AmpliconDS).')
    group_output = parser.add_argument_group('Resources')  # Resources
    group_output.add_argument('-w', '--workflow-path-pattern', required=True, help='The path to the jflow client of the workflows. The tag WF_NAME is used as placeholder to adapt the command at the workflow. Example: /usr/local/bioinfo/WF_NAME/current/bin/jflow_cli.py.')
    group_output.add_argument('-d', '--design-path-pattern', required=True, help='The path to the folder containing the designs definitions folders. The tag WF_NAME is used as placeholder to adapt the command at the workflow. Example: /usr/local/bioinfo/WF_NAME/current/resources.')
    group_output.add_argument('-g', '--genome', required=True, help='The assembly version and the path to the reference genome used. The two information are separated by ":" like in the following example: GRCh37:/bank/Homo_sapiens/GRCH37/dna.fasta.')
    group_mail = parser.add_argument_group('Mail')  # Mail
    group_mail.add_argument('--mail-recipients', nargs='+', help='Mails of the recipients.')
    group_mail.add_argument('--mail-sender', help='Mail of the sender (it may be a not  existing adrsesses).')
    group_mail.add_argument('--mail-smtp', help='The SMTP server used to send mail.')
    args = parser.parse_args()
    genome = {
        'assembly': args.genome.split(":")[0],
        'sequences': args.genome.split(":")[1]
    }
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
                        in_basecalls_folder = os.path.join(in_run_folder, "Data", "Intensities", "BaseCalls")
                        protocol = getProtocol(in_basecalls_folder)
                        protocol["design"].resources_folder = os.path.join(args.design_path_pattern, protocol["design"].name)
                        out_folder_by_wf = dict()
                        if not os.path.exists(out_run_folder):
                            os.mkdir(out_run_folder)
                        # Launch workflows
                        step = "Launch analysis workflows"
                        for curr_wf in getWorkflows(protocol):
                            out_wf_folder = os.path.join(out_run_folder, curr_wf)
                            out_folder_by_wf[curr_wf] = out_wf_folder
                            if os.path.exists(out_wf_folder):
                                log.warning('The workflow "{}" has already be processed for run "{}".'.format(curr_wf, in_run_folder))
                            else:
                                os.mkdir(out_wf_folder)
                                cmd_analysis = getRunCmd(curr_wf, in_basecalls_folder, out_wf_folder, genome)
                                exec_cmd(cmd_analysis, log)
                        # Copy raw data
                        step = "Copy raw"
                        if args.storage_folder is not None:
                            out_run_storage = os.path.join(args.storage_folder, run_id)
                            if not os.path.exists(out_run_storage):  # The run folder has not been copied in storage folder
                                cmd_copy_raw = ["rsync", "--recursive", "--perms", "--times", in_run_folder + os.sep, out_run_storage]
                                exec_cmd(cmd_copy_raw, log)
                                # Remove copy of the Illumina's workflows
                                storage_basecalls_folder = os.path.join(out_run_storage, "Data", "Intensities", "BaseCalls")
                                alignments_files = glob.glob(os.path.join(storage_basecalls_folder, "Alignment*"))
                                for curr_file in alignments_files:
                                    if os.path.isdir(curr_file):
                                        cmd_rm_analysis = ["rm", "-r", curr_file]
                                        exec_cmd(cmd_rm_analysis, log)
                        # Copy Illumina's workflows
                        step = "Copy Illumina's workflows"
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
                                    if analysis_basename not in ["AmpliconDS"]: out_folder_by_wf[analysis_basename] = out_wf_folder
                                    log_file = os.path.join(out_wf_folder, "CompletedJobInfo.xml")
                                    if not os.path.exists(log_file):  # The results have not been copied in workflow folder
                                        cmd_copy_analysis = ["rsync", "--recursive", "--perms", "--times", in_wf_folder + os.sep, out_wf_folder]
                                        exec_cmd(cmd_copy_analysis, log)
                        # Save in database
                        step = "Database storage"
                        out_run_storage = in_run_folder
                        if args.storage_folder is not None:
                            out_run_storage = os.path.join(args.storage_folder, run_id)
                        cmd_load = getLoadCmd(protocol["design"], out_run_storage, out_folder_by_wf)
                        exec_cmd(cmd_load, log)
                        # Log end process
                        with open(completed_analyses_file, "w") as FH:
                            FH.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                        if send_mail:
                            sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "end", in_run_folder)
                        # Send to archive
                        step = "Send to archive"
                        if args.archivage_profiles_folder is not None:
                            cmd_archivage = getArchiveCmd(args.archivage_profiles_folder, run_id, list(out_folder_by_wf.keys()))
                            exec_cmd(cmd_archivage, log, True)
                            if send_mail:
                                sendMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "[NGS] Sequencer data archivage", "The run {} has been archived.".format(in_run_folder))
                    except Exception:
                        log.error("The post-process on run {} has failed.".format(run_id))
                        if send_mail:
                            sendStatusMail(args.mail_smtp, args.mail_sender, args.mail_recipients, "fail", in_run_folder)
                        with open(failed_analyses_file, "w") as FH:
                            FH.write('{}: error in step "{}"\n\nStacktrace:\n{}'.format(
                                datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
                                step,
                                traceback.print_exc()
                            )
                    log.info("End post-process on run {}.".format(run_id))
        time.sleep(args.roll_time)
