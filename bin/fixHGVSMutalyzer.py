#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import os
import sys
import logging
import argparse
import requests
import urllib.parse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.annotVcf import AnnotVCFIO
from anacore.sv import HashedSVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class RunMutalyzerLegend:
    """Class to manipulate field legend from runMutalyzer[Light]."""

    def __init__(self, legend_data):
        """
        Build and return an instance of RunMutalyzerLegend.

        :param legend_data: Legend field from runMutalyzer[Light].
        :type legend_data: list
        :return: The new instance.
        :rtype: RunMutalyzerLegend
        """
        self.data = legend_data

    def getIdByName(self):
        """
        Return element accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.

        :return: element accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.
        :rtype: dict
        """
        id_by_name = {}
        for elt in self.data:
            if "id" in elt and "name" in elt:
                id_by_name[elt["name"]] = elt["id"]
        return id_by_name

    def getProtBytr(self):
        """
        Return protein accession by transcript base accessione.

        :return: Protein accession by transcript base accession.
        :rtype: dict
        """
        prot_by_tr = {}
        id_by_name = self.getIdByName()
        for elt_name, elt_id in id_by_name.items():
            if "_i" in elt_name:
                prot_base_acc = elt_id.split(".")[0]
                tr_base_acc = id_by_name[elt_name.replace("_i", "_v")].split(".")[0]
                prot_by_tr[tr_base_acc] = prot_base_acc
        return prot_by_tr


def getHGVSgFromRec(record, acc_by_chrom=None, annotations_field="ANN"):
    """
    Return HGVSg from annotated VCF record.

    :param record: Annotated VCF record.
    :type record: anacore.annotVcf.VCFRecord
    :param acc_by_chrom: Chromosome RefSeq accession by chromosome name.
    :type acc_by_chrom: dict
    :param annotations_field: Field used to store annotations.
    :type annotations_field: str
    :return: HGVSg.
    :rtype: str
    """
    hgvs_g = set()
    for annot in record.info[annotations_field]:
        if "HGVSg" in annot and annot["HGVSg"] is not None and annot["HGVSg"] != "":
            hgvs_g.add(annot["HGVSg"])
    if len(hgvs_g) > 1:
        raise Exception("The variant {} is describes with several HGVSg: {}.".format(record.getName(), sorted(hgvs_g)))
    elif len(hgvs_g) != 0:
        hgvs_g = list(hgvs_g)[0]
        if not hgvs_g.startswith("NC"):
            chr_acc = acc_by_chrom[record.chrom]
            hgvs_g = hgvs_g.replace(record.chrom + ":", chr_acc + ":")
    else:  # Empty HGVSg
        hgvs_g = None
    return hgvs_g


def parseDesc(desc_data, id_by_name):
    """
    Return HGVS by RefSeq accession from proteinDescriptions or transcriptDescriptions of runMutalyzer[Light].

    :param desc_data: proteinDescriptions or transcriptDescriptions from runMutalyzer[Light].
    :type desc_data: list
    :param id_by_name: Protein/Transcript RefSeq accession by name. Example: {"KIT_v001": "NM_000222.2", "KIT_i001": "NP_000213.1"}.
    :type id_by_name: dict
    :return: HGVS by RefSeq accession.
    :rtype: dict
    """
    HGVS_by_elt = {}
    for elt in desc_data:
        match_renamed = re.search("(.+)\((.+)\):(.)\.(.+)", elt)
        if match_renamed:
            ref_acc, elt_name, hgvs_type, change = match_renamed.groups()
            elt_acc = elt_name
            if elt_name in id_by_name:
                elt_acc = id_by_name[elt_name]
            elt_base_acc = elt_acc.split(".")[0]
            if change == "?":
                HGVS_by_elt[elt_base_acc] = ""
            else:
                HGVS_by_elt[elt_base_acc] = "{}:{}.{}".format(elt_acc, hgvs_type, change)
        else:
            match_std = re.search("(.+):(.)\.(.+)", elt)
            if match_std:
                elt_acc, hgvs_type, change = match_std.groups()
                elt_base_acc = elt_acc.split(".")[0]
                if change == "?":
                    HGVS_by_elt[elt_base_acc] = ""
                else:
                    HGVS_by_elt[elt_base_acc] = "{}:{}.{}".format(elt_acc, hgvs_type, change)
            else:
                raise Exception("The feature description {} cannot be parsed.".format(elt))
    return HGVS_by_elt


def getHGVSByTr(res_data):
    """
    Return HGVSg, HGVSc/n and HGVSp by transcript base RefSeq accession from runMutalyzer[Light].

    :param res_data: Results from runMutalyzer[Light] (required fields: legend, genomicDescription, transcriptDescriptions and proteinDescriptions).
    :type res_data: list
    :return: HGVSg, HGVSc/n and HGVSp by transcript base RefSeq accession.
    :rtype: dict
    """
    legend = RunMutalyzerLegend(res_data["legend"])
    id_by_name = legend.getIdByName()
    prot_by_tr = legend.getProtBytr()
    new_HGVSg = res_data["genomicDescription"]
    HGVSc_by_tr = parseDesc(res_data["transcriptDescriptions"], id_by_name)
    HGVSp_by_prot = parseDesc(res_data["proteinDescriptions"], id_by_name)
    HGVS_by_tr = {}
    for tr_base_ac, HGVSc in HGVSc_by_tr.items():
        HGVSp = ""
        if tr_base_ac not in prot_by_tr:
            if tr_base_ac[1] != "R":  # Correspond to mRNA and the link with prot does not exist
                raise Exception("The protein ID for the transcript {} cannot be found in mutalyzer data: {}".format(tr_base_ac, res_data))
        else:
            prot_base_acc = prot_by_tr[tr_base_ac]
            HGVSp = HGVSp_by_prot[prot_base_acc]
        HGVS_by_tr[tr_base_ac] = {
            "HGVSg": new_HGVSg,
            "HGVSc": HGVSc,
            "HGVSp": HGVSp
        }
    return HGVS_by_tr


class LoggerAction(argparse.Action):
    """Manages logger level parameters (The value "INFO" becomes logging.info and so on)."""

    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Fix or add HGVSg, HGVSc and HGVSp on variants annotations. The HGVS used are based on biocommons/hgvs.')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-m', '--mutalyzer-url', default="https://mutalyzer.nl", help='URL to the mutalizer server. [Default: %(default)s]')
    parser.add_argument('-p', '--proxy-url', help='URL to the proxy server if the http(s) connexions are only allowed through a proxy.')
    parser.add_argument('-s', '--assembly-version', default="GRCh38", help='Human genome assembly version used in alignment, variants calling and variants annotation. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_input.add_argument('-c', '--input-assembly-accessions', help='Path to the file describing link between chromosome name and RefSeq accession (format: TSV). The header must contain: sequence_id<tab>...<tab>RefSeq_accession where sequence_id is the name of the chromosome in CHROM column of the VCF. [Default: the chromosome name in VCF is the RefSeq accession]')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Get accession by chromosome ID
    acc_by_chrom = {}
    if args.input_assembly_accessions:
        with HashedSVIO(args.input_assembly_accessions, title_starter=None) as FH:
            for record in FH:
                acc_by_chrom[record["sequence_id"]] = record["RefSeq_accession"]

    # Write
    nb_records = {"analysed": 0, "fixed_HGVSg": 0, "fixed_HGVSc": 0, "fixed_HGVSp": 0, "contains_colloc_annot": 0}
    with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotations_field) as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotations_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            if "HGVSg" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSg")
            if "HGVSc" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSc")
            if "HGVSp" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSp")
            FH_out.writeHeader()
            # Records
            # ~ print("HGVS_type", "Variant", "HGVS_g", "Feature", "Old_HGVS", "New_HGVS", "Fixed", sep="\t")
            param_assembly = urllib.parse.quote(args.assembly_version, safe='')
            for record in FH_in:
                nb_records["analysed"] += 1
                if len(record.alt) > 1:
                    raise Exception("The record {} is multi-allelic.".format(record.getName()))
                # Get info from mutalyzer
                old_HGVSg = getHGVSgFromRec(record, acc_by_chrom, FH_in.annot_field)
                if old_HGVSg is None:
                    log.warning("The variant {} does not contain any HGVSg.".format(record.getName()))
                else:
                    param_hgvsg = urllib.parse.quote(old_HGVSg, safe='')
                    param_fields = urllib.parse.quote(",".join(["legend", "proteinDescriptions", "transcriptDescriptions", "genomicDescription"]), safe='')
                    url_request = '{}/json/runMutalyzerLight?build={};variant={};extra={}'.format(args.mutalyzer_url, param_assembly, param_hgvsg, param_fields)
                    log.debug(url_request)
                    response = requests.get(
                        url_request,
                        proxies=(None if args.proxy_url is None else {"https": args.proxy_url, "http": args.proxy_url})
                    )
                    if response.status_code != 200:
                        raise Exception("Request {} has failed.".format(url_request))
                    res_data = response.json()
                    # Store all HGVS by transcript base accession
                    mutalyzer_tr = {elt["id"].split(".")[0] for elt in res_data["legend"] if "id" in elt}
                    annot_tr = {annot["Feature"].split(".")[0] for annot in record.info[args.annotations_field]}
                    if len(annot_tr - mutalyzer_tr) != 0:
                        log.warning("All the transcripts annotated for variant {} cannot be found in used version of mutalyzer. Missing transcripts: {}".format(record.getName(), sorted(annot_tr - mutalyzer_tr)))
                    HGVS_by_tr = getHGVSByTr(res_data)
                    # Update annotations
                    is_fixed_HGVSg = False
                    is_fixed_HGVSc = False
                    is_fixed_HGVSp = False
                    contains_colloc_annot = False
                    new_annot = []
                    for annot in record.info[args.annotations_field]:
                        if annot["Allele"] != record.alt[0]:  # Annotation come from a collocated alternative allele
                            contains_colloc_annot = True
                        else:  # Annotation come from the alternative allele
                            tr_base_acc = annot["Feature"].split(".")[0]
                            if tr_base_acc in HGVS_by_tr:
                                # Trace
                                old = {
                                    "g": "" if "HGVSg" not in annot or annot["HGVSg"] is None else annot["HGVSg"],
                                    "c": "" if "HGVSc" not in annot or annot["HGVSc"] is None else annot["HGVSc"],
                                    "p": "" if "HGVSp" not in annot or annot["HGVSp"] is None else annot["HGVSp"]
                                }
                                # HGVSg
                                annot["HGVSg"] = HGVS_by_tr[tr_base_acc]["HGVSg"]
                                if old["g"] != annot["HGVSg"]:
                                    is_fixed_HGVSg = True
                                    # ~ print("HGVSg", record.getName(), annot["HGVSg"], annot["Feature"], old["g"], annot["HGVSg"], True, sep="\t")
                                # HGVSc
                                annot["HGVSc"] = HGVS_by_tr[tr_base_acc]["HGVSc"]
                                if old["c"] != annot["HGVSc"]:
                                    is_fixed_HGVSc = True
                                # ~ print("HGVSc", record.getName(), annot["HGVSg"], annot["Feature"], old["c"], annot["HGVSc"], old["c"] != annot["HGVSc"], sep="\t")
                                # HGVSp
                                annot["HGVSp"] = HGVS_by_tr[tr_base_acc]["HGVSp"]
                                if old["p"] != annot["HGVSp"]:
                                    is_fixed_HGVSp = True
                                # ~ print("HGVSp", record.getName(), annot["HGVSg"], annot["Feature"], old["p"], annot["HGVSp"], old["p"] != annot["HGVSp"].replace("(", "").replace(")", ""), sep="\t")
                            new_annot.append(annot)
                    record.info[args.annotations_field] = new_annot
                    # Trace results
                    if is_fixed_HGVSg:
                        nb_records["fixed_HGVSg"] += 1
                    if is_fixed_HGVSc:
                        nb_records["fixed_HGVSc"] += 1
                    if is_fixed_HGVSp:
                        nb_records["fixed_HGVSp"] += 1
                    if contains_colloc_annot:
                        nb_records["contains_colloc_annot"] += 1
                FH_out.write(record)
    # Log
    if nb_records["contains_colloc_annot"] != 0:
        log.warning("{}/{} variants contain collocated annotation removed during the process.".format(nb_records["contains_colloc_annot"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on HGVSg.".format(nb_records["fixed_HGVSg"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSc.".format(nb_records["fixed_HGVSc"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSp.".format(nb_records["fixed_HGVSp"], nb_records["analysed"]))
    log.info("End of job")
