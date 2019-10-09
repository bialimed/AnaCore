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
__version__ = '0.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import re
import os
import sys
import yaml
import copy
import logging
import argparse
import hgvs.posedit
import hgvs.location
import hgvs.assemblymapper
import hgvs.sequencevariant
import hgvs.dataproviders.uta

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
def getAssemblyMapper(assembly_version, input_UTA_config=None):
    """
    Return UTA assemby mapper.

    :param assembly_version: Human genome assembly version.
    :type assembly_version: str
    :param input_UTA_config: Path to the configuration file for connexion to the UTA database (format: YAML). It must contain dialect, login, password, host, port and database. [Default: Connection to UTA public database]
    :type input_UTA_config: str
    :return: UTA assembly mapper
    :rtype: hgvs.assemblymapper.AssemblyMapper
    """
    uta_url = None
    if input_UTA_config:
        with open(input_UTA_config) as FH_in:
            uta_cfg = yaml.safe_load(FH_in)
            uta_url = "{}://{}:{}@{}:{}/uta/{}".format(
                uta_cfg.get("dialect"),
                uta_cfg.get("login"),
                uta_cfg.get("password"),
                uta_cfg.get("host"),
                uta_cfg.get("port"),
                uta_cfg.get("database")
            )
    hdp = hgvs.dataproviders.uta.connect(uta_url)
    assembly_mapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=assembly_version, alt_aln_method='splign', replace_reference=True)
    return assembly_mapper


def getHGVSGFromVCFRec(record, chr_accession=None, alt=None):
    """
    Return HGVS sequence variant (HGVSg) from a VCFRecord.

    :param record: The record coming from a VCF file.
    :type record: anacore.vcf.VCFRecord
    :param chr_accession: The chromosome RefSeq accession used in returned HGVS. [Default: use the chromosme name in VCFRecord]
    :type chr_accession: str
    :param alt: The alternative ariant used in returned HGVS. [Default: use the alternative allele in VCFRecord]
    :type alt: str
    :return: HGVS sequence variant (HGVSg).
    :rtype: hgvs.sequencevariant.SequenceVariant
    """
    chr_id = record.chrom if chr_accession is None else chr_accession
    alt = record.alt[0] if alt is None else alt
    change = hgvs.edit.NARefAlt(
        ref=(None if "-" == record.ref else record.ref),
        alt=(None if "-" == alt else alt)
    )
    start = record.pos
    end = start
    if len(record.ref) > 1:
        end += len(record.ref) - 1
    if record.ref == "-":
        start -= 1
        end = start + 1
    interval = hgvs.location.Interval(
        start=hgvs.location.SimplePosition(start),
        end=hgvs.location.SimplePosition(end)
    )
    return hgvs.sequencevariant.SequenceVariant(
        ac=chr_id,
        type="g",
        posedit=hgvs.posedit.PosEdit(pos=interval, edit=change)
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Fix or add HGVSg, HGVSc and HGVSp on variants annotations. The HGVS used are based on biocommons/hgvs.')
    parser.add_argument('-s', '--assembly-version', default="GRCh38", help='Human genome assembly version used in alignment, variants calling and variants annotation. [Default: %(default)s]')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_input.add_argument('-u', '--input-UTA-config', help='Path to the configuration file for connexion to the UTA database (format: YAML). It must contain dialect, login, password, host, port and database. [Default: Connection to UTA public database]')
    group_input.add_argument('-c', '--input-assembly-accessions', help='Path to the file describing link between chromosome name and RefSeq accession (format: TSV). The header must contain: sequence_id<tab>...<tab>RefSeq_accession where sequence_id is the name of the chromosome in CHROM column of the VCF. [Default: the chromosome name in VCF is the RefSeq accession]')
    group_input.add_argument('-r', '--input-sequence-repository', help='Path to the HGVS seqrepo. [Default: Connection to the public database or the use the repository specified in environment variable $HGVS_SEQREPO_DIR]')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get accession by chromosome ID
    acc_by_chrom = {}
    if args.input_assembly_accessions:
        with HashedSVIO(args.input_assembly_accessions, title_starter=None) as FH:
            for record in FH:
                acc_by_chrom[record["sequence_id"]] = record["RefSeq_accession"]

    # Connect to HGVS mapper
    if args.input_sequence_repository is not None:
        os.environ["HGVS_SEQREPO_DIR"] = args.input_sequence_repository
    hgvs_mapper = getAssemblyMapper(args.assembly_version, args.input_UTA_config)

    # Write
    nb_records = {"analysed": 0, "fixed_HGVSg": 0, "fixed_HGVSc": 0, "fixed_HGVSp": 0}
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
            for record in FH_in:
                nb_records["analysed"] += 1
                if len(record.alt) > 1:
                    raise Exception("The record {} is multi-allelic.".format(record.getName()))
                chr_acc = acc_by_chrom[record.chrom] if args.input_assembly_accessions else record.chrom
                hgvs_by_allele = {}
                std_record = None
                if record.isIndel():
                    std_record = copy.deepcopy(record)
                    std_record.normalizeSingleAllele()
                # Get HGVSg by annotation allele
                for annot in record.info[args.annotations_field]:
                    if annot["Allele"] not in hgvs_by_allele:
                        hgvs_g = None
                        if std_record is None:  # is substit
                            hgvs_g = getHGVSGFromVCFRec(record, chr_acc, annot["Allele"])
                        else:  # is inDel
                            hgvs_g = getHGVSGFromVCFRec(std_record, chr_acc, annot["Allele"])
                        hgvs_by_allele[annot["Allele"]] = {
                            "HGVSg": hgvs_g,
                            "transcripts": {tr.split(".")[0]: tr for tr in hgvs_mapper.relevant_transcripts(hgvs_g)},
                        }
                # Update annotations
                is_fixed_HGVSg = False
                is_fixed_HGVSc = False
                is_fixed_HGVSp = False
                for annot in record.info[args.annotations_field]:
                    # Trace
                    old = {
                        "g": "" if "HGVSg" not in annot or annot["HGVSg"] is None else annot["HGVSg"],
                        "c": "" if "HGVSc" not in annot or annot["HGVSc"] is None else annot["HGVSc"],
                        "p": "" if "HGVSp" not in annot or annot["HGVSp"] is None else annot["HGVSp"]
                    }
                    # HGVSg
                    hgvs_g = hgvs_by_allele[annot["Allele"]]["HGVSg"]
                    str_hgvs_g = str(hgvs_g)
                    if not(old["g"] == "" and str_hgvs_g.endswith(":g.?")):
                        if str_hgvs_g.endswith("="):
                            match = re.search("\d([ATGCN]+>[ATGCN]+)$", old["c"])
                            str_hgvs_g = str_hgvs_g.replace("=", match.group(1))
                        annot["HGVSg"] = str_hgvs_g
                        if old["g"] != str_hgvs_g:
                            is_fixed_HGVSg = True
                    # HGVSc or HGVSp
                    transcripts = hgvs_by_allele[annot["Allele"]]["transcripts"]
                    curr_tr_id = annot["Feature"].split(".")[0]
                    if curr_tr_id in transcripts:
                        annot["HGVSg"] = str(hgvs_g)
                        annot["HGVSc"] = ""
                        annot["HGVSp"] = ""
                        try:
                            hgvs_c = hgvs_mapper.g_to_c(hgvs_g, transcripts[curr_tr_id])
                            str_hgvs_c = str(hgvs_c)
                            if not(old["c"] == "" and str_hgvs_c.endswith(":c.?")):
                                if str_hgvs_c.endswith("="):
                                    match = re.search("\d([ATGCN]+>[ATGCN]+)$", old["c"])
                                    str_hgvs_c = str_hgvs_c.replace("=", match.group(1))
                                annot["HGVSc"] = str_hgvs_c
                                if old["c"] != str_hgvs_c:
                                    is_fixed_HGVSc = True
                            hgvs_p = hgvs_mapper.c_to_p(hgvs_c)
                            str_hgvs_p = str(hgvs_p)
                            if not(old["p"] == "" and str_hgvs_p.endswith(":p.?")):
                                annot["HGVSp"] = str_hgvs_p
                                if old["p"] != str_hgvs_p:
                                    is_fixed_HGVSp = True
                        except hgvs.exceptions.HGVSUsageError:
                            pass
                # Trace results
                if is_fixed_HGVSg:
                    nb_records["fixed_HGVSg"] += 1
                if is_fixed_HGVSc:
                    nb_records["fixed_HGVSc"] += 1
                if is_fixed_HGVSp:
                    nb_records["fixed_HGVSp"] += 1
                FH_out.write(record)
    # Log
    log.info("{}/{} variants have been fixed on HGVSg.".format(nb_records["fixed_HGVSg"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSc.".format(nb_records["fixed_HGVSc"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSp.".format(nb_records["fixed_HGVSp"], nb_records["analysed"]))
    log.info("End of job")
