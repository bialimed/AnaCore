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
__version__ = '2.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.annotVcf import AnnotVCFIO, getAlleleRecord, VCFRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def isSameAlt(var_record, var_annot):
    """
    Evaluate an annotation of a variant and return False if it concerns an collocated variant and not the same variant.

    :param var_record: The variant. It must contain only one alternative allele.
    :type var_record: anacore.vcf.VCFRecord
    :param var_annot: One annotation of the variant (for example coming from ANN INFO field produced by VEP).
    :type var_annot: dict
    :return: True if the var_record and var_annot concern the same alternative allele.
    :rtype: boolean
    """
    record_alt_str = var_record.alt[0].replace(".", "").replace("-", "").upper()
    annot_alt_str = var_annot["Allele"].replace(".", "").replace("-", "").upper()
    is_self_variant = (record_alt_str == annot_alt_str)
    if not is_self_variant:
        std_record = VCFRecord(
            var_record.chrom,
            var_record.pos,
            None,
            var_record.ref,
            var_record.alt
        )
        std_record.standardizeSingleAllele()
        record_alt_str = std_record.alt[0].replace(".", "").replace("-", "").upper()
        is_self_variant = (record_alt_str == annot_alt_str)
    return is_self_variant


def getAnnotSummary(allele_record, initial_alt, annot_field="ANN", pop_prefixes=None, pathogenicity_fields=None, logger=None):
    """
    Return a summary of the diffrent annotations of the variant. This summary is about identical known variants (xref), AF in populations (pop_AF), annotations of the variant and annotations of the collocated variants.

    :param allele_record: The variant.
    :type allele_record: anacore.vcf.VCFRecord
    :param initial_alt: The alternative variant before any transformation (standardization, ...). It must be the alternative variant directly extract from the annotated VCF.
    :type initial_alt: str
    :param annot_field: Field used to store annotations.
    :type annot_field: str
    :param pop_prefixes: The prefixes used to determine database name in population allele frequency fields. Example: the prefix "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF.
    :type pop_prefixes: list
    :param pathogenicity_fields: The titles of fields used to store pathogenicity predictor results. Example: SIFT, PolyPhen, CADD_PHRED.
    :type pathogenicity_fields: str
    :param logger: The logger object.
    :type loggger: logging.Logger
    :return: First the dentical known variants (e.g. {"cosmic": ["COSM14", "COSM15"], "dbSNP":[]}), second AF in populations (e.g. [{"source":"1KG", "name":"Global", "AF":0.85}]), third annotations of the variant and fourth annotations of the collocated variants.
    :rtype: list
    :warnings: The allele_record must only contains one variant.
    """
    xref = {"cosmic": set(), "dbSNP": set(), "HGMD": set(), "Unknown": set()}
    pop_AF = dict()
    variant_annot = list()
    collocated_annot = list()
    for annot in allele_record.info[annot_field]:
        is_self_variant = isSameAlt(allele_record, annot)
        # Similar knowns variants
        if is_self_variant and annot["Existing_variation"] is not None:
            for db_id in annot["Existing_variation"].split("&"):
                if db_id.startswith("rs"):
                    xref["dbSNP"].add(db_id)
                elif db_id.startswith("COSM"):
                    xref["cosmic"].add(db_id)
                elif db_id.startswith("CM"):
                    xref["HGMD"].add(db_id)
                else:
                    xref["Unknown"].add(db_id)
                    if logger is not None:
                        logger.warning('The database using the variant ID "{}" is not managed by "{}".'.format(db_id, sys.argv[0]))
        # Allele frequency in populations
        if is_self_variant:
            for key in annot:
                if key.endswith("_AF") and annot[key] is not None:
                    source, name = getPopInfo(key, pop_prefixes, logger)
                    pop_id = source + "_" + name
                    for curr_AF in annot[key].split("&"):
                        if pop_id not in pop_AF:
                            pop_AF[pop_id] = {
                                "source": source,
                                "name": name,
                                "AF": float(curr_AF)
                            }
                        else:
                            if pop_AF[pop_id]["AF"] != float(curr_AF):
                                raise Exception(
                                    'The allele frequency for the variant {} in population {} is reported several times with different values in "{}".'.format(
                                        allele_record.getName(), pop_id, args.input_variants
                                    )
                                )
        # Annotations
        annot_container = variant_annot if is_self_variant else collocated_annot
        annot_container.append({
            "subject": {"symbol": annot["SYMBOL"], "feature": annot["Feature"], "feature_type": annot["Feature_type"]},
            "changes": {"HGVSc": annot["HGVSc"], "HGVSp": annot["HGVSp"]},
            "conseq": annot["Consequence"],
            "pathogenicity": getPathogenicityPredictors(annot, pathogenicity_fields)
        })
    xref = {db: list(xref[db]) for db in xref}
    pop_AF = list(pop_AF.values())
    return xref, pop_AF, variant_annot, collocated_annot


def getPathogenicityPredictors(annot, pathogenicity_fields=None):
    """
    Return by predictor the predicted pathogenicity.

    :param annot: The information from an annotation feature.
    :type annot: dict
    :param pathogenicity_fields: The titles of fields used to store pathogenicity predictor results. Example: SIFT, PolyPhen, CADD_PHRED.
    :type pathogenicity_fields: str
    :return: by predictor the predicted pathogenicity.
    :rtype: dict
    """
    pathogenicity_fields = ["CLIN_SIG", "CADD_PHRED", "MetaLR_rankscore", "VEST3_rankscore"] if pathogenicity_fields is None else pathogenicity_fields
    rename = {
        "CLIN_SIG": "ClinVar",
        "CADD_PHRED": "CADD_phred",
        "MetaLR_rankscore": "MetaLR",
        "VEST3_rankscore": "VEST3"
    }
    score_by_predictor = {}
    for key in pathogenicity_fields:
        if key in annot and annot[key] != "":
            source = key
            if key in rename:
                source = rename[key]
            score_by_predictor[source] = annot[key]
    return score_by_predictor


def getPopInfo(annot_key, pop_prefixes=None, logger=None):
    """
    Return source and name of a sub-population studied in large genomic programs (1KG, ExAC, ...) from the annotation tag.

    :param annot_key: The tag used in annotation to store AF of the variant in sub-population.
    :type annot_key: str
    :param pop_prefixes: The prefixes used to determine database name in population allele frequency fields (example: "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF).
    :type pop_prefixes: list
    :param logger: The logger object.
    :type loggger: logging.Logger
    :return: The source (project name) and the name of the sub-population.
    :rtype: list
    """
    pop_lower_prefixes = ["exac", "gnomad", "1kg", "esp"] if pop_prefixes is None else [elt.lower() for elt in pop_prefixes]
    source = None
    name = None
    if annot_key == "AF":
        source = "1KG"
        name = "Global"
    elif "_" in annot_key and annot_key.lower().split("_")[0] in pop_lower_prefixes:  # ExAC_AF, ExAC_Adj_AF, ExAC_AFR_AF, ExAC_AMR_AF, ...
        source = annot_key.split("_")[0]
        name = "Global" if annot_key.count('_') == 1 else annot_key.split("_")[1]
    elif annot_key.count('_') == 1:  # AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, ...
        source = "1KG"
        name = annot_key.split("_")[0]
    else:
        if logger is not None:
            logger.warning('The population information stored with tag "{}" cannot be used by "{}".'.format(annot_key, sys.argv[0]))
    return source, name


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Converts VCF annotated with VEP in JSON format.')
    parser.add_argument('-r', '--assembly-id', help='ID of the reference used for the variant calling (example: GRCh38.p12).')
    parser.add_argument('-t', '--pathogenicity-fields', default=["CLIN_SIG", "CADD_PHRED", "MetaLR_rankscore", "VEST3_rankscore"], nargs='+', help='The titles of fields used to store pathogenicity predictor results (example: SIFT, PolyPhen, CADD_PHRED). [Default: %(default)s]')
    parser.add_argument('-p', '--populations-prefixes', default=["exac", "gnomad", "1kg", "esp"], nargs='+', help='The prefixes used to determine database name in population allele frequency fields (example: "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF). [Default: %(default)s]')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-source', default=None, help='Add source of the calling in support information.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file file containing variants and annotated with VEP v88+ (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='The path to the file outputted file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Convert VCF in python dict
    json_data = list()
    with AnnotVCFIO(args.input_variants, "r", args.annotation_field) as FH_vcf:
        for record in FH_vcf:
            for idx_alt in range(len(record.alt)):
                allele_record = getAlleleRecord(FH_vcf, record, idx_alt)
                allele_record.standardizeSingleAllele()
                curr_json = dict()
                # Coord information
                curr_json["coord"] = {
                    "region": allele_record.chrom,
                    "pos": allele_record.pos,
                    "ref": allele_record.ref,
                    "alt": allele_record.alt[0],
                    "assembly": (None if args.assembly_id is None else args.assembly_id)
                }
                # Support information
                curr_json["supports"] = [{
                    "filters": allele_record.filter,
                    "quality": allele_record.qual,
                    "libraries": [{"alt_depth": allele_record.getAD(library)[0], "depth": allele_record.getDP(library), "name": library} for library in FH_vcf.samples],
                    "source": args.calling_source
                }]
                if FH_vcf.annot_field in allele_record.info:
                    # Identical known variants, AF in populations and annotations
                    curr_json["xref"], curr_json["pop_AF"], curr_json["annot"], curr_json["collocated_annot"] = getAnnotSummary(
                        allele_record,
                        record.alt[idx_alt],
                        args.annotation_field,
                        args.populations_prefixes,
                        args.pathogenicity_fields,
                        log
                    )
                json_data.append(curr_json)

    # Write output file
    with open(args.output_variants, "w") as FH_out:
        FH_out.write(
            json.dumps(json_data, default=lambda o: o.__dict__, sort_keys=True)
        )
    log.info("END of process.")
