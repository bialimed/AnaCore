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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import warnings
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.VEPvcf import VEPVCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAnnotSummary(allele_record, initial_alt):
    """
    @summary: Returns a summary of the diffrent annotations of the variant. This summary is about identical known variants (xref), AF in populations (pop_AF), annotations of the variant and annotations of the colocated variants.
    @param allele_record: [VCFRecord] The variant.
    @param initial_alt: [str] The alternative variant before any transformation (standardization, ...). It must be the alternative variant directly extract from the annotated VCF.
    @return: [list] Identical known variants (e.g. {"cosmic": ["COSM14", "COSM15"], "dbSNP":[]}), AF in populations (e.g. [{"source":"1KG", "name":"Global", "AF":0.85}]), annotations of the variant and annotations of the colocated variants.
    @warnings: The allele_record must only contains one variant.
    """
    xref = {"cosmic": set(), "dbSNP": set(), "HGMD": set(), "Unknown": set()}
    pop_AF = dict()
    variant_annot = list()
    colocated_annot = list()
    for annot in allele_record.info["CSQ"]:
        is_self_variant = (annot["Allele"] == initial_alt)
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
                    warnings.warn('The database using the variant ID "{}" is not managed by "{}".'.format(db_id, sys.argv[0]))
        # Allele frequency in populations
        if is_self_variant:
            for key in annot:
                if key.endswith("_AF") and annot[key] is not None:
                    source, name = getPopInfo(key)
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
        annot_container = variant_annot if is_self_variant else colocated_annot
        annot_container.append({
            "subject": {"symbol": annot["SYMBOL"], "feature": annot["Feature"], "feature_type": annot["Feature_type"]},
            "change": {"HGVSc": annot["HGVSc"], "HGVSp": annot["HGVSp"]},
            "conseq": {
                "tag": annot["Consequence"],
                "polyphen": annot["PolyPhen"],
                "sift": annot["SIFT"],
                "clinvar": annot["CLIN_SIG"]
            }
        })
    xref = {db: list(xref[db]) for db in xref}
    pop_AF = list(pop_AF.values())
    return xref, pop_AF, variant_annot, colocated_annot


def getPopInfo(annot_key):
    """
    @summary: Returns source and name of a sub-population studied in large genomic programs (1KG, ExAC, ...) from the annotation tag.
    @param annot_key: [str] The tag used in annotation to store AF of the variant in sub-population.
    @return: [list] The source (project name) and the name of the sub-population.
    """
    source = None
    name = None
    if annot_key == "AF":
        source = "1KG"
        name = "Global"
    elif annot_key.lower().startswith("exac_") or annot_key.lower().startswith("gnomAD_") or annot_key.lower().startswith("1kg_") or annot_key.lower().startswith("esp_"):  # ExAC_AF, ExAC_Adj_AF, ExAC_AFR_AF, ExAC_AMR_AF, ...
        source = annot_key.split("_")[0]
        name = "Global" if annot_key.count('_') == 1 else annot_key.split("_")[1]
    elif annot_key.count('_') == 1:  # AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, ...
        source = "1KG"
        name = annot_key.split("_")[0]
    else:
        warnings.warn('The population information stored with tag "{}" cannot be used by "{}".'.format(annot_key, sys.argv[0]))
    return source, name


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Converts VCF annotated with VEP in JSON format.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file file containing variants and annotated with VEP v88+ (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='The path to the file outputted file (format: JSON).')
    args = parser.parse_args()

    # Convert VCF in python dict
    json_data = list()
    with VEPVCFIO(args.input_variants) as FH_vcf:
        for record in FH_vcf:
            for idx_alt in range(len(record.alt)):
                allele_record = getAlleleRecord(FH_vcf, record, idx_alt)
                allele_record.standardizeSingleAllele()
                curr_json = dict()
                # Core information
                curr_json["core"] = {
                    "region": allele_record.chrom,
                    "pos": allele_record.pos,
                    "ref": allele_record.ref,
                    "alt": allele_record.alt[0]
                }
                # Support information
                curr_json["support"] = {
                    "filters": allele_record.filter,
                    "qual": allele_record.qual,
                    "libraries": [{"alt_depth": allele_record.getAD(library)[0], "depth": allele_record.getDP(library), "name": library} for library in FH_vcf.samples]
                }
                if "CSQ" in allele_record.info:
                    # Identical known variants, AF in populations and annotations
                    curr_json["xref"], curr_json["pop_AF"], curr_json["annot"], curr_json["colocated_annot"] = getAnnotSummary(allele_record, record.alt[idx_alt])
                json_data.append(curr_json)

    # Write output file
    with open(args.output_variants, "w") as FH_out:
        FH_out.write(
            json.dumps(json_data, default=lambda o: o.__dict__, sort_keys=True)
        )
