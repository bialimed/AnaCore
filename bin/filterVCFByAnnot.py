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
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.annotVcf import AnnotVCFIO, getAlleleRecord, HeaderFilterAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGeneByNM(gene_to_id_file, trim_version=False):
    """
    Return gene name by RNA_id.

    :param gene_to_id_file: Path to the file describing the link between genes and RNA_id (format: TSV). Each line has the following format: <GENE>\t<RNA_ID>.
    :type gene_to_id_file: str
    :param trim_version: With True the version number is removed from the id.
    :type trim_version: bool
    """
    gene_by_NM = dict()
    with open(gene_to_id_file) as FH_ref:
        for line in FH_ref:
            if not line.startswith("#"):
                gene, NM = [field.strip() for field in line.split("\t")]
                if trim_version:
                    NM = NM.split(".")[0]
                gene_by_NM[NM] = gene
    return gene_by_NM


def tagAnnotRNA(annot, kept_id, trim_version=False):
    """
    Add tag ANN.RBA in annot["FILTER"] if RNA used for annotation is not in selected RNA.

    :param annot: One variant annotation (extracted from list variant.info["ANN"]).
    :type annot: dict
    :param kept_id: The IDs of the selected RNA.
    :type kept_id: dict
    :param trim_version: If True the version number of the RNA is not tacking into account. Otherwise only the selected RNA with the specific version provided in ID are not tagged.
    :type trim_version: bool
    """
    RNA_id = annot["Feature"]
    if RNA_id is not None and trim_version:
        RNA_id = RNA_id.split(".")[0]
    if RNA_id not in kept_id:
        annot["FILTER"].add("ANN.RNA")


def tagAnnotCSQ(annot, valid_consequences):
    """
    Add tag ANN.CSQ in annot["FILTER"] if no consequences are in valid_consequences.

    :param annot: One variant annotation (extracted from list variant.info["ANN"]).
    :type annot: dict
    :param valid_consequences: An annotation with one of these consequences is not tagged.
    :type valid_consequences: list
    """
    is_filtered = True
    if annot["Consequence"] is not None:
        consequences = annot["Consequence"].split("&")  # For multi-consequence allele: Example: 'start_lost&NMD_transcript_variant'
        for current_csq in consequences:
            if current_csq in valid_consequences:
                is_filtered = False
    if is_filtered:
        annot["FILTER"].add("ANN.CSQ")


def tagCollocated(annot, alt_in_annot_format):
    """
    Add tag ANN.COLLOC in annot["FILTER"] if annotation come from collocated variant and not the same variant.

    :param annot: One variant annotation (extracted from list variant.info["ANN"]).
    :type annot: dict
    :param alt_in_annot_format: The alternative allele of the variant in same format as annotation allele.
    :type alt_in_annot_format: str
    """
    if annot["Allele"] != alt_in_annot_format:
        annot["FILTER"].add("ANN.COLLOC")


def tagAnnotPolymophism(annot, checked_pop_tags, min_AF=0.01):
    """
    Add tag ANN.popAF in annot["FILTER"] if the annotation come from a polymorphism.

    :param annot: One variant annotation (extracted from list variant.info["ANN"]).
    :type annot: dict
    :param checked_pop_tags: The evaluated populations.
    :type checked_pop_tags: list
    :param min_AF: If the AF in a evaluated population is superior than this value the variant is take as polymorphism in this population.
    :type min_AF: float
    """
    is_polymorphism = False
    for pop_freq in checked_pop_tags:
        if pop_freq in annot and annot[pop_freq] is not None:
            frequencies = annot[pop_freq].split("&")  # For multi-reference SNP. Example: {'Existing_variation': 'rs7367494&rs34188929', 'Gene': 'ENSG00000184908', 'Allele': 'T', ..., 'ExAC_SAS_AF': '0.9253&0.9253', 'AF': '0.8299', 'EAS_AF': '0.9980&0.9980'}
            for current_freq in frequencies:
                if float(current_freq) >= min_AF:
                    is_polymorphism = True
    if is_polymorphism:
        annot["FILTER"].add("ANN.popAF")


def writeHeader(FH_in, FH_out, args):
    """
    Write VCF header.

    :param FH_in: File handle to model file.
    :type FH_in: file object
    :param FH_out: File handle to output file.
    :type FH_out: file object
    :param args: Scripts arguments.
    :type args: NameSpace
    """
    FH_out.copyHeader(FH_in)
    if "FILTER" not in FH_out.ANN_titles:
        FH_out.ANN_titles.append("FILTER")
    FH_out.filter["popAF"] = HeaderFilterAttr("popAF", "The variant is present with more of " + str(args.polym_threshold * 100) + "% in one of the following population: '" + "' ".join(args.polym_populations) + "'.")
    FH_out.filter["ANN.popAF"] = HeaderFilterAttr("ANN.popAF", "The variant is present with more of " + str(args.polym_threshold * 100) + "% in one of the following population: '" + "' ".join(args.polym_populations) + "'.")  # Is distinct of popAF because annotations can contain collocated variants
    if args.input_selected_RNA is not None:
        FH_out.filter["ANN.RNA"] = HeaderFilterAttr("ANN.RNA", "The annotation RNA is not one of the selected ({}).".format(args.input_selected_RNA))
    FH_out.filter["CSQ"] = HeaderFilterAttr("CSQ", "The variant has no consequence corresponding at one in the following list: '" + "' ".join(args.kept_consequences) + "'.")
    FH_out.filter["ANN.CSQ"] = HeaderFilterAttr("ANN.CSQ", "The annotation consequence does not correspond at one in the following list: '" + "' ".join(args.kept_consequences) + "'.")
    FH_out.filter["ANN.COLLOC"] = HeaderFilterAttr("ANN.COLLOC", "The annotation corresponds to a collocated variant.")
    FH_out.writeHeader()


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filters variants and their annotations on annotations. In "remove" mode the annotations are deleted if they not fit criteria and the variant is removed if none of his annotations fit criterias.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag "CSQ" and/or "popAF" is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]')
    group_filter.add_argument('-p', '--polym-populations', default=["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "AA_AF", "EA_AF", "ExAC_AF", "ExAC_Adj_AF", "ExAC_AFR_AF", "ExAC_AMR_AF", "ExAC_EAS_AF", "ExAC_FIN_AF", "ExAC_NFE_AF", "ExAC_OTH_AF", "ExAC_SAS_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"], help='Populations frequencies used as reference for polymorphism detection. [Default: %(default)s]')
    group_filter.add_argument('-l', '--polym-threshold', type=float, default=0.01, help='Minimum frequency in population to tag allele as polymorphism. [Default: %(default)s]')
    group_filter.add_argument('-k', '--kept-consequences', default=["TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant"], nargs='+', help='The variants without one of these consequences are tagged as CSQ (see http://www.ensembl.org/info/genome/variation/predicted_data.html). [Default: %(default)s]')
    group_filter.add_argument('-r', '--input-selected-RNA', help='The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.')
    group_filter.add_argument('-w', '--rna-without-version', action='store_true', help='With this option the version number of the reference RNA is not used in filter.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file containing annotated variants (format: VCF). CAUTION: The annotations produced by VEP must be corrected to be coherent between ALT and INFO.annotation.Allele (by default VEP normalizes INFO.annotation.Allele but not ALT).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the filtered file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    kept_ID = {}
    if args.input_selected_RNA is not None:
        kept_ID = getGeneByNM(args.input_selected_RNA, args.rna_without_version)
    with AnnotVCFIO(args.input_variants, "r", args.annotation_field) as FH_in:
        with AnnotVCFIO(args.output_variants, "w") as FH_out:
            # Header
            writeHeader(FH_in, FH_out, args)
            # Records
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord(FH_in, record, alt_idx)
                    # Evaluates annotations
                    annot_pass = False
                    record_is_filtered_on_polym = False
                    record_is_filtered_on_csq = True
                    for annot_idx, annot in enumerate(alt_record.info[FH_in.annot_field]):
                        old_filters = set()
                        if "FILTER" in annot and annot["FILTER"] is not None:
                            if annot["FILTER"] != "PASS":
                                old_filters = set(annot["FILTER"].split(","))
                        annot["FILTER"] = set()
                        # Not same variant
                        tagCollocated(annot, alt)
                        # Polymorphism
                        tagAnnotPolymophism(annot, args.polym_populations, args.polym_threshold)
                        if "ANN.COLLOC" not in annot["FILTER"] and "ANN.popAF" in annot["FILTER"]:  # The variant is not a collocated and is polymorphism
                            record_is_filtered_on_polym = True
                        # Reference RNA
                        if args.input_selected_RNA is not None:
                            tagAnnotRNA(annot, kept_ID, args.rna_without_version)
                        # Consequences on RNA
                        tagAnnotCSQ(annot, args.kept_consequences)
                        if "ANN.COLLOC" not in annot["FILTER"] and "ANN.RNA" not in annot["FILTER"] and "ANN.CSQ" not in annot["FILTER"]:
                            record_is_filtered_on_csq = False
                        # Manage FILTER tag
                        if len(annot["FILTER"]) == 0:
                            annot_pass = True
                        new_filters = annot["FILTER"] | old_filters
                        annot["FILTER"] = "&".join(sorted(new_filters))
                        if annot["FILTER"] == "":
                            annot["FILTER"] = "PASS"
                    # Filter record
                    if args.mode == "tag":
                        if alt_record.filter is None or len(alt_record.filter) == 0 or alt_record.filter[0] == "PASS":
                            alt_record.filter = list()
                        if not annot_pass:
                            if record_is_filtered_on_csq:
                                alt_record.filter.append("CSQ")
                            if record_is_filtered_on_polym:
                                alt_record.filter.append("popAF")
                        elif len(alt_record.filter) == 0:
                            alt_record.filter.append("PASS")
                        FH_out.write(alt_record)
                    elif annot_pass:
                        # Delete filtered annot
                        delete_annot_ix = list()
                        for annot_idx, annot in enumerate(alt_record.info[FH_in.annot_field]):
                            if annot["FILTER"] != "PASS":
                                delete_annot_ix.append(annot_idx)
                        for curr_idx in reversed(delete_annot_ix):
                            del(alt_record.info[FH_in.annot_field][curr_idx])
                        # Write record
                        FH_out.write(alt_record)
    log.info("End of job")
