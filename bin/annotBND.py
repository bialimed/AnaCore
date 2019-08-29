#!/usr/bin/env python3
#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.gtf import loadModel
from anacore.genomicRegion import Intron
from anacore.annotVcf import AnnotVCFIO, VCFIO
from anacore.region import Region, splittedByRef


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGeneAnnot(record, genes_by_chr):
    """
    Return genomic items overlapped by the BND record.

    :param record: The BND record.
    :type record: anacore.vcf.VCFRecord
    :param genes_by_chr: By chromosomes a tree where nodes are genes, transcripts, protein, exons and CDS.
    :type genes_by_chr: dict
    :return: The list of annotations (one annotation by overlapped transcript).
    :rtype: list
    """
    variant_region = Region(record.pos, None, None, record.chrom, record.getName())
    annotations = []
    if record.chrom not in genes_by_chr:
        log.warn("The region {} where the breakpoint {} has been detected does not contain any annotation.".format(record.chrom, record.id))
    else:
        container_genes = genes_by_chr[record.chrom].getContainers(variant_region)
        for curr_gene in container_genes:
            container_transcripts = curr_gene.children.getContainers(variant_region)
            if len(container_transcripts) == 0:
                log.warn("The breakpoint {} is contained by gene {} but by 0 of these transcripts.".format(variant_region, curr_gene))
            else:
                for curr_transcript in container_transcripts:
                    curr_annot = {
                        "SYMBOL": curr_gene.name,
                        "Gene": curr_gene.annot["id"],
                        "Feature": curr_transcript.annot["id"],
                        "Feature_type": "Transcript",
                        "STRAND": None,
                        "EXON": None,
                        "INTRON": None,
                        "CDS_position": None,
                        "Protein_position": None,
                        "Codon_position": None
                    }
                    if curr_transcript.strand is not None:
                        curr_annot["STRAND"] = ("1" if curr_transcript.strand == "+" else "-1")
                    subregion, subregion_idx = curr_transcript.getSubFromRefPos(variant_region.start)
                    if issubclass(subregion.__class__, Intron):  # On intron
                        curr_annot["INTRON"] = "{}/{}".format(
                            subregion_idx,
                            len(curr_transcript.children) - 1
                        )
                    else:  # On exon
                        curr_annot["EXON"] = "{}/{}".format(
                            subregion_idx,
                            len(curr_transcript.children)
                        )
                        if len(curr_transcript.proteins) > 1:
                            log.error(
                                "The management of several proteins for one transcript is not implemented. The transcript {} contains several proteins {}.".format(curr_transcript, curr_transcript.proteins),
                                exec_info=True
                            )
                        if len(curr_transcript.proteins) > 0:
                            curr_annot["CDS_position"] = curr_transcript.proteins[0].getNtPosFromRefPos(variant_region.start)
                            if curr_annot["CDS_position"] is not None:
                                curr_annot["Protein_position"], curr_annot["Codon_position"] = curr_transcript.proteins[0].getPosOnRegion(variant_region.start)
                    annotations.append(curr_annot)
    return annotations


def annotGeneShard(record):
    """
    Add which shard of genes are in fusion (up or down).

    :param record: The annotated BND record. The BND is previously annotated by genomic regions (see getRegionGeneAnnot()).
    :type record: anacore.vcf.VCFRecord
    """
    # Get position relative to the break
    is_before_break = []
    for bnd_idx, alt in enumerate(record.alt):
        if alt.startswith("[") or alt.startswith("]"):
            is_before_break.append(False)
        else:
            is_before_break.append(True)
    if len(set(is_before_break)) > 1:
        record_name = record.id if record.id is not None else record.getName()
        log.error(
            "The breakend {} has several fusion partners with different break's configuration.".format(record_name),
            exec_info=True
        )
    is_before_break = is_before_break[0]
    # Set BND_stream
    for annot in record.info["ANN"]:
        annot["GENE_SHARD"] = None
        if annot["STRAND"] is not None:
            annot["GENE_SHARD"] = "down"
            if annot["STRAND"] == "1":
                if is_before_break:
                    annot["GENE_SHARD"] = "up"
            else:
                if not is_before_break:
                    annot["GENE_SHARD"] = "up"


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Annotate BND in a VCF with content of a GTF.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-annotations', required=True, help='Path to the file containing the annotations of genes and transcript for the reference used in variant calling. (format: GTF).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the file containing variants. (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the annotated file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load annotations
    log.info("Load model from {}.".format(args.input_annotations))
    genes = loadModel(args.input_annotations, "genes")
    genes_by_chr = splittedByRef(genes)

    # Annot variants
    log.info("Annot variants in {}.".format(args.input_variants))
    with AnnotVCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.ANN_titles = ["SYMBOL", "Gene", "Feature", "Feature_type", "STRAND", "EXON", "INTRON", "CDS_position", "Protein_position", "GENE_SHARD"]
            FH_out.info[args.annotation_field] = {
                "type": str,
                "type_tag": "String",
                "number": None,
                "number_tag": ".",
                "description": "Consequence annotations. Format: " + "|".join(FH_out.ANN_titles)
            }
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                if record.info["SVTYPE"] == "BND":
                    record.info[args.annotation_field] = getGeneAnnot(record, genes_by_chr)
                    annotGeneShard(record)
                FH_out.write(record)
    log.info("End of job")
