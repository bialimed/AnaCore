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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import warnings
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.bed import BEDIO
from anacore.region import Region, RegionList, splittedByRef
from anacore.genomicRegion import Transcript, Gene, Exon
from anacore.gff import GFF3IO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGeneByRefTr(in_ref_tr, trim_version=True):
    """
    Return gene object by RNA_id.

    :param in_ref_tr: Path to the file describing the link between genes and RNA_id (format: TSV). Each line has the following format: <GENE>\t<RNA_ID>.
    :type in_ref_tr: str
    :param trim_version: With True the version number is removed from the id.
    :type trim_version: boolean
    :return: The genes by RNA_id.
    :rtype: dict
    """
    gene_by_tr = dict()
    with open(in_ref_tr) as FH_tr:
        for line in FH_tr:
            if not line.startswith("#"):
                gene, tr_id = [field.strip() for field in line.split("\t")]
                if tr_id != "":
                    if trim_version:
                        tr_id = tr_id.split(".")[0]
                    gene_by_tr[tr_id] = Gene(name=gene)
    return gene_by_tr


def getTranscriptAnnot(in_annot, gene_by_tr):
    """
    Get genomic model (genes, transcripts and exons) for the selected transcripts.

    :param in_annot: Path to the genomic annotations (format: GFF3).
    :type in_annot: str
    :param gene_by_tr: Gene by selected transcripts.
    :type gene_by_tr: dict
    :return: The list of selected transcripts.
    :rtype: anacore.region.RegionList
    """
    tr_by_id = dict()
    with GFF3IO(in_annot) as FH_annot:
        for record in FH_annot:
            if record.type == "mRNA" and "transcript_id" in record.annot:
                tr_id = record.annot["transcript_id"]
                tr_id = tr_id.split(".")[0]  # Remove transcript version
                if tr_id in gene_by_tr:  # Transcript is in panel
                    if tr_id not in tr_by_id:
                        tr_by_id[tr_id] = Transcript(
                            record.start,
                            record.end,
                            record.strand,
                            record.seq_id,
                            tr_id,
                            {},
                            gene_by_tr[tr_id]
                        )
            if record.type == "exon" and "transcript_id" in record.annot:
                tr_id = record.annot["transcript_id"]
                tr_id = tr_id.split(".")[0]  # Remove transcript version
                if tr_id in gene_by_tr:  # Transcript is in panel
                    # Store the exon
                    tr_by_id[tr_id].addChild(
                        Exon(
                            record.start,
                            record.end,
                            record.strand,
                            record.seq_id
                        )
                    )
    if len(gene_by_tr) != len(tr_by_id):
        raise Exception("The following transcripts are missing in {}: {}".format(
            args.input_annotation,
            set(gene_by_tr.keys()).difference(set(tr_by_id.keys()))
        ))
    return RegionList(tr_by_id.values())


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Rename the BED's entries with annotations information (gene and exon) for the selected transcripts.")
    parser.add_argument('-c', '--is-thick-based', action='store_true', help='With this option the annotation is only based on thick coordinates of the entry if they exist.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-t', '--input-reference-tr', required=True, help='Path to the file containing the link between selected genes and transcripts (format: TSV). Each line contains 2 fields: gene_name and transcript_id.')
    group_input.add_argument('-a', '--input-annotation', required=True, help='Path to the file containing the genomic annotations (format: GFF3).')
    group_input.add_argument('-r', '--input-regions', required=True, help='Path to the file where the regions will be renamed (format: BED).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-regions', default="renamed.bed", help='Path to the file containing the renamed regions (format: BED). [Default: %(default)s]')
    args = parser.parse_args()

    # Get transcripts
    gene_by_tr = getGeneByRefTr(args.input_reference_tr)
    selected_transcripts = getTranscriptAnnot(args.input_annotation, gene_by_tr)
    tr_by_chr = splittedByRef(selected_transcripts)
    # Write renamed regions
    out_nb_col = BEDIO.getMaxNbCol(args.input_regions)
    if out_nb_col == 3:
        out_nb_col = 4
    with BEDIO(args.input_regions) as FH_regions:
        with BEDIO(args.output_regions, "w", out_nb_col) as FH_out:
            for record_idx, record in enumerate(FH_regions):
                target = Region(record.start, record.end, record.strand, record.chrom)
                if args.is_thick_based and record.thickStart is not None and record.thickEnd is not None:
                    target.start = record.thickStart
                    target.end = record.thickEnd
                overlapped_tr = list()
                if record.chrom in tr_by_chr:
                    overlapped_tr = tr_by_chr[record.chrom].getOverlapped(target)
                if len(overlapped_tr) > 1:
                    warnings.warn("The region {} overlaps several transcripts ({}).".format(target, [str(tr) for tr in overlapped_tr]))
                if len(overlapped_tr) >= 1:
                    overlapped_exons = overlapped_tr[0].children.getOverlapped(target)
                    features = list()
                    for curr_feature in overlapped_exons:
                        features.append("ex{}".format(curr_feature.annot["siblings_idx"]))
                    record.name = "{}_{}_{}_{}".format(
                        overlapped_tr[0].parent.name,
                        "-".join(features),
                        ("" if target.strand is None else target.strand),
                        "trgt" + str(record_idx)
                    )
                else:
                    record.name = "{}__{}_{}".format(
                        "noFeature",
                        ("" if target.strand is None else target.strand),
                        "trgt" + str(record_idx)
                    )
                FH_out.write(record)
