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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse
import warnings
from anacore.region import Region, RegionList, splittedByRef, iterOverlappedByRegion
from anacore.genomicRegion import Intron
from anacore.bed import getAreas
from anacore.gtf import loadModel
from anacore.vcf import VCFIO
from anacore.gff import GFF3IO, GFF3Record


########################################################################
#
# FUNCTIONS
#
########################################################################
def getTargets(in_aln, in_targets=None):
    """
    Return the list of targeted regions.

    :param in_aln: Path to the alignment file (format: SAM/BAM).
    :type in_aln: str
    :param in_targets: Path to the targeted regions (format: BED). They must not contains any overlap.
    :type in_targets: str
    :return: List of targeted regions.
    :rtype: anacore.region.RegionList
    """
    selected_regions = RegionList()
    if in_targets is None:
        with pysam.AlignmentFile(in_aln, "rb") as FH_bam:
            for ref_info in FH_bam.header["SQ"]:
                selected_regions.append(
                    Region(1, ref_info["LN"], "+", ref_info["SN"], ref_info["SN"])
                )
    else:
        selected_regions = getAreas(in_targets)
        # Check lack of overlap
        selected_regions = sorted(selected_regions, key=lambda x: (x.reference.name, x.start, x.end))
        prev_region = selected_regions[0]
        for curr_region in selected_regions[1:]:
            if curr_region.reference.name == prev_region.reference.name:
                if prev_region.end >= curr_region.start:
                    raise Exception("The regions {} and {} contains an overlap.".format(prev_region, curr_region))
            prev_region = curr_region
    return selected_regions


def shallowFromAlignment(aln_path, selected_regions, depth_mode, min_depth, log):
    """
    Return the list of shallow regions from the alignment file.

    :param aln_path: Path to the alignment file (format: SAM/BAM).
    :type aln_path: str
    :param selected_regions: Targeted regions. They must not contains any overlap between them.
    :type selected_regions: anacore.region.RegionList
    :param depth_mode: How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once).
    :type depth_mode: str
    :param min_depth: All the locations with a depth under this value are reported in shallows areas.
    :type min_depth: int
    :param log: Logger of the script.
    :type log: logging.Logger
    :return: List of shallow regions.
    :rtype: anacore.region.RegionList
    """
    shallow = RegionList()
    nb_selected_regions = len(selected_regions)
    idx_in_part = 1
    with pysam.AlignmentFile(aln_path, "rb") as FH_bam:
        for idx_region, region in enumerate(selected_regions):
            if idx_in_part > nb_selected_regions / 10:
                idx_in_part = 0
                log.info("Processed regions {}/{}.".format(idx_region + 1, nb_selected_regions))
            idx_in_part += 1
            prev_opened = {"start": None, "end": None}
            for pileupcolumn in FH_bam.pileup(region.reference.name, region.start - 1, region.end - 1, max_depth=100000000):
                if pileupcolumn.reference_pos + 1 >= region.start and pileupcolumn.reference_pos + 1 <= region.end:
                    curr_reads_depth = 0
                    curr_frag = set()
                    for pileupread in pileupcolumn.pileups:
                        if pileupcolumn.reference_pos + 1 < region.start or pileupcolumn.reference_pos + 1 > region.end:
                            raise Exception("The reference position {}:{} is out of target {}.".format(region.reference.name, pileupcolumn.reference_pos, region))
                        if not pileupread.alignment.is_secondary and not pileupread.alignment.is_duplicate and not pileupread.is_refskip:
                            curr_reads_depth += 1
                            curr_frag.add(pileupread.alignment.query_name)
                    curr_depth = curr_reads_depth
                    if depth_mode == "fragment":
                        curr_depth = len(curr_frag)
                    if min_depth > curr_depth:
                        if prev_opened["start"] is None:
                            prev_opened = {"start": pileupcolumn.reference_pos, "end": pileupcolumn.reference_pos}
                        else:
                            if prev_opened["end"] == pileupcolumn.reference_pos - 1:
                                prev_opened["end"] = pileupcolumn.reference_pos
                            else:
                                shallow.append(
                                    Region(prev_opened["start"] + 1, prev_opened["end"] + 1, "+", region.reference)
                                )
                                prev_opened = {"start": pileupcolumn.reference_pos, "end": pileupcolumn.reference_pos}
            if prev_opened["start"] is not None:
                shallow.append(
                    Region(prev_opened["start"] + 1, prev_opened["end"] + 1, "+", region.reference)
                )
    return shallow


def variantsRegionFromVCF(vcf_path):
    """
    Return the region object corresponding to the known variants in a VCF.

    :param vcf_path: Path to the variants file (format: VCF).
    :type vcf_path: str
    :return: List of variants regions.
    :rtype: anacore.region.RegionList
    .. warning:: The insertion regions are extended to one nucleotid before start and one nucleotid after end.
    """
    variants_region = RegionList()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with VCFIO(vcf_path) as FH_in:
            for record in FH_in:
                variants_region.append(
                    Region(
                        int(record.refStart()),  # Insertions between two bases impact the two bases
                        int(record.refEnd() + 0.5),  # Insertions between two bases impact the two bases
                        None,
                        record.chrom,
                        record.id,
                        {
                            "id": record.id,
                            "gene": ("" if "GENE" not in record.info else record.info["GENE"]),
                            "HGVSp": ("" if "HGVSp" not in record.info else record.info["AA"]),
                            "HGVSc": ("" if "CDS" not in record.info else record.info["CDS"]),
                        }
                    )
                )
    return variants_region


def setVariantsByOverlap(queries, variants):
    """
    Annotate each query by the list of variants overlapping them.

    :param queries: Regions to annotate.
    :type queries: anacore.region.Region
    :param variants: The list of variants where overlapped variants will be searched.
    :type variants: anacore.region.RegionList
    """
    variants_by_chr = splittedByRef(variants)
    queries_by_chr = splittedByRef(queries)
    for chrom, curr_query, overlapped_subjects in iterOverlappedByRegion(queries_by_chr, variants_by_chr):
        curr_query.annot["VAR"] = []
        for sbjct in overlapped_subjects:
            curr_query.annot["VAR"].append(sbjct)


def setTranscriptsAnnotByOverlap(queries, transcripts):
    """
    Annotate each query by the information coming from the transcripts overlapping them.

    :param region: Regions to annotate.
    :type region: anacore.region.Region
    :param transcripts: The list of transcripts where overlapped transcripts will be searched.
    :type transcripts: anacore.region.RegionList
    """
    transcripts_by_chr = splittedByRef(transcripts)
    queries_by_chr = splittedByRef(queries)
    for chrom, curr_query, overlapped_subjects in iterOverlappedByRegion(queries_by_chr, transcripts_by_chr):
        curr_query.annot["ANN"] = getTranscriptsAnnot(curr_query, overlapped_subjects)


def getTranscriptsAnnot(region, transcripts):
    """
    Return for each overlapped transcript the location of start and end of the query region on the transcript and the protein.

    :param region: The query region.
    :type region: anacore.region.Region
    :param transcripts: List of transcripted overlapped by the query region.
    :type transcripts: anacore.region.RegionList
    :return: List of annotations (one by transcript).
    :rtype: list
    """
    annotations = []
    for curr_tr in transcripts:
        curr_annot = {
            "SYMBOL": curr_tr.parent.name,
            "Gene": curr_tr.parent.annot["id"],
            "Feature": curr_tr.annot["id"],
            "Feature_type": "Transcript",
            "STRAND": None,
            "start_EXON": None,
            "start_INTRON": None,
            "start_Protein_position": None,
            "end_EXON": None,
            "end_INTRON": None,
            "end_Protein_position": None
        }
        # Overlap on upstream
        overlap_start = {
            "tr_ref_pos": region.start,
            "tr_sub_idx": None,
            "tr_sub_type": None,
            "prot_pos": None
        }
        if region.start < curr_tr.start:  # The region starts before the transcript and overlap the transcript
            overlap_start["tr_ref_pos"] = curr_tr.start
            overlap_start["tr_sub_type"] = "EXON"
            overlap_start["tr_sub_idx"] = "{}/{}".format(
                (1 if curr_tr.strand != "-" else len(curr_tr.children)),
                len(curr_tr.children)
            )
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region):
                protein = curr_tr.proteins[0]
                overlap_start["prot_pos"] = (1 if curr_tr.strand != "-" else protein.aaLength())
        else:
            subregion, subregion_idx = curr_tr.getSubFromRefPos(region.start)
            if issubclass(subregion.__class__, Intron):  # The region starts in an intron
                # Get first pos of next exon
                downstream_exon_idx = subregion_idx + 1 - 1  # 0-based
                if curr_tr.strand == "-":
                    downstream_exon_idx = subregion_idx - 1  # 0-based
                overlap_start["tr_ref_pos"] = curr_tr.children[downstream_exon_idx].start
                overlap_start["tr_sub_type"] = "INTRON"
                overlap_start["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children) - 1)
            else:
                overlap_start["tr_sub_type"] = "EXON"
                overlap_start["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children))
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region):
                protein = curr_tr.proteins[0]
                if overlap_start["tr_ref_pos"] < protein.start:  # The region overlap an UTR
                    overlap_start["prot_pos"] = (1 if curr_tr.strand != "-" else protein.aaLength())
                else:
                    overlap_start["prot_pos"] = protein.getPosOnRegion(overlap_start["tr_ref_pos"])[0]
        # Overlap on downstream
        overlap_end = {
            "tr_ref_pos": region.end,
            "tr_sub_idx": None,
            "tr_sub_type": None,
            "prot_pos": None
        }
        if region.end > curr_tr.end:  # The region ends after the transcript and overlap the transcript
            overlap_end["tr_ref_pos"] = curr_tr.end
            overlap_end["tr_sub_type"] = "EXON"
            overlap_end["tr_sub_idx"] = "{}/{}".format(
                (len(curr_tr.children) if curr_tr.strand != "-" else 1),
                len(curr_tr.children)
            )
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region):
                protein = curr_tr.proteins[0]
                overlap_end["prot_pos"] = (protein.aaLength() if curr_tr.strand != "-" else 1)
        else:
            subregion, subregion_idx = curr_tr.getSubFromRefPos(region.end)
            if issubclass(subregion.__class__, Intron):  # The region ends in an intron
                # Get last pos of previous exon
                upstream_exon_idx = subregion_idx - 1  # 0-based
                if curr_tr.strand == "-":
                    upstream_exon_idx = subregion_idx + 1 - 1  # 0-based
                overlap_end["tr_ref_pos"] = curr_tr.children[upstream_exon_idx].end
                overlap_end["tr_sub_type"] = "INTRON"
                overlap_end["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children) - 1)
            else:
                overlap_end["tr_sub_type"] = "EXON"
                overlap_end["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children))
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region):
                protein = curr_tr.proteins[0]
                if overlap_end["tr_ref_pos"] > protein.end:  # The region overlap an UTR
                    overlap_end["prot_pos"] = (protein.aaLength() if curr_tr.strand != "-" else 1)
                else:
                    overlap_end["prot_pos"] = protein.getPosOnRegion(overlap_end["tr_ref_pos"])[0]
        # Store info in annotations
        start = overlap_start
        end = overlap_end
        curr_annot["STRAND"] = "1"
        if curr_tr.strand == "-":
            curr_annot["STRAND"] = "-1"
            start = overlap_end
            end = overlap_start
        curr_annot["start_" + start["tr_sub_type"]] = start["tr_sub_idx"]
        curr_annot["start_Protein_position"] = start["prot_pos"]
        curr_annot["end_" + end["tr_sub_type"]] = end["tr_sub_idx"]
        curr_annot["end_Protein_position"] = end["prot_pos"]
        annotations.append(curr_annot)
    return annotations


def writeOutput(out_path, shallow):
    """
    Write shallow areas there annotations in a GFF3 file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param shallow: The list of shallow areas.
    :type shallow: anacore.region.RegionList
    """
    with GFF3IO(out_path, "w") as FH_out:
        for curr_shallow in sorted(shallow, key=lambda x: (x.reference.name, x.start, x.end)):
            record = GFF3Record(
                curr_shallow.reference.name,
                "depthAnalysis",
                "experimental_feature",
                curr_shallow.start,
                curr_shallow.end
            )
            if args.input_annotations is not None:
                for idx, annot in enumerate(curr_shallow.annot["ANN"]):
                    fields = []
                    for k, v in sorted(annot.items()):
                        fields.append("{}:{}".format(k, v))
                    record.annot["ann_{}".format(idx + 1)] = "|".join(fields)
            if len(args.inputs_variants) > 0:
                for idx, var_region in enumerate(curr_shallow.annot["VAR"]):
                    fields = []
                    for k, v in sorted(var_region.annot.items()):
                        fields.append("{}:{}".format(k, v))
                    record.annot["var_{}".format(idx + 1)] = "|".join(fields)
            FH_out.write(record)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Extract shallow areas from the alignment are annotate them with genomic features and known variants.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-m', '--depth-mode', choices=["read", "fragment"], default="fragment", help='How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once). [Default: %(default)s]')
    parser.add_argument('-d', '--min-depth', type=int, default=30, help='All the locations with a depth under this value are reported in shallows areas. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-b', '--input-aln', required=True, help='Path to the variants files (format: BAM).')
    group_input.add_argument('-t', '--input-targets', help='Path to the targeted regions (format: BED). They must not contains any overlap. [Default: all positions defined in the alignment file header]')
    group_input.add_argument('-a', '--input-annotations', help='Path to the file defining transcripts, genes and proteins locations (format: GTF). This file allow to annotate locations on genes and proteins located on shallows areas. [Default: The shallows areas are not annotated]')
    group_input.add_argument('-s', '--inputs-variants', nargs="+", help='Path(es) to the file(s) defining known variants (format: VCF). This file allow to annotate variant potentially masked because they are on shallows areas. [Default: The variants on shallows areas are not reported]')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-shallow', default="shallow_areas.gff3", help='Path to the file containing shallow areas and there annotations. (format: GFF3). [Default: %(default)s]')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load selected regions
    log.info("Load targeted regions.")
    selected_regions = getTargets(args.input_aln, args.input_targets)

    # Find shallow areas
    log.info("Find shallow areas.")
    shallow = shallowFromAlignment(args.input_aln, selected_regions, args.depth_mode, args.min_depth, log)

    # Annotate shallow areas
    if args.input_annotations is not None:
        log.info("Load annotations from {}.".format(args.input_annotations))
        transcripts = loadModel(args.input_annotations, "transcripts")
        log.info("Annotate shallow areas.")
        setTranscriptsAnnotByOverlap(shallow, transcripts)

    # Retrieved known variants potentialy masked in shallow areas
    for curr_input in args.inputs_variants:
        log.info("Load variants from {}.".format(curr_input))
        variant_regions = variantsRegionFromVCF(curr_input)
        log.info("List potentialy masked mutations.")
        setVariantsByOverlap(shallow, variant_regions)

    # Write output
    log.info("Write output.")
    writeOutput(args.output_shallow, shallow)
