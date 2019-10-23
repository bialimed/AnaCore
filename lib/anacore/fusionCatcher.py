# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing fusionCatcher's output. Each record is converted to VCFRecord."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import uuid
from anacore.sv import HashedSVIO
from anacore.vcf import VCFRecord, HeaderInfoAttr, HeaderFormatAttr


def getCoordDictFromCoordStr(coord):
    """
    Return dict from elements of coordinates defined by fusionCatcher: "CHROM:POS:STRAND" to {"chrom": "CHROM", "pos": "POS", "strand":"STRAND"}.

    :param coord: Coordinates defined by fusionCatcher in format "CHROM:POS:STRAND" (see: Fusion_point_for_gene_1(5end_fusion_partner), Fusion_point_for_gene_2(3end_fusion_partner)).
    :type coord: str
    :return: Coordinates in dict format {"chrom": "CHROM", "pos": "POS", "strand":"STRAND"}.
    :rtype: dict
    """
    matches = re.match("(.+):(\d+):(.)", coord)
    if matches is None:
        raise Exception("The coordinates fields {} cannot be parsed.".format(coord))
    return {
        "chrom": matches.group(1),
        "pos": matches.group(2),
        "strand": matches.group(3)
    }


def getAltFromCoord(first_coord, second_coord):
    """
    Return VCF alternative field value for the breakend of the 5' and 3' shard of fusion.

    :param first_coord: Coordinates of the breakend for the first 5' shard in fusion.
    :type first_coord: dict
    :param second_coord: Coordinates of the breakend for the first 3' shard in fusion.
    :type second_coord: dict
    :return: Alternative for the first shard in fusion and alternative for the first shard in fusion.
    :rtype: (dict, dict))
    """
    first_alt = None
    second_alt = None
    if first_coord["strand"] == "+":
        if second_coord["strand"] == "-":  # First shard is + and second is -
            first_alt = "N]{}:{}]".format(second_coord["chrom"], second_coord["pos"])
            second_alt = "N]{}:{}]".format(first_coord["chrom"], first_coord["pos"])
        else:   # First shard is + and second is +
            first_alt = "N[{}:{}[".format(second_coord["chrom"], second_coord["pos"])
            second_alt = "]{}:{}]N".format(first_coord["chrom"], first_coord["pos"])
    else:
        if second_coord["strand"] == "-":  # First shard is - and second is -
            first_alt = "]{}:{}]N".format(second_coord["chrom"], second_coord["pos"])
            second_alt = "N[{}:{}[".format(first_coord["chrom"], first_coord["pos"])
        else:  # First shard is - and second is +
            first_alt = "[{}:{}[N".format(second_coord["chrom"], second_coord["pos"])
            second_alt = "[{}:{}[N".format(first_coord["chrom"], first_coord["pos"])
    return first_alt, second_alt


def setFusionCatcherVCFHeader(vcf_io, annotation_field="FCANN"):
    """
    Set header of an AnnotVCFIO to write FusionCatcher records in VCF output. This function pre-set ANN_titles, info and format used by FusionCatcherIO records.

    :param vcf_io: The VCF used to write FusionCatcher records.
    :type vcf_io: anacore.annotVcf.AnnotVCFIO
    :param annotation_field: INFO field used for store annotations.
    :type annotation_field: str
    """
    # INFO
    vcf_io.info = {
        "MATEID": HeaderInfoAttr("MATEID", type="String", number="A", description="ID of mate breakend."),
        "SVTYPE": HeaderInfoAttr("SVTYPE", type="String", number="1", description="Type of structural variant."),
        "SOURCES": HeaderInfoAttr("SOURCES", type="String", number=".", description="Aligning method used for mapping the reads and finding the fusion genes. Here are two methods used which are: (i) BOWTIE = only Bowtie aligner is used for mapping the reads on the genome and exon-exon fusion junctions, (ii) BOWTIE+BLAT = Bowtie aligner is used for mapping reads on the genome and BLAT is used for mapping reads for finding the fusion junction, (iii) BOWTIE+STAR = Bowtie aligner is used for mapping reads on the genome and STAR is used for mapping reads for finding the fusion junction, (iv) BOWTIE+BOWTIE2 = Bowtie aligner is used for mapping reads on the genome and Bowtie2 is used for mapping reads for finding the fusion junction."),
        "RNA_FIRST": HeaderInfoAttr("RNA_FIRST", type="Flag", number="0", description="For RNA fusions, this break-end is 5' in the fusion transcript."),
        "RNA_CONTIG": HeaderInfoAttr("RNA_CONTIG", type="String", number="1", description="The sequence of the breakend spanning contig."),
        annotation_field: HeaderInfoAttr(annotation_field, type="String", number=".", description="Consequence annotations. Format: SYMBOL|Gene|EXON|Effect.")
    }
    # ANN_titles
    vcf_io.ANN_titles = ["SYMBOL", "Gene", "EXON", "Effect"]
    # FORMAT
    vcf_io.format = {
        "PR": HeaderFormatAttr("PR", type="Integer", number="1", description="Count of pairs of reads supporting the fusion (including also the multimapping reads)."),
        "SR": HeaderFormatAttr("SR", type="Integer", number="1", description="Count of unique reads (i.e. unique mapping positions) mapping on the fusion junction. Shortly, here are counted all the reads which map on fusion junction minus the PCR duplicated reads."),
        "CM": HeaderFormatAttr("CM", type="Integer", number="1", description="Count of reads mapping simultaneously on both genes which form the fusion gene. This is an indication how similar are the DNA/RNA sequences of the genes forming the fusion gene (i.e. what is their homology because highly homologous genes tend to appear show as candidate fusion genes). In case of completely different sequences of the genes involved in forming a fusion gene then here it is expected to have the value zero."),
        "LA": HeaderFormatAttr("LA", type="Integer", number="1", description="Longest anchor (hangover) found among the unique reads mapping on the fusion junction.")
    }


class FusionCatcherIO(HashedSVIO):
    """Class to read and write fusions in fusionCatcher TSV output format (see https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md#output-data). The records are in anacore.vcf.VCFRecord format and a fusion is represented by 2 records (breakend for the 5' shard and breakend for the 3' shard)."""

    def __init__(self, filepath, mode="r", sample_name="sample", annot_field="FCANN"):
        """
        Build and return an instance of FusionCatcherIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param sample_name: The name of the sample.
        :type sample_name: str
        :param annot_field: The tag for the field used to store annotations.
        :type annot_field: str
        :return: The new instance.
        :rtype: FusionCatcherIO
        """
        super().__init__(filepath, mode, "\t")
        self.annot_field = annot_field
        self.sample_name = sample_name
        if self.titles is None:
            self.titles = [
                "Gene_1_symbol(5end_fusion_partner)",
                "Gene_2_symbol(3end_fusion_partner)",
                "Fusion_description",
                "Counts_of_common_mapping_reads",
                "Spanning_pairs",
                "Spanning_unique_reads",
                "Longest_anchor_found",
                "Fusion_finding_method",
                "Fusion_point_for_gene_1(5end_fusion_partner)",
                "Fusion_point_for_gene_2(3end_fusion_partner)",
                "Gene_1_id(5end_fusion_partner)",
                "Gene_2_id(3end_fusion_partner)",
                "Exon_1_id(5end_fusion_partner)",
                "Exon_2_id(3end_fusion_partner)",
                "Fusion_sequence",
                "Predicted_effect"
            ]

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line.
        :rtype: dict
        """
        fusion_record = super()._parseLine()
        return self._dictToBNDRecords(fusion_record)

    def _dictToBNDRecords(self, fusion_record):
        """
        Return the two breakends records from the fusion record.

        :param fusion_record: The fusion record.
        :type fusion_record: dict
        :return: The breakend record for the 5' shard of fusion and the breakend record for the 3' shard.
        :rtype: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        """
        filter_field = None if fusion_record["Fusion_description"] == "" else [elt.strip() for elt in fusion_record["Fusion_description"].split(",")]
        sample_field = {
            self.sample_name: {
                "PR": fusion_record["Spanning_pairs"],
                "SR": fusion_record["Spanning_unique_reads"],
                "CM": fusion_record["Counts_of_common_mapping_reads"],
                "LA": fusion_record["Longest_anchor_found"]
            }
        }
        format_field = sorted(sample_field[self.sample_name].keys())
        info_sources = None if fusion_record["Fusion_finding_method"] == "" else [elt.strip() for elt in fusion_record["Fusion_finding_method"].split(";")]
        first_id = str(uuid.uuid4())
        first_coord = getCoordDictFromCoordStr(fusion_record["Fusion_point_for_gene_1(5end_fusion_partner)"])
        second_id = str(uuid.uuid4())
        second_coord = getCoordDictFromCoordStr(fusion_record["Fusion_point_for_gene_2(3end_fusion_partner)"])
        first_alt, second_alt = getAltFromCoord(first_coord, second_coord)
        # First breakend
        first_info = {
            "MATEID": [second_id],
            "SVTYPE": "BND",
            "SOURCES": info_sources,
            "RNA_FIRST": True,
            self.annot_field: [{
                "SYMBOL": fusion_record["Gene_1_symbol(5end_fusion_partner)"],
                "Gene": fusion_record["Gene_1_id(5end_fusion_partner)"],
                "EXON": fusion_record["Exon_1_id(5end_fusion_partner)"],
                "Effect": fusion_record["Predicted_effect"]
            }]
        }
        first_bnd = VCFRecord(
            first_coord["chrom"],
            first_coord["pos"],
            first_id,
            "N",  # ref
            [first_alt],  # alt
            pFilter=filter_field,
            info=first_info,
            pFormat=format_field,
            samples=sample_field
        )
        # Second breakends
        second_info = {
            "MATEID": [first_id],
            "SVTYPE": "BND",
            "SOURCES": info_sources,
            "RNA_CONTIG": fusion_record["Fusion_sequence"],
            self.annot_field: [{
                "SYMBOL": fusion_record["Gene_2_symbol(3end_fusion_partner)"],
                "Gene": fusion_record["Gene_2_id(3end_fusion_partner)"],
                "EXON": fusion_record["Exon_2_id(3end_fusion_partner)"],
                "Effect": fusion_record["Predicted_effect"]
            }]
        }
        second_bnd = VCFRecord(
            second_coord["chrom"],
            second_coord["pos"],
            second_id,
            "N",  # ref
            [second_alt],  # alt
            pFilter=filter_field,
            info=second_info,
            pFormat=format_field,
            samples=sample_field
        )
        return [first_bnd, second_bnd]

    # def _BNDRecordsToDict(self, bnd_records):
    #     return None ########################################################################
    #
    # def recordsToLine(self, bnd_records):
    #     """
    #     Return the two breakends as fusion in SV format.
    #
    #     :param bnd_records: The two records of the fusion.
    #     :type bnd_records: list.
    #     :return: The SV line corresponding to the fusion.
    #     :rtype: str
    #     """
    #     fusion_record = self._BNDRecordsToDict(bnd_records)
    #     return super().recordToLine(fusion_record)
