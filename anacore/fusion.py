# -*- coding: utf-8 -*-
"""
Classes and functions for reading/writing fusions callers's output and VCF containing fusions. Each record is converted to VCFRecord.

:Example:

    Read results form STAR-Fusion

    .. highlight:: python
    .. code-block:: python

        # File content>
        # #FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots
        # MN1--BEND2	90	39	ONLY_REF_SPLICE	MN1^ENSG00000169184.6	chr22:27796763:-	BEND2^ENSG00000177324.14	chrX:18195442:-	YES_LDAS	38.9246	GT	1.9656	AG	1.7232	["INTERCHROMOSOMAL[chr22--chr9]"]

        from anacore.fusion import FusionFileReader

        spl_name = "spl1"
        print("Breakpoint_1\\tGene_brkpt_1\\tBreakpoint_2\\tGene_brkpt_2\\tSpanning_pair\\tSplit_read")
        with FusionFileReader.factory("test.tsv", sample_name=spl_name) as reader:
            for first, second in reader:
                print(
                    "{}:{}".format(first.chrom, frist.pos),
                    ",".join([elt["SYMBOL"] for elt in first.info["ANN"]]),
                    "{}:{}".format(second.chrom, second.pos),
                    ",".join([elt["SYMBOL"] for elt in second.info["ANN"]]),
                    first.sample["spl_name"]["PR"],
                    first.sample["spl_name"]["SR"],
                    sep="\\t"
                )

        # Result>
        # Breakpoint_1\tGene_brkpt_1\tBreakpoint_2\tGene_brkpt_2\tSpanning_pair\tSplit_read
        # chr22:27796763\tMN1\tchrX:18195442\tBEND2\t39\t99

    Read VCF

    .. highlight:: python
    .. code-block:: python

        from anacore.fusion import FusionFileReader

        # File content>
        # ##fileformat=VCFv4.3
        # ##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
        # ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
        # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
        # 2	321682	bnd_V	T	]13:123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
        # 4	888888	.	T	G	.	PASS	.
        # 13	123456	bnd_U	C	C[2:321682[,C[17:198983[	6	PASS	SVTYPE=BND;MATEID=bnd_V,bnd_Z
        # 17	198983	bnd_Z	A	]13:123456]A	6	PASS	SVTYPE=BND;MATEID=bnd_U

        print("Breakpoint_1\\tBreakpoint_1_name\\tBreakpoint_2\\tBreakpoint_2_name")
        with FusionFileReader.factory("test.vcf") as reader:
            for first, second in reader:
                print(
                    "{}:{}".format(first.chrom, frist.pos),
                    first.id
                    "{}:{}".format(second.chrom, second.pos),
                    second.id
                    sep="\\t"
                )

        # Result>
        # Breakpoint_1\tBreakpoint_1_name\tBreakpoint_2\tBreakpoint_2_name
        # 2:321682\tbnd_V\t13:123456\tbnd_U_0
        # 13:123456\tbnd_U_1\t17:198983\tbnd_Z
"""
"""
Alt representation in VCF specification:

    First in RNA and strand + => alt: N|
                    ------>
                           |
        >>>>>>>>>>>>>>     |     >>>>>>>>>>>>
        ...................|.................

                    ------>
                           |
        <<<<<<<<<<<<<<     |     <<<<<<<<<<<<
        ...................|.................

    Second in RNA and strand + => alt: |N
                            ------>
                           |
        >>>>>>>>>>>>>>     |     >>>>>>>>>>>>
        ...................|.................

                            ------>
                           |
        <<<<<<<<<<<<<<     |     <<<<<<<<<<<<
        ...................|.................

    First in RNA and strand - => alt: |N
                            <------
                           |
        >>>>>>>>>>>>>>     |     >>>>>>>>>>>>
        ...................|.................

                            <------
                           |
        <<<<<<<<<<<<<<     |     <<<<<<<<<<<<
        ...................|.................

    Second in RNA and strand - => alt: N|
                    <------
                           |
        >>>>>>>>>>>>>>     |     >>>>>>>>>>>>
        ...................|.................

                    <------
                           |
        <<<<<<<<<<<<<<     |     <<<<<<<<<<<<
        ...................|.................
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.6.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
import json
import re
import uuid
from anacore.abstractFile import isGzip
from anacore.annotVcf import AnnotVCFIO
from anacore.sv import HashedSVIO, SVIO
from anacore.vcf import decodeInfoValue, getAlleleRecord, VCFIO, VCFRecord, HeaderInfoAttr, HeaderFilterAttr, HeaderFormatAttr


def getBNDInterval(record):
    """
    Return the start and end for BND record (CIPOS is taken into account).

    :param record: The breakend record.
    :type record: anacore.vcf.VCFRecord
    :return: Start and end for BND record (CIPOS is taken into account).
    :rtype: (int, int)
    """
    start = record.pos
    end = start
    if "CIPOS" in record.info:
        start += record.info["CIPOS"][0]
        end += record.info["CIPOS"][1]
        if "IMPRECISE" in record.info:
            start -= 1
    return start, end


def getCoordDictFromCoordStr(coord):
    """
    Return elements from coordinates defined by a fusion caller ("CHROM:POS:STRAND").

    :param coord: Coordinates defined by fusion caller in format "CHROM:POS:STRAND".
    :type coord: str
    :return: Coordinates in dict format {"chrom": "CHROM", "pos": "POS", "strand":"STRAND"}.
    :rtype: dict
    """
    matches = re.match("(.+):(\d+):([+-.])", coord)
    if matches is None:
        raise Exception("The coordinates fields {} cannot be parsed.".format(coord))
    return {
        "chrom": matches.group(1),
        "pos": int(matches.group(2)),
        "strand": matches.group(3)
    }


def getAltFromCoord(first_coord, second_coord):
    """
    Return VCF alternative field value for the breakends from the 5' and 3' shard of fusion.

    :param first_coord: Coordinates of the breakend for the first 5' shard in fusion.
    :type first_coord: dict
    :param second_coord: Coordinates of the breakend for the first 3' shard in fusion.
    :type second_coord: dict
    :return: Alternative for the first shard in fusion and alternative for the second shard in fusion.
    :rtype: (str, str)
    """
    first_alt = None
    second_alt = None
    if first_coord["strand"] == "." or second_coord["strand"] == ".":
        raise Exception("getAltFromCoord cannot process coordinates with unknown strand.")
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


def getStrand(breakend, is_first=None):
    """
    Return strand from breakend record.

    :param breakend: Breakend defining a part of an RNA fusion.
    :type breakend: anacore.vcf.VCFRecord
    :param is_first: The shard is the first in RNA.
    :type is_first: boolean
    :return: Strand from breakend record.
    :rtype: str
    """
    if is_first is None:
        is_first = True if "RNA_FIRST" in breakend.info else False
    strand = None
    if breakend.alt[0].startswith("[") or breakend.alt[0].startswith("]"):
        strand = "-" if is_first else "+"
    else:
        strand = "+" if is_first else "-"
    return strand


def getCoordStr(breakend, is_first=None):
    """
    Return coordinates elements (chrom, pos and strand) a from breakend record.

    :param breakend: Breakend defining a part of an RNA fusion.
    :type breakend: anacore.vcf.VCFRecord
    :param is_first: The shard is the first in RNA.
    :type is_first: boolean
    :return: Coordinates in dict format {"chrom": "CHROM", "pos": "POS", "strand":"STRAND"}.
    :rtype: dict
    """
    return {"chrom": breakend.chrom, "pos": breakend.pos, "strand": getStrand(breakend, is_first)}


class FusionFileReader(object):
    """Factory to identify and return the dedicated to the software used to produce the fusions file."""

    @staticmethod
    def factory(filepath, *args, **kwargs):
        """
        Return the parser dedicated to the software used to produce the fusions file.

        :param filepath: Path to the file (format: TSV).
        :type filepath: str
        :param args: Additional arguments.
        :type args: list
        :param kwargs: Additional keyword arguments.
        :type kwargs: dict
        :return: The parser dedicated to the software used to produce the fusions file.
        :rtype: ArribaIO or FusionCatcherIO or STARFusionIO
        """
        if STARFusionIO.isValid(filepath):
            return STARFusionIO(filepath, *args, **kwargs)
        elif FusionCatcherIO.isValid(filepath):
            return FusionCatcherIO(filepath, *args, **kwargs)
        elif ArribaIO.isValid(filepath):
            return ArribaIO(filepath, *args, **kwargs)
        elif BreakendVCFIO.isValid(filepath):
            return BreakendVCFIO(filepath, *args, **kwargs)
        else:
            raise IOError("The file {} does not have a valid format for 'FusionFileReader'.".format(filepath))


class FusionCatcherIO(HashedSVIO):
    """Class to read and write fusions in fusionCatcher TSV output format (see https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md#output-data). The records are in anacore.vcf.VCFRecord format and a fusion is represented by 2 records (breakend for the 5' shard and breakend for the 3' shard)."""

    titles = [
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

    def __init__(self, filepath, mode="r", annot_field="FCANN", sample_name="sample"):
        """
        Build and return an instance of FusionCatcherIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param annot_field: The tag for the field used to store annotations.
        :type annot_field: str
        :param sample_name: The name of the sample.
        :type sample_name: str
        :return: The new instance.
        :rtype: FusionCatcherIO
        """
        super().__init__(filepath, mode, "\t")
        self.annot_field = annot_field
        self.sample_name = sample_name
        if self.titles is None:
            self.titles = FusionCatcherIO.titles

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
        filters_field = [] if fusion_record["Fusion_description"] == "" else sorted([elt.strip() for elt in fusion_record["Fusion_description"].split(",")])
        if fusion_record["Spanning_unique_reads"] == "0":
            filters_field.append("Imprecise")
            filters_field.sort()
        sample_field = {
            self.sample_name: {
                "PR": int(fusion_record["Spanning_pairs"]),
                "SR": int(fusion_record["Spanning_unique_reads"]),
                "CM": int(fusion_record["Counts_of_common_mapping_reads"]),
                "LA": int(fusion_record["Longest_anchor_found"])
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
            pFilter=filters_field,
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
            pFilter=filters_field,
            info=second_info,
            pFormat=format_field,
            samples=sample_field
        )
        return [first_bnd, second_bnd]

    def write(self, breakends):
        """
        Write record line in file.

        :param record: The first breakend and the second breakend in pair.
        :type record: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        """
        super().write(self._BNDRecordsToDict(breakends))

    def _BNDRecordsToDict(self, breakends):
        """
        Return the FusionCatcher dict from the two breakends of a fusion.

        :param breakends: The two breakends of a fusion.
        :type breakends: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        :return: The FusionCatcher dict.
        :rtype: dict
        """
        first, second = breakends
        first_coord = getCoordStr(first, True)
        second_coord = getCoordStr(second, False)
        return {
            "Gene_1_symbol(5end_fusion_partner)": first.info[self.annot_field][0]["SYMBOL"],
            "Gene_2_symbol(3end_fusion_partner)": second.info[self.annot_field][0]["SYMBOL"],
            "Fusion_description": ",".join([elt for elt in first.filter if elt != "Imprecise"]),
            "Counts_of_common_mapping_reads": first.samples[self.sample_name]["CM"],
            "Spanning_pairs": first.samples[self.sample_name]["PR"],
            "Spanning_unique_reads": first.samples[self.sample_name]["SR"],
            "Longest_anchor_found": first.samples[self.sample_name]["LA"],
            "Fusion_finding_method": ";".join(first.info["SOURCES"]),
            "Fusion_point_for_gene_1(5end_fusion_partner)": "{}:{}:{}".format(
                first_coord["chrom"],
                first_coord["pos"],
                first_coord["strand"]
            ),
            "Fusion_point_for_gene_2(3end_fusion_partner)": "{}:{}:{}".format(
                second_coord["chrom"],
                second_coord["pos"],
                second_coord["strand"]
            ),
            "Gene_1_id(5end_fusion_partner)": first.info[self.annot_field][0]["Gene"],
            "Gene_2_id(3end_fusion_partner)": second.info[self.annot_field][0]["Gene"],
            "Exon_1_id(5end_fusion_partner)": first.info[self.annot_field][0]["EXON"],
            "Exon_2_id(3end_fusion_partner)": second.info[self.annot_field][0]["EXON"],
            "Fusion_sequence": second.info["RNA_CONTIG"],
            "Predicted_effect": second.info[self.annot_field][0]["Effect"]
        }

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in FusionCatcher output format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in FusionCatcher output format.
        :rtype: bool
        """
        is_valid = False
        try:
            with SVIO(filepath) as reader:
                mandatory_fields = set(FusionCatcherIO.titles) - {"Predicted_effect"}
                if len(mandatory_fields - set(reader.titles)) == 0:
                    is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid

    @staticmethod
    def setVCFHeader(vcf_io, annotation_field="FCANN"):
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


starfusion_filters = {
    "GTEx_recurrent_StarF2019": HeaderFilterAttr(
        "GTEx_recurrent_StarF2019",
        "Fusions found recurrently in GTEx as mined using STAR-Fusion v1.5.0"
    ),
    "BodyMap": HeaderFilterAttr(
        "BodyMap",
        "Fusions found by STAR-Fusion as applied to the Illumina Human Body Map reference data"
    ),
    "DGD_PARALOGS": HeaderFilterAttr(
        "DGD_PARALOGS",
        "Duplicated genes as per the Duplicated Genes Database"
    ),
    "HGNC_GENEFAM": HeaderFilterAttr(
        "HGNC_GENEFAM",
        "HGNC gene family membership as per ftp://ftp.ebi.ac.uk/pub/databases/genenames/genefam_list.txt.gz"
    ),
    "Greger_Normal": HeaderFilterAttr(
        "Greger_Normal",
        "Fusion transcripts (mostly from tandem genes) detected based on analysis of RNA-Seq from 1000 genomes project samples. List derived from Greger et al. PLOS One, 2014"
    ),
    "Babiceanu_Normal": HeaderFilterAttr(
        "Babiceanu_Normal",
        "Recurrent chimeric fusion RNAs in non-cancer tissues and cells as per Babiceanu et al. NAR, 2016"
    ),
    "ConjoinG": HeaderFilterAttr(
        "ConjoinG",
        "Fused transcripts derived from the Conjoined Gene Database"
    ),
    "Imprecise": HeaderFilterAttr(
        "Imprecise",
        "RNA fusion candidates for which no spanning contig was found."
    )
}


class STARFusionIO(HashedSVIO):
    """Class to read and write fusions in STAR-Fusion TSV output format (see https://github.com/STAR-Fusion/STAR-Fusion/wiki#Outputs). The records are in anacore.vcf.VCFRecord format and a fusion is represented by 2 records (breakend for the 5' shard and breakend for the 3' shard)."""

    titles = [
        "FusionName",
        "JunctionReadCount",
        "SpanningFragCount",
        "SpliceType",
        "LeftGene",
        "LeftBreakpoint",
        "RightGene",
        "RightBreakpoint",
        "JunctionReads",
        "SpanningFrags",
        "LargeAnchorSupport",
        "FFPM",
        "LeftBreakDinuc",
        "LeftBreakEntropy",
        "RightBreakDinuc",
        "RightBreakEntropy",
        "annots"
    ]

    def __init__(self, filepath, mode="r", annot_field="FCANN", sample_name="sample"):
        """
        Build and return an instance of STARFusionIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param annot_field: The tag for the field used to store annotations.
        :type annot_field: str
        :param sample_name: The name of the sample.
        :type sample_name: str
        :return: The new instance.
        :rtype: STARFusionIO
        """
        super().__init__(filepath, mode=mode, separator="\t", title_starter="#")
        self.annot_field = annot_field
        self.sample_name = sample_name
        if self.titles is None:
            self.titles = STARFusionIO.titles

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
        breakends = []
        large_anchor_support = "1" if fusion_record["LargeAnchorSupport"] == "YES_LDAS" else "0"
        info_annot_tags = [elt.strip() for elt in json.loads(fusion_record["annots"])]
        filters_field = sorted(list(set(starfusion_filters) & set(info_annot_tags)))
        if fusion_record["JunctionReadCount"] == "0":
            filters_field.append("Imprecise")
            filters_field.sort()
        if len(filters_field) == 0:
            filters_field = ["PASS"]
        coord_list = []
        for side, mate_side in [("Left", "Right"), ("Right", "Left")]:
            id = str(uuid.uuid4())
            coord = getCoordDictFromCoordStr(fusion_record[side + "Breakpoint"])
            coord_list.append(coord)
            # INFO
            gene_symbol, gene_id = fusion_record[side + "Gene"].split("^")
            info = {
                "SVTYPE": "BND",
                "BREAK_DINUC": fusion_record[side + "BreakDinuc"],
                "BREAK_ENTROPY": fusion_record[side + "BreakEntropy"],
                "SPLICE_TYPE": fusion_record["SpliceType"],
                self.annot_field: [{
                    "SYMBOL": gene_symbol,
                    "Gene": gene_id,
                    "Tags": "&".join(info_annot_tags)
                }]
            }
            # SAMPLES
            samples_field = {
                self.sample_name: {
                    "PR": int(fusion_record["SpanningFragCount"]),
                    "SR": int(fusion_record["JunctionReadCount"]),
                    "hasLAS": large_anchor_support,
                    "FFPM": float(fusion_record["FFPM"])
                }
            }
            if side == "Left":
                info["RNA_FIRST"] = True
                if "JunctionReads" in fusion_record:
                    samples_field[self.sample_name]["JRL"] = fusion_record["JunctionReads"]
                if "SpanningFrags" in fusion_record:
                    samples_field[self.sample_name]["SFL"] = fusion_record["SpanningFrags"]
            # Record
            breakends.append(
                VCFRecord(
                    coord["chrom"],
                    coord["pos"],
                    id,
                    "N",  # ref
                    [],  # alt
                    info=info,
                    pFilter=filters_field,
                    pFormat=sorted(samples_field[self.sample_name].keys()),
                    samples=samples_field
                )
            )
        left_alt, right_alt = getAltFromCoord(*coord_list)
        breakends[0].alt = [left_alt]
        breakends[1].alt = [right_alt]
        breakends[0].info["MATEID"] = [breakends[1].id]
        breakends[1].info["MATEID"] = [breakends[0].id]
        return breakends

    def write(self, breakends):
        """
        Write record line in file.

        :param record: The first breakend and the second breakend in pair.
        :type record: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        """
        super().write(self._BNDRecordsToDict(breakends))

    def _BNDRecordsToDict(self, breakends):
        """
        Return the STAR-fusion dict from the two breakends of a fusion.

        :param breakends: The two breakends of a fusion.
        :type breakends: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        :return: The STAR-fusion dict.
        :rtype: dict
        """
        first, second = breakends
        first_coord = getCoordStr(first, True)
        second_coord = getCoordStr(second, False)
        return {
            "FusionName": "{}--{}".format(
                first.info[self.annot_field][0]["SYMBOL"],
                second.info[self.annot_field][0]["SYMBOL"],
            ),
            "JunctionReadCount": first.samples[self.sample_name]["SR"],
            "SpanningFragCount": first.samples[self.sample_name]["PR"],
            "SpliceType": first.info["SPLICE_TYPE"],
            "LeftGene": "{}^{}".format(
                first.info[self.annot_field][0]["SYMBOL"],
                first.info[self.annot_field][0]["Gene"]
            ),
            "LeftBreakpoint": "{}:{}:{}".format(
                first_coord["chrom"],
                first_coord["pos"],
                first_coord["strand"]
            ),
            "RightGene": "{}^{}".format(
                second.info[self.annot_field][0]["SYMBOL"],
                second.info[self.annot_field][0]["Gene"]
            ),
            "RightBreakpoint": "{}:{}:{}".format(
                second_coord["chrom"],
                second_coord["pos"],
                second_coord["strand"]
            ),
            "JunctionReads": "" if "JRL" not in first.samples[self.sample_name] else first.samples[self.sample_name]["JRL"],
            "SpanningFrags": "" if "SFL" not in first.samples[self.sample_name] else first.samples[self.sample_name]["SFL"],
            "LargeAnchorSupport": "YES_LDAS" if first.samples[self.sample_name]["hasLAS"] == "1" else "NO_LDAS",
            "FFPM": first.samples[self.sample_name]["FFPM"],
            "LeftBreakDinuc": first.info["BREAK_DINUC"],
            "LeftBreakEntropy": first.info["BREAK_ENTROPY"],
            "RightBreakDinuc": second.info["BREAK_DINUC"],
            "RightBreakEntropy": second.info["BREAK_ENTROPY"],
            "annots": json.dumps(first.info[self.annot_field][0]["Tags"].split("&"))
        }

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in STAR-Fusion output format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in STAR-Fusion output format.
        :rtype: bool
        """
        is_valid = False
        try:
            with SVIO(filepath, title_starter="#") as reader:
                mandatory_fields = set(STARFusionIO.titles) - {"JunctionReads", "SpanningFrags"}
                if len(mandatory_fields - set(reader.titles)) == 0:  # All mandatory fields of Arriba are in reader
                    is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid

    @staticmethod
    def setVCFHeader(vcf_io, annotation_field="FCANN"):
        """
        Set header of an AnnotVCFIO to write STAR-Fusion records in VCF output. This function pre-set ANN_titles, info and format used by STARFusionIO records.

        :param vcf_io: The VCF used to write STAR-Fusion records.
        :type vcf_io: anacore.annotVcf.AnnotVCFIO
        :param annotation_field: INFO field used for store annotations.
        :type annotation_field: str
        """
        # FILTERS
        vcf_io.filter = starfusion_filters
        # INFO
        vcf_io.info = {
            "MATEID": HeaderInfoAttr("MATEID", type="String", number="A", description="ID of mate breakend."),
            "SVTYPE": HeaderInfoAttr("SVTYPE", type="String", number="1", description="Type of structural variant."),
            "RNA_FIRST": HeaderInfoAttr("RNA_FIRST", type="Flag", number="0", description="For RNA fusions, this breakend is 5' in the fusion transcript."),  # In the fusion transcript, the 'LeftGene' is always 5' relative to the 'RightGene' (https://groups.google.com/forum/#!topic/star-fusion/9rBLT-v5JHI)
            "BREAK_DINUC": HeaderInfoAttr("BREAK_DINUC", type="String", number="1", description="Dinucleotides flanking the breakend of the fragment excluded by the fusion."),
            "BREAK_ENTROPY": HeaderInfoAttr("BREAK_ENTROPY", type="Float", number="1", description="Shannon entropy of the 15 exonic bases flanking the breakpoint. The maximum entropy is 2, representing highest complexity. The lowest would be zero (involving a 15 base mononucleotide run). Low entropy sites should generally be treated as less confident breakpoints."),
            "SPLICE_TYPE": HeaderInfoAttr("SPLICE_TYPE", type="String", number="1", description="Whether the proposed breakpoint occurs at reference exon junctions as provided by the reference transcript structure annotations (ex. gencode)."),
            annotation_field: HeaderInfoAttr(annotation_field, type="String", number=".", description="Annotation generated by FusionAnnotator (see: https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki). Format: SYMBOL|Gene|Tags")
        }
        # ANN_titles
        vcf_io.ANN_titles = ["SYMBOL", "Gene", "Tags"]
        # FORMAT
        vcf_io.format = {
            "PR": HeaderFormatAttr("PR", type="Integer", number="1", description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (SpanningFragCount)."),
            "SR": HeaderFormatAttr("SR", type="Integer", number="1", description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction (JunctionReadCount)."),
            "hasLAS": HeaderFormatAttr("hasLAS", type="String", number="1", description="This column indicates whether there are split reads that provide 'long' (set to length of 25 bases) alignments on both sides of the putative breakpoint (LargeAnchorSupport)."),
            "FFPM": HeaderFormatAttr("FFPM", type="Float", number="1", description="Normalized measures of the fusion-supporting rna-seq fragments (fusion fragments per million total reads)."),
            "JRL": HeaderFormatAttr("JRL", type="String", number="1", description="RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction."),
            "SFL": HeaderFormatAttr("SFL", type="String", number="1", description="RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.")
        }


arriba_event_level_filters = {
    "merge_adjacent": HeaderFilterAttr(
        "merge_adjacent",
        "When multiple alternative alignments with equal quality are possible, STAR does not always choose the same one. This leads to multiple breakpoints in close proximity, all of which arise from one and the same event. This filter merges all breakpoints in a window of 5 bp to the breakpoint with the highest number of supporting reads."
    ),
    "non_coding_neighbors": HeaderFilterAttr(
        "non_coding_neighbors",
        "Genes which are not well studied suffer from incomplete annotation. Many exons are annotated as separate genes even though they might actually be part of one and the same gene. Predicted genes named RP11-... are common examples for this. When poorly understood genes lie next to each other on the same strand, this would frequently lead to false positive predictions of deletions, because the transcripts that span both genes give rise to reads, which resemble focal deletions. This filter discards deletions, which are predicted between two neighboring genes, if both genes are non-coding or one breakpoint is intergenic."
    ),
    "intragenic_exonic": HeaderFilterAttr(
        "intragenic_exonic",
        "Since exons usually make up only a small fraction of a gene, it is more likely that a genomic rearrangement starts and ends in intronic regions. On the transcriptomic level, this manifests as breakpoints at splice-sites or in introns. Many candidates found by STAR have both breakpoints within exons of the same gene. This is particularly true for intragenic events, which are prone to in vitro artifacts. This filter removes intragenic events, if both breakpoints are in exons and more than 80% of the region between the breakpoints is intronic, such that it should be very unlikely that both breakpoints are located inside exons."
    ),
    "min_support": HeaderFilterAttr(
        "min_support",
        "This filter discards all events with few reads."
    ),
    "relative_support": HeaderFilterAttr(
        "relative_support",
        "With increasing expression of a gene the number of chimeric reads mapping to the gene increases, too. Arriba assumes a polynomial relationship between the number of events and the number of supporting reads for a given event. This assumption is based on empirical evidence. An expected value (e-value) is calculated for every event, reflecting how many events with the given number of supporting reads are expected by chance. This filter selects only those events with a low e-value, i.e., that have a high number of supporting reads relative to the overall number of events in a gene. Multiple covariates are taken into account, such as whether a breakpoint is at a splice-site or whether the event is intragenic. The parameters of the polynomial relationship are fixed and were estimated from several hundred RNA-Seq samples of various cancer types."
    ),
    "intronic": HeaderFilterAttr(
        "intronic",
        "This filter discards an event, when none of the supporting reads overlaps an exon (i.e., all alignments are in introns or intergenic regions). True events most often have at least one breakpoint in an exon or at an exon boundary."
    ),
    "known_fusions": HeaderFilterAttr(
        "known_fusions",
        "Events with too few supporting reads but known in cancer."
    ),  # Protection
    "pcr_fusions": HeaderFilterAttr(
        "pcr_fusions",
        "In some tissues certain genes are expressed at very high levels, for example hemoglobin and fibrinogen in blood or collagens in connective tissue. Presumably, the abundance of fragments from such genes increases the chance of unrelated molecules sticking together, which can serve as undesired primers for PCR or may cause the reverse transcriptase enzyme to switch templates. These processes generate a large amount of chimeric fragments in vitro. Such artifactual fusions can be recognized as an extraordinary number of events with breakpoints within exons (rather than at exon boundaries, which is more common for true predictions). This filter eliminates events with genes that are highly expressed (top 0.2%) and have an unbalanced number of split-reads vs. discordant mates or that have an excessive amount of intra-exonic breakpoints."
    ),
    "spliced": HeaderFilterAttr(
        "spliced",
        "This filter recovers events discarded due to a low number of supporting reads, given that both breakpoints of the event are at splice-sites and there is at least one additional event linking the same pair of genes."
    ),  # Protection
    "select_best": HeaderFilterAttr(
        "select_best",
        "If there are multiple breakpoints detected between the same pair of genes, this filter discards all but the most credible one. Events with split reads in both genes are preferred over events with only discordant mates, because in the latter case, the precise breakpoint is unknown. Moreover, events with a higher number of supporting reads are favored."
    ),
    "many_spliced": HeaderFilterAttr(
        "many_spliced",
        "Occassionally, a genomic rearrangement produces multiple alternatively spliced transcripts, all of which are low expressed and therefore discarded by the filters relative_support or min_support. This filter recovers events that were discarded due to too few supporting reads, under the circumstances that there are many events between a pair of genes with at least one of the breakpoints at splice-sites. All events between the pair of genes with at least one spliced breakpoint are recovered."
    ),  # Protection
    "no_genomic_support": HeaderFilterAttr(
        "no_genomic_support",
        "This filter removes events with low confidence, if they are not confirmed by structural variant calls obtained from whole-genome sequencing."
    ),
    "blacklist": HeaderFilterAttr(
        "blacklist",
        "This filter removes events with breakpoint coordinates matching entries in the blacklist. If an event has no split-reads (but only discordant mates), the precise breakpoint coordinates are unknown. In this case the event is discarded, if the breakpoints are within a range of the insert size around the blacklisted coordinates."
    ),
    "short_anchor": HeaderFilterAttr(
        "short_anchor",
        "A chimeric read aligns to some part in one of the fused genes and to some part to the other gene. The anchor of a read is the longer aligned segment. It is more likely to be aligned correctly. True positives often have anchors in both fused genes, whereas alignment artifacts are frequently characterized by only a small segment aligning to one of the genes and all anchors to the other."
    ),
    "end_to_end": HeaderFilterAttr(
        "end_to_end",
        "Theoretically, it should be impossible to observe fusions which only retain the 3' ends of the fusion partners, because the promoter would be missing. In reality, end-to-end fused transcripts are observed, albeit rarely. Most of them are false positives. Therefore, this filter discards such events, unless they have a lot of supporting reads, i.e., split-reads in both fusion partners or split-reads and discordant mates."
    ),
    "no_coverage": HeaderFilterAttr(
        "no_coverage",
        "For intronic and intragenic breakpoints as well as read-through fusions, this filter checks, if there is some coverage in the vicinity of the breakpoint in the normal alignments. Only reads that are not chimeric are considered. When there are no non-chimeric reads near the breakpoint, this is indicative of an alignment artifact."
    ),
    "homologs": HeaderFilterAttr(
        "homologs",
        "This filter discards events between genes that have high sequence homology, which frequently leads to erroneous alignments. Homology is quantified by counting the number of shared 16-mers."
    ),
    "mismappers": HeaderFilterAttr(
        "mismappers",
        "Many alignment artifacts are caused by an excessive number of reference mismatches in close proximity due to sequencing errors or adjacent SNPs. STAR tends to clip reads when it encounters a cluster of mismatches. The clipped segment is then used for a chimeric alignment, which occassionally aligns elsewhere in the genome. The filter mismappers performs a sensitive realignment of both segments. If both segments can be aligned to the same gene (while allowing more mismatches than STAR does), the chimeric alignment is considered to be an artifact. When an important proportion of the supporting reads are classified as being aligned incorrectly, the event is discarded."
    ),
    "genomic_support": HeaderFilterAttr(
        "genomic_support",
        "This filter recovers events which were discarded by previous filters due to few supporting reads, but which can be explained by genomic rearrangements as evidenced by structural variant calls obtained from whole-genome sequencing data. Arriba considers structural variant calls to match with breakpoints seen in transcriptomic data, when the breakpoints are less then 100 kb apart and the orientation of the genomic and transcriptomic breakpoints are identical."
    ),
    "isoforms": HeaderFilterAttr(
        "isoforms",
        "This filter searches for additional isoforms for those gene pairs that are predicted to be fused in proper orientation. Typically there is a major isoform which is expressed at a high level and a few low expressed isoforms."
    ),
    "Imprecise": HeaderFilterAttr(
        "Imprecise",
        "RNA fusion candidates for which no spanning contig was found."
    ),
    "Unstranded": HeaderFilterAttr(
        "Unstranded",
        "RNA fusion is not stranded."
    )
}


class ArribaIO(HashedSVIO):
    """Class to read and write fusions in Arriba TSV output format (see https://github.com/suhrig/arriba/blob/master/documentation/output-files.md). The records are in anacore.vcf.VCFRecord format and a fusion is represented by 2 records (breakend for the 5' shard and breakend for the 3' shard)."""

    titles = [
        "gene1",
        "gene2",
        "strand1(gene/fusion)",
        "strand2(gene/fusion)",
        "breakpoint1",
        "breakpoint2",
        "site1",
        "site2",
        "type",
        "direction1",
        "direction2",
        "split_reads1",
        "split_reads2",
        "discordant_mates",
        "coverage1",
        "coverage2",
        "confidence",
        "closest_genomic_breakpoint1",
        "closest_genomic_breakpoint2",
        "filters",
        "fusion_transcript",
        "reading_frame",
        "peptide_sequence",
        "read_identifiers"
    ]

    def __init__(self, filepath, mode="r", annot_field="FCANN", sample_name="sample"):
        """
        Build and return an instance of ArribaIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param annot_field: The tag for the field used to store annotations.
        :type annot_field: str
        :param sample_name: The name of the sample.
        :type sample_name: str
        :return: The new instance.
        :rtype: ArribaIO
        """
        super().__init__(filepath, mode=mode, separator="\t", title_starter="#")
        self.annot_field = annot_field
        self.sample_name = sample_name
        if self.titles is None:
            self.titles = ArribaIO.titles

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
        is_unstranded = fusion_record["strand1(gene/fusion)"][-1] == "." or fusion_record["strand2(gene/fusion)"][-1] == "."
        event_level_filters = set(arriba_event_level_filters)
        filters_from_file = {elt.strip() for elt in fusion_record["filters"].split(",")}
        filters_field = sorted(list(filters_from_file & event_level_filters))
        read_level_filters = sorted(list(filters_from_file - event_level_filters))
        is_imprecise = False if int(fusion_record["split_reads1"]) + int(fusion_record["split_reads2"]) > 0 else True
        if is_imprecise:
            filters_field.append("Imprecise")
        if is_unstranded:
            filters_field.append("Unstranded")
        if is_imprecise or is_unstranded:
            filters_field.sort()
        breakends = []
        coord_list = []
        for side, mate_side in [("1", "2"), ("2", "1")]:
            id = str(uuid.uuid4())
            strand = fusion_record["strand" + side + "(gene/fusion)"].split("/")[1]
            if is_unstranded:
                strand = "+"
                if side == "1":
                    if fusion_record["direction1"] == "downstream":  # The first in fusion go to the breakpoint: ---->| or |<----
                        strand = "-"
                else:
                    if fusion_record["direction2"] == "upstream":  # The second in fusion go from the breakpoint: <----| or |---->
                        strand = "-"
            curr_breakpoint = fusion_record["breakpoint" + side] + ":" + strand
            coord = getCoordDictFromCoordStr(curr_breakpoint)
            coord_list.append(coord)
            # INFO
            has_frameshift = {"in-frame": "false", "out-of-frame": "true", ".": "."}
            info = {
                "SVTYPE": "BND",
                self.annot_field: [{
                    "SYMBOL": fusion_record["gene" + side],
                    "STRAND": fusion_record["strand" + side + "(gene/fusion)"].split("/")[0],
                    "Site": fusion_record["site" + side],
                    "Type": fusion_record["type"],
                    "GENE_SHARD": fusion_record["direction" + side].replace("stream", ""),
                    "FRAMESHIFT": has_frameshift[fusion_record["reading_frame"]],
                    "Protein_contig": fusion_record["peptide_sequence"].replace("|", "@")
                }]
            }
            if is_imprecise:
                info["IMPRECISE"] = True
            if is_unstranded:
                info["UNSTRANDED"] = True
            if fusion_record["closest_genomic_breakpoint" + side] != ".":
                info["GBP"] = fusion_record["closest_genomic_breakpoint" + side]
            # SAMPLES
            samples_field = {
                self.sample_name: {
                    "PR": int(fusion_record["discordant_mates"]),
                    "SR": int(fusion_record["split_reads1"]) + int(fusion_record["split_reads2"]),
                    "SR1": int(fusion_record["split_reads1"]),
                    "SR2": int(fusion_record["split_reads2"]),
                    "CFD": fusion_record["confidence"],
                    "RFIL": read_level_filters
                }
            }
            if fusion_record["coverage" + side] != ".":
                samples_field[self.sample_name]["DPS"] = int(fusion_record["coverage" + side])
            if fusion_record["read_identifiers"] != ".":
                samples_field[self.sample_name]["SRL"] = fusion_record["read_identifiers"].split(",")
            if side == "1":
                info["RNA_FIRST"] = True
            else:
                if fusion_record["fusion_transcript"] != ".":
                    info["RNA_CONTIG"] = fusion_record["fusion_transcript"].replace("|", "@")
            # Record
            breakends.append(
                VCFRecord(
                    coord["chrom"],
                    coord["pos"],
                    id,
                    "N",  # ref
                    [],  # alt
                    pFilter=["PASS"] if len(filters_field) == 0 else filters_field,
                    info=info,
                    pFormat=sorted(samples_field[self.sample_name].keys()),
                    samples=samples_field
                )
            )
        left_alt, right_alt = getAltFromCoord(*coord_list)
        breakends[0].alt = [left_alt]
        breakends[1].alt = [right_alt]
        breakends[0].info["MATEID"] = [breakends[1].id]
        breakends[1].info["MATEID"] = [breakends[0].id]
        return breakends

    def write(self, breakends):
        """
        Write record line in file.

        :param record: The first breakend and the second breakend in pair.
        :type record: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        """
        super().write(self._BNDRecordsToDict(breakends))

    def _BNDRecordsToDict(self, breakends):
        """
        Return the Arriba dict from the two breakends of a fusion.

        :param breakends: The two breakends of a fusion.
        :type breakends: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
        :return: The Arriba dict.
        :rtype: dict
        """
        first, second = breakends
        first_coord = getCoordStr(first, True)
        second_coord = getCoordStr(second, False)
        fs_to_reading_frame = {".": ".", "false": "in-frame", "true": "out-of-frame"}
        arriba_filters = ""
        if first.filter:
            excluded = {"Imprecise", "Unstranded"}
            arriba_filters = ",".join(sorted(set(first.filter) - excluded))
        if "RFIL" in first.samples[self.sample_name] and first.samples[self.sample_name]["RFIL"]:
            arriba_filters += ",".join(first.samples[self.sample_name]["RFIL"])
        return {
            "gene1": first.info[self.annot_field][0]["SYMBOL"],
            "gene2": second.info[self.annot_field][0]["SYMBOL"],
            "strand1(gene/fusion)": "{}/{}".format(
                first.info[self.annot_field][0]["STRAND"],
                "." if "UNSTRANDED" in first.info else first_coord["strand"]
            ),
            "strand2(gene/fusion)": "{}/{}".format(
                second.info[self.annot_field][0]["STRAND"],
                "." if "UNSTRANDED" in first.info else second_coord["strand"]
            ),
            "breakpoint1": "{}:{}".format(first_coord["chrom"], first_coord["chrom"]),
            "breakpoint2": "{}:{}".format(second_coord["chrom"], second_coord["chrom"]),
            "site1": first.info[self.annot_field][0]["Site"],
            "site2": second.info[self.annot_field][0]["Site"],
            "type": second.info[self.annot_field][0]["Type"],
            "direction1": first.info[self.annot_field][0]["GENE_SHARD"] + "stream",
            "direction2": second.info[self.annot_field][0]["GENE_SHARD"] + "stream",
            "split_reads1": first.samples[self.sample_name]["SR1"],
            "split_reads2": first.samples[self.sample_name]["SR2"],
            "discordant_mates": first.samples[self.sample_name]["PR"],
            "coverage1": first.samples[self.sample_name]["DPS"],
            "coverage2": second.samples[self.sample_name]["DPS"],
            "confidence": first.samples[self.sample_name]["CFD"],
            "closest_genomic_breakpoint1": "." if "GBP" not in first.info else first.info["GBP"],
            "closest_genomic_breakpoint2": "." if "GBP" not in second.info else second.info["GBP"],
            "filters": arriba_filters,
            "fusion_transcript": "." if "RNA_CONTIG" not in first.info else first.info["RNA_CONTIG"].replace("@", "|"),
            "reading_frame": fs_to_reading_frame[second.info[self.annot_field][0]["FRAMESHIFT"]],
            "peptide_sequence": "." if first.info[self.annot_field][0]["Protein_contig"] is None else first.info[self.annot_field][0]["Protein_contig"].replace("@", "|"),
            "read_identifiers": "." if "SRL" not in first.samples[self.sample_name] or len(first.samples[self.sample_name]["SRL"]) == 0 else first.samples[self.sample_name]["SRL"],
        }

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in Arriba output format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in Arriba output format.
        :rtype: bool
        """
        is_valid = False
        try:
            with SVIO(filepath, title_starter="#") as reader:
                if len(set(ArribaIO.titles) - set(reader.titles)) == 0:  # All mandatory fields of Arriba are in reader
                    is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid

    @staticmethod
    def setVCFHeader(vcf_io, annotation_field="FCANN"):
        """
        Set header of an AnnotVCFIO to write Arriba records in VCF output. This function pre-set ANN_titles, info and format used by ArribaIO records.

        :param vcf_io: The VCF used to write Arriba records.
        :type vcf_io: anacore.annotVcf.AnnotVCFIO
        :param annotation_field: INFO field used for store annotations.
        :type annotation_field: str
        """
        # INFO
        vcf_io.info = {
            "GBP": HeaderInfoAttr("GBP", type="String", number="1", description="The coordinates of the genomic breakpoint which is closest to the transcriptomic breakpoint."),
            "IMPRECISE": HeaderInfoAttr("IMPRECISE", type="Flag", number="0", description="Imprecise structural variation."),
            "MATEID": HeaderInfoAttr("MATEID", type="String", number="A", description="ID of mate breakend."),
            "RNA_CONTIG": HeaderInfoAttr("RNA_CONTIG", type="String", number="1", description="The transcript sequence assembled from the supporting reads of the most highly expressed transcript."),
            "RNA_FIRST": HeaderInfoAttr("RNA_FIRST", type="Flag", number="0", description="For RNA fusions, this break-end is 5' in the fusion transcript."),  # In the fusion transcript, the 'LeftGene' is always 5' relative to the 'RightGene' (https://groups.google.com/forum/#!topic/star-fusion/9rBLT-v5JHI)
            "SVTYPE": HeaderInfoAttr("SVTYPE", type="String", number="1", description="Type of structural variant."),
            "UNSTRANDED": HeaderInfoAttr("UNSTRANDED", type="Flag", number="0", description="Unstranded structural variation."),
            annotation_field: HeaderInfoAttr(annotation_field, type="String", number=".", description="Consequence annotations. Format: SYMBOL|STRAND|Site|Type|GENE_SHARD|FRAMESHIFT|Protein_contig")
        }
        # ANN_titles
        vcf_io.ANN_titles = ["SYMBOL", "STRAND", "Site", "Type", "GENE_SHARD", "FRAMESHIFT", "Protein_contig"]
        # FILTERS
        vcf_io.filter = arriba_event_level_filters
        # FORMAT
        vcf_io.format = {
            "PR": HeaderFormatAttr("PR", type="Integer", number="1", description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment."),
            "SR": HeaderFormatAttr("SR", type="Integer", number="1", description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction."),
            "SR1": HeaderFormatAttr("SR1", type="Integer", number="1", description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on first shard."),
            "SR2": HeaderFormatAttr("SR2", type="Integer", number="1", description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on second shard."),
            "CFD": HeaderFormatAttr("CFD", type="String", number="1", description="Each prediction is assigned one of the confidences low, medium, or high. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations."),
            "DPS": HeaderFormatAttr("DPS", type="Integer", number="1", description="Coverage near breakpoint. The coverage is calculated as the number of fragments near the breakpoint on the side of the breakpoint that is retained in the fusion transcript. Note that the coverage calculation counts all fragments (even duplicates)."),
            "RFIL": HeaderFormatAttr("RFIL", type="String", number=".", description="Filters which removed one or more of the supporting reads. The number of filtered reads is given in parantheses after the name of the filter. If a filter discarded the event as a whole (all reads), the number of filtered reads is missing."),
            "SRL": HeaderFormatAttr("SRL", type="String", number=".", description="The names of the supporting reads.")
        }


class BreakendVCFIO(AnnotVCFIO):
    """Read and write VCF file containing breakends. Each iteration return a couple of breakends (the first and the second in fusion)."""

    def __init__(self, filepath, mode="r", annot_field=None):
        """
        Return instance of BreakendVCFIO.

        :param filepath: The filepath.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a', 'i'). The mode 'i' allow to open file in indexed mode to fetch by region (tabix).
        :type mode: str
        :return: The new instance.
        :rtype: anacore.fusionVcf.BreakendVCFIO
        """
        self.sv_rec_by_id = {}  # By record ID the byte position of start and end of line
        self.pick_reader = None
        if mode in {"r", "i"}:
            if isGzip(filepath):
                self.pick_reader = gzip.open(filepath, "rt")
            else:
                self.pick_reader = open(filepath, "r")
            self.loadIndex()
        self._iter_already_processed = set()
        super().__init__(filepath, mode, annot_field)

    def loadIndex(self):
        """Parse file and store in index the byte positions (start and end) of each record by record ID."""
        self.sv_rec_by_id = {}
        next_start = 0
        line_nb = 1
        line = self.pick_reader.readline()
        while line:
            start = next_start
            next_start = self.pick_reader.tell()
            if not line.startswith("#"):
                end = next_start - 1
                id = line.split(None, 3)[2].strip()
                self.sv_rec_by_id[id] = (start, end, line_nb)
            line = self.pick_reader.readline()
            line_nb += 1

    def close(self):
        """Close file handle."""
        super().close()
        if self.pick_reader is not None:
            self.pick_reader.close()

    def get(self, id):
        """
        Return record with provided ID.

        :param id: record ID.
        :type id: str
        :return: The record.
        :rtype: anacore.vcf.VCFRecord
        """
        # Get current pos
        current = {
            "line": self.current_line,
            "line_nb": self.current_line_nb,
        }
        # Get new record
        start, end, line_nb = self.sv_rec_by_id[id]
        self.pick_reader.seek(start)
        self.current_line = self.pick_reader.read(end - start + 1)
        self.current_line_nb = line_nb
        record = self._parseLine()
        # Jump to previous record
        self.current_line = current["line"]
        self.current_line_nb = current["line_nb"]
        # Return
        return record

    def getSub(self, chr, start, end):
        """
        Return generator on records overlapping the specified region.

        .. warning::
            This method can only be used with an instance of anacore.VCFIO opened
            with mode "i".

        :param chr: Chromosome name of the selected region.
        :type chr: str
        :param start: Start of the selected region (1-based).
        :type start: int
        :param end: End of the selected region (1-based).
        :type end: int
        :return: Fusions (record and his mates) overlappind the specified region.
        :rtype: generator for (anacore.vcf.VCFRecord, [anacore.vcf.VCFRecord, ...])
        """
        for record in super().getSub(chr, start, end):
            mates = []
            for mate_id in record.info["MATEID"]:
                mates.append(self.get(mate_id))
            yield (record, mates)

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a new record (it is not a comment, an header line or an already processed fusion).

        :param line: The evaluated line.
        :type line: str.
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        is_record = super().isRecordLine(line)
        if is_record:
            id = line.split("\t", 3)[2].strip()
            mate_info = line.split("\t", 8)[7].strip()
            match = re.search("MATEID\=([^;]+);?", mate_info)
            if match:
                mate_id = decodeInfoValue(match.groups()[0])
                fusion_id = " @@ ".join(sorted([id, mate_id]))
                if fusion_id in self._iter_already_processed:  # Partial pre-filtering for records without multi-alternatives
                    is_record = False
            else:
                is_record = False
        return is_record

    def __iter__(self):
        for record in super().__iter__():
            for mate_idx, mate_id in enumerate(record.info["MATEID"]):
                mate = self.get(mate_id)
                fusion_id = " @@ ".join(sorted([record.id, mate.id]))
                if fusion_id not in self._iter_already_processed:
                    self._iter_already_processed.add(fusion_id)
                    # Change ID
                    record_new_id = record.id
                    alt_record = record
                    if len(record.alt) > 1:
                        record_new_id += "_" + str(mate_idx)  # Record must be splitted for each mate
                        alt_record = getAlleleRecord(self, record, mate_idx)
                        alt_record.info["MATEID"] = [record.info["MATEID"][mate_idx]]
                    mate_new_id = mate.id
                    alt_mate = mate
                    if len(mate.alt) > 1:
                        record_idx = mate.info["MATEID"].index(record.id)
                        mate_new_id += "_" + str(record_idx)  # Record must be splitted for each mate
                        alt_mate = getAlleleRecord(self, mate, record_idx)
                        alt_mate.info["MATEID"] = [mate.info["MATEID"][record_idx]]
                    alt_record.id = record_new_id
                    alt_mate.info["MATEID"] = [record_new_id]
                    alt_mate.id = mate_new_id
                    alt_record.info["MATEID"] = [mate_new_id]
                    # Order by transcript
                    first = alt_record
                    second = alt_mate
                    if "RNA_FIRST" in second.info:
                        first = alt_mate
                        second = alt_record
                    # Return
                    yield first, second

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in VCF output format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in VCF output format.
        :rtype: bool
        """
        is_valid = False
        try:
            with VCFIO(filepath) as reader:
                if "SVTYPE" in reader.info and "MATEID" in reader.info:
                    is_valid = True
        except FileNotFoundError:
            raise
        except Exception:
            pass
        return is_valid

    def write(self, first, second):
        """
        Write fusion record in VCF.

        :param first: The fisrt breakend.
        :type first: anacore.vcf.VCFRecord
        :param second: The second breakend.
        :type second: anacore.vcf.VCFRecord
        """
        super().write(first)
        super().write(second)
