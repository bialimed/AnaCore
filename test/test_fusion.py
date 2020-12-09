#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.5.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import pysam
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.vcf import VCFRecord, HeaderFilterAttr, HeaderInfoAttr
from anacore.fusion import *


####################################################################
#
# FUNCTIONS
#
####################################################################
class FusionFileReaderTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_fusions.tsv")
        self.test_cases = [
            {
                "expected": FusionCatcherIO,
                "content": 'Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect'
            },
            {
                "expected": FusionCatcherIO,
                "content": 'Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence'
            },
            {
                "expected": STARFusionIO,
                "content": '#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots'
            },
            {
                "expected": STARFusionIO,
                "content": '#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots'
            },
            {
                "expected": ArribaIO,
                "content": '#gene1	gene2	strand1(gene/fusion)	strand2(gene/fusion)	breakpoint1	breakpoint2	site1	site2	type	direction1	direction2	split_reads1	split_reads2	discordant_mates	coverage1	coverage2	confidence	closest_genomic_breakpoint1	closest_genomic_breakpoint2	filters	fusion_transcript	reading_frame	peptide_sequence	read_identifiers'
            },
            {
                "expected": BreakendVCFIO,
                "content": '''##fileformat=VCFv4.3
##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|Gene|Feature|Feature_type|STRAND|EXON|INTRON|CDS_position|Protein_position|GENE_SHARD">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=IDSRC,Number=.,Type=String,Description="ID of breakend by source">
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this break-end is 5prim in the fusion transcript.">
##INFO=<ID=SRC,Number=.,Type=String,Description="Fusions callers where the breakend is identified.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=s0_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local translocation breakend",Source="manta">
##INFO=<ID=s0_CIEND,Number=2,Type=Integer,Description="Confidence interval around END",Source="manta">
##INFO=<ID=s0_CIGAR,Number=A,Type=String,Description="CIGAR alignment for each alternate indel allele",Source="manta">
##INFO=<ID=s0_CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS",Source="manta">
##INFO=<ID=s0_END,Number=1,Type=Integer,Description="End position of the variant described in this record",Source="manta">
##INFO=<ID=s0_EVENT,Number=1,Type=String,Description="ID of event associated to breakend",Source="manta">
##INFO=<ID=s0_HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints",Source="manta">
##INFO=<ID=s0_HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical homology at event breakpoints",Source="manta">
##INFO=<ID=s0_IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation",Source="manta">
##INFO=<ID=s0_LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length",Source="manta">
##INFO=<ID=s0_MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote translocation mate breakend",Source="manta">
##INFO=<ID=s0_MATE_REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at the other breakend",Source="manta">
##INFO=<ID=s0_REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at this breakend",Source="manta">
##INFO=<ID=s0_RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length",Source="manta">
##INFO=<ID=s0_RNA_CONTIG,Number=1,Type=String,Description="The sequence of the breakend spanning contig",Source="manta">
##INFO=<ID=s0_RNA_CONTIG_ALN,Number=2,Type=Integer,Description="Length of the spanning contig alignment on each breakend",Source="manta">
##INFO=<ID=s0_RNA_FwRvReads,Number=2,Type=Integer,Description="For RNA fusions, number of stranded reads supporting forward or reverse direction of transcription",Source="manta">
##INFO=<ID=s0_RNA_Reads,Number=1,Type=Integer,Description="The number of reads and pairs that potentially support this candidate before refinement and scoring",Source="manta">
##INFO=<ID=s0_RNA_STRANDED,Number=0,Type=Flag,Description="For RNA fusions, the direction of transcription is known",Source="manta">
##INFO=<ID=s0_SVINSLEN,Number=.,Type=Integer,Description="Length of insertion",Source="manta">
##INFO=<ID=s0_SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion",Source="manta">
##INFO=<ID=s0_SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles",Source="manta">
##INFO=<ID=s0_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="manta">
##INFO=<ID=s1_BREAK_DINUC,Number=1,Type=String,Description="Dinucleotides flanking the breakend of the fragment excluded by the fusion.",Source="STAR_Fusion">
##INFO=<ID=s1_BREAK_ENTROPY,Number=1,Type=Float,Description="Shannon entropy of the 15 exonic bases flanking the breakpoint. The maximum entropy is 2, representing highest complexity. The lowest would be zero (involving a 15 base mononucleotide run). Low entropy sites should generally be treated as less confident breakpoints.",Source="STAR_Fusion">
##INFO=<ID=s1_FCANN,Number=.,Type=String,Description="Annotation generated by FusionAnnotator (see: https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki). Format: SYMBOL|Gene|Tags",Source="STAR_Fusion">
##INFO=<ID=s1_SPLICE_TYPE,Number=1,Type=String,Description="Whether the proposed breakpoint occurs at reference exon junctions as provided by the reference transcript structure annotations (ex. gencode).",Source="STAR_Fusion">
##INFO=<ID=s1_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="STAR_Fusion">
##INFO=<ID=s2_FCANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|STRAND|Site|Type|GENE_SHARD|FRAMESHIFT|Protein_contig",Source="Arriba">
##INFO=<ID=s2_GBP,Number=1,Type=String,Description="The coordinates of the genomic breakpoint which is closest to the transcriptomic breakpoint.",Source="Arriba">
##INFO=<ID=s2_IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation.",Source="Arriba">
##INFO=<ID=s2_RNA_CONTIG,Number=1,Type=String,Description="The transcript sequence assembled from the supporting reads of the most highly expressed transcript.",Source="Arriba">
##INFO=<ID=s2_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="Arriba">
##INFO=<ID=s3_FCANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|Gene|EXON|Effect.",Source="FusionCatcher">
##INFO=<ID=s3_RNA_CONTIG,Number=1,Type=String,Description="The sequence of the breakend spanning contig.",Source="FusionCatcher">
##INFO=<ID=s3_SOURCES,Number=.,Type=String,Description="Aligning method used for mapping the reads and finding the fusion genes. Here are two methods used which are: (i) BOWTIE = only Bowtie aligner is used for mapping the reads on the genome and exon-exon fusion junctions, (ii) BOWTIE+BLAT = Bowtie aligner is used for mapping reads on the genome and BLAT is used for mapping reads for finding the fusion junction, (iii) BOWTIE+STAR = Bowtie aligner is used for mapping reads on the genome and STAR is used for mapping reads for finding the fusion junction, (iv) BOWTIE+BOWTIE2 = Bowtie aligner is used for mapping reads on the genome and Bowtie2 is used for mapping reads for finding the fusion junction.",Source="FusionCatcher">
##INFO=<ID=s3_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="FusionCatcher">
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
##FILTER=<ID=s0_Imprecise,Description="RNA fusion candidates for which no spanning contig was found",Source="manta">
##FILTER=<ID=s0_Local,Description="RNA call covering short genomic distance",Source="manta">
##FILTER=<ID=s0_LowEvidence,Description="RNA fusion calls without both split read and spanning pair support",Source="manta">
##FILTER=<ID=s1_Babiceanu_Normal,Description="Recurrent chimeric fusion RNAs in non-cancer tissues and cells as per Babiceanu et al. NAR, 2016",Source="STAR_Fusion">
##FILTER=<ID=s1_BodyMap,Description="Fusions found by STAR-Fusion as applied to the Illumina Human Body Map reference data",Source="STAR_Fusion">
##FILTER=<ID=s1_ConjoinG,Description="Fused transcripts derived from the Conjoined Gene Database",Source="STAR_Fusion">
##FILTER=<ID=s1_DGD_PARALOGS,Description="Duplicated genes as per the Duplicated Genes Database",Source="STAR_Fusion">
##FILTER=<ID=s1_GTEx_recurrent_StarF2019,Description="Fusions found recurrently in GTEx as mined using STAR-Fusion v1.5.0",Source="STAR_Fusion">
##FILTER=<ID=s1_Greger_Normal,Description="Fusion transcripts (mostly from tandem genes) detected based on analysis of RNA-Seq from 1000 genomes project samples. List derived from Greger et al. PLOS One, 2014",Source="STAR_Fusion">
##FILTER=<ID=s1_HGNC_GENEFAM,Description="HGNC gene family membership as per ftp://ftp.ebi.ac.uk/pub/databases/genenames/genefam_list.txt.gz",Source="STAR_Fusion">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Count of pairs of reads supporting the fusion">
##FORMAT=<ID=PRSRC,Number=.,Type=Integer,Description="Count of pairs of reads supporting the fusion by source">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Count of reads mapping on the fusion junction">
##FORMAT=<ID=SRSRC,Number=.,Type=Integer,Description="Count of reads mapping on the fusion junction by source">
##FORMAT=<ID=s0_PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed",Source="manta">
##FORMAT=<ID=s0_SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed",Source="manta">
##FORMAT=<ID=s1_FFPM,Number=1,Type=Float,Description="Normalized measures of the fusion-supporting rna-seq fragments (fusion fragments per million total reads).",Source="STAR_Fusion">
##FORMAT=<ID=s1_JRL,Number=1,Type=String,Description="RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.",Source="STAR_Fusion">
##FORMAT=<ID=s1_PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (SpanningFragCount).",Source="STAR_Fusion">
##FORMAT=<ID=s1_SFL,Number=1,Type=String,Description="RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.",Source="STAR_Fusion">
##FORMAT=<ID=s1_SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction (JunctionReadCount).",Source="STAR_Fusion">
##FORMAT=<ID=s1_hasLAS,Number=1,Type=String,Description="This column indicates whether there are split reads that provide long (set to length of 25 bases) alignments on both sides of the putative breakpoint (LargeAnchorSupport).",Source="STAR_Fusion">
##FORMAT=<ID=s2_CFD,Number=1,Type=String,Description="Each prediction is assigned one of the confidences low, medium, or high. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations.",Source="Arriba">
##FORMAT=<ID=s2_DPS,Number=1,Type=Integer,Description="Coverage near breakpoint. The coverage is calculated as the number of fragments near the breakpoint on the side of the breakpoint that is retained in the fusion transcript. Note that the coverage calculation counts all fragments (even duplicates).",Source="Arriba">
##FORMAT=<ID=s2_PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.",Source="Arriba">
##FORMAT=<ID=s2_RFIL,Number=.,Type=String,Description="Filters which removed one or more of the supporting reads. The number of filtered reads is given in parantheses after the name of the filter. If a filter discarded the event as a whole (all reads), the number of filtered reads is missing.",Source="Arriba">
##FORMAT=<ID=s2_SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.",Source="Arriba">
##FORMAT=<ID=s2_SR1,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on first shard.",Source="Arriba">
##FORMAT=<ID=s2_SR2,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on second shard.",Source="Arriba">
##FORMAT=<ID=s2_SRL,Number=.,Type=String,Description="The names of the supporting reads.",Source="Arriba">
##FORMAT=<ID=s3_CM,Number=1,Type=Integer,Description="Count of reads mapping simultaneously on both genes which form the fusion gene. This is an indication how similar are the DNA/RNA sequences of the genes forming the fusion gene (i.e. what is their homology because highly homologous genes tend to appear show as candidate fusion genes). In case of completely different sequences of the genes involved in forming a fusion gene then here it is expected to have the value zero.",Source="FusionCatcher">
##FORMAT=<ID=s3_LA,Number=1,Type=Integer,Description="Longest anchor (hangover) found among the unique reads mapping on the fusion junction.",Source="FusionCatcher">
##FORMAT=<ID=s3_PR,Number=1,Type=Integer,Description="Count of pairs of reads supporting the fusion (including also the multimapping reads).",Source="FusionCatcher">
##FORMAT=<ID=s3_SR,Number=1,Type=Integer,Description="Count of unique reads (i.e. unique mapping positions) mapping on the fusion junction. Shortly, here are counted all the reads which map on fusion junction minus the PCR duplicated reads.",Source="FusionCatcher">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A'''
            }
        ]
        self.test_exceptions = [
            "",
            'Gene	Count',
            'ATGCATTGGCAT'
        ]

    def tearDown(self):
        # Clean temporary files
        if os.path.exists(self.tmp_in):
            os.remove(self.tmp_in)

    def test_factory(self):
        for curr_test in self.test_cases:
            with open(self.tmp_in, "w") as writer:
                writer.write(curr_test["content"])
            self.assertEqual(
                curr_test["expected"],
                FusionFileReader.factory(self.tmp_in).__class__
            )

    def test_factoryExceptions(self):
        for curr_content in self.test_exceptions:
            with open(self.tmp_in, "w") as writer:
                writer.write(curr_content)
            with self.assertRaises(IOError) as context:
                FusionFileReader.factory(self.tmp_in)

    def test_factoryParams(self):
        for curr_test in self.test_cases:
            if curr_test["expected"] != BreakendVCFIO:
                with open(self.tmp_in, "w") as writer:
                    writer.write(curr_test["content"])
                with FusionFileReader.factory(self.tmp_in, "r", "AnnotField", sample_name="splA") as reader:
                    self.assertEqual(
                        curr_test["expected"],
                        reader.__class__
                    )
                    self.assertEqual(
                        "splA",
                        reader.sample_name
                    )
                    self.assertEqual(
                        "AnnotField",
                        reader.annot_field
                    )


class FusionCatcherIOTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        # Temporary files
        self.tmp_in_tsv = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out_tsv = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Create in TSV
        content = '''Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect
RREB1	MKL2	oncogene,cancer,exon-exon	0	9	10	30	BOWTIE	6:7211709:+	16:14245528:+	ENSG00000124782	ENSG00000186260	ENSE00003408308	ENSE00001202958	GACTGGAGACCCACATGGAGACCCATTCAGATAACCCACTAAG*GCCACTCAATGACAAAAATAGTAACAGTGGGAATTCAGCTTTG	in-frame'''
        with open(self.tmp_in_tsv, "w") as writer:
            writer.write(content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_tsv, self.tmp_out_tsv]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_ReadWrite(self):
        # Read, convert to VCF, covert to dict and write
        with FusionCatcherIO(self.tmp_out_tsv, "w") as writer:
            with FusionCatcherIO(self.tmp_in_tsv) as reader:
                for first_bnd, second_bnd in reader:
                    writer.write((first_bnd, second_bnd))
        # Eval
        expected = None
        with open(self.tmp_out_tsv) as expec_reader:
            expected = expec_reader.readlines()
        observed = None
        with open(self.tmp_out_tsv) as obs_reader:
            observed = obs_reader.readlines()
        self.assertEqual(expected, observed)


class STARFusionIOTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        # Temporary files
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_in_abridged = os.path.join(tmp_folder, unique_id + "_in_abridged.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Create in TSV
        content = '''#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots
MN1--BEND2	90	39	ONLY_REF_SPLICE	MN1^ENSG00000169184.6	chr22:27796763:-	BEND2^ENSG00000177324.14	chrX:18195442:-	M70265:74:000000000-B4F4J:1:2116:12294:7575,M70265:74:000000000-B4F4J:1:2115:7184:4292,M70265:74:000000000-B4F4J:1:1104:5050:21345,M70265:74:000000000-B4F4J:1:1119:9539:18627,M70265:74:000000000-B4F4J:1:2115:8678:9709,M70265:74:000000000-B4F4J:1:2119:19014:22293,M70265:74:000000000-B4F4J:1:2105:16442:9076,M70265:74:000000000-B4F4J:1:1117:12801:7721,M70265:74:000000000-B4F4J:1:1111:14401:18304,M70265:74:000000000-B4F4J:1:1112:22240:14651,M70265:74:000000000-B4F4J:1:1108:20469:22260,M70265:74:000000000-B4F4J:1:2102:11105:8918,M70265:74:000000000-B4F4J:1:2108:19409:16466,M70265:74:000000000-B4F4J:1:1113:20196:5091,M70265:74:000000000-B4F4J:1:1105:10766:13643,M70265:74:000000000-B4F4J:1:1103:19101:13070,M70265:74:000000000-B4F4J:1:2107:24460:10546,M70265:74:000000000-B4F4J:1:1106:19866:14418,M70265:74:000000000-B4F4J:1:1106:5053:13766,M70265:74:000000000-B4F4J:1:2113:20643:17000,M70265:74:000000000-B4F4J:1:2102:25631:11285,M70265:74:000000000-B4F4J:1:1119:27699:19981,M70265:74:000000000-B4F4J:1:2113:4962:18403,M70265:74:000000000-B4F4J:1:2106:25372:20356,M70265:74:000000000-B4F4J:1:2108:16489:16372,M70265:74:000000000-B4F4J:1:1105:14999:8683,M70265:74:000000000-B4F4J:1:1107:13181:5229,M70265:74:000000000-B4F4J:1:1103:9694:18262,M70265:74:000000000-B4F4J:1:1106:7893:19397,M70265:74:000000000-B4F4J:1:1119:28664:16235,M70265:74:000000000-B4F4J:1:1119:25900:13031,M70265:74:000000000-B4F4J:1:2117:5341:12350,M70265:74:000000000-B4F4J:1:2108:21016:10605,M70265:74:000000000-B4F4J:1:1111:6184:20952,M70265:74:000000000-B4F4J:1:1111:7115:5525,M70265:74:000000000-B4F4J:1:1110:27378:16651,M70265:74:000000000-B4F4J:1:1108:14118:21000,M70265:74:000000000-B4F4J:1:2110:8046:7787,M70265:74:000000000-B4F4J:1:2112:29162:12358,M70265:74:000000000-B4F4J:1:1114:12762:7330,M70265:74:000000000-B4F4J:1:1112:6478:12309,M70265:74:000000000-B4F4J:1:2114:8353:4989,M70265:74:000000000-B4F4J:1:1110:12165:19870,M70265:74:000000000-B4F4J:1:2119:26398:14585,M70265:74:000000000-B4F4J:1:1104:19297:14091,M70265:74:000000000-B4F4J:1:1114:12488:10493,M70265:74:000000000-B4F4J:1:2110:19477:14905,M70265:74:000000000-B4F4J:1:1111:29396:14060,M70265:74:000000000-B4F4J:1:1112:5997:8955,M70265:74:000000000-B4F4J:1:1114:21233:24893,M70265:74:000000000-B4F4J:1:1117:16546:12881,M70265:74:000000000-B4F4J:1:2102:15683:16355,M70265:74:000000000-B4F4J:1:1102:7568:7922,M70265:74:000000000-B4F4J:1:1104:24497:10850,M70265:74:000000000-B4F4J:1:1112:12330:15917,M70265:74:000000000-B4F4J:1:2104:23488:20289,M70265:74:000000000-B4F4J:1:2115:18539:5639,M70265:74:000000000-B4F4J:1:2108:27004:11103,M70265:74:000000000-B4F4J:1:1108:7280:6580,M70265:74:000000000-B4F4J:1:1111:13093:9434,M70265:74:000000000-B4F4J:1:1105:22026:13215,M70265:74:000000000-B4F4J:1:1103:12233:24724,M70265:74:000000000-B4F4J:1:1114:21833:4508,M70265:74:000000000-B4F4J:1:1102:5182:21387,M70265:74:000000000-B4F4J:1:2102:20234:8226,M70265:74:000000000-B4F4J:1:1104:5911:7048,M70265:74:000000000-B4F4J:1:1104:19911:7039,M70265:74:000000000-B4F4J:1:2101:23116:15940,M70265:74:000000000-B4F4J:1:1102:28272:19048,M70265:74:000000000-B4F4J:1:2107:3513:10128,M70265:74:000000000-B4F4J:1:2101:3654:8461,M70265:74:000000000-B4F4J:1:1119:19110:18013,M70265:74:000000000-B4F4J:1:1107:18400:1636,M70265:74:000000000-B4F4J:1:2115:15123:16288,M70265:74:000000000-B4F4J:1:2110:13483:19384,M70265:74:000000000-B4F4J:1:2101:9347:11364,M70265:74:000000000-B4F4J:1:1103:14624:5192,M70265:74:000000000-B4F4J:1:2103:25588:18903,M70265:74:000000000-B4F4J:1:1117:15805:5600,M70265:74:000000000-B4F4J:1:1119:13495:19891,M70265:74:000000000-B4F4J:1:1108:7984:5492,M70265:74:000000000-B4F4J:1:1117:9253:13094,M70265:74:000000000-B4F4J:1:1112:8428:11232,M70265:74:000000000-B4F4J:1:2103:22447:6340,M70265:74:000000000-B4F4J:1:1111:26219:5015,M70265:74:000000000-B4F4J:1:1112:14939:4411,M70265:74:000000000-B4F4J:1:2109:18888:3219,M70265:74:000000000-B4F4J:1:2110:22338:16190,M70265:74:000000000-B4F4J:1:2114:11325:2713,M70265:74:000000000-B4F4J:1:2114:6567:13905	M70265:74:000000000-B4F4J:1:1111:7154:15877,M70265:74:000000000-B4F4J:1:2117:22762:16468,M70265:74:000000000-B4F4J:1:1112:3545:11649,M70265:74:000000000-B4F4J:1:1102:13385:24356,M70265:74:000000000-B4F4J:1:1103:21865:16107,M70265:74:000000000-B4F4J:1:2109:19125:3230,M70265:74:000000000-B4F4J:1:1113:11062:16439,M70265:74:000000000-B4F4J:1:2115:11207:19666,M70265:74:000000000-B4F4J:1:2102:1922:14639,M70265:74:000000000-B4F4J:1:2105:29283:14469,M70265:74:000000000-B4F4J:1:1109:15243:9627,M70265:74:000000000-B4F4J:1:1108:22583:13263,M70265:74:000000000-B4F4J:1:1103:23453:4405,M70265:74:000000000-B4F4J:1:1113:9694:20456,M70265:74:000000000-B4F4J:1:2103:19632:8012,M70265:74:000000000-B4F4J:1:1102:18799:15205,M70265:74:000000000-B4F4J:1:1102:19075:11592,M70265:74:000000000-B4F4J:1:2114:28568:18618,M70265:74:000000000-B4F4J:1:1112:21951:11679,M70265:74:000000000-B4F4J:1:1111:22886:13062,M70265:74:000000000-B4F4J:1:1115:17822:21329,M70265:74:000000000-B4F4J:1:1116:21890:11280,M70265:74:000000000-B4F4J:1:1104:10237:9204,M70265:74:000000000-B4F4J:1:1112:3494:8607,M70265:74:000000000-B4F4J:1:2103:5353:4657,M70265:74:000000000-B4F4J:1:2101:18273:12696,M70265:74:000000000-B4F4J:1:2109:13674:17280,M70265:74:000000000-B4F4J:1:1103:18734:3105,M70265:74:000000000-B4F4J:1:2102:17275:5848,M70265:74:000000000-B4F4J:1:1112:24048:22285,M70265:74:000000000-B4F4J:1:1115:21982:9046,M70265:74:000000000-B4F4J:1:1109:27605:16863,M70265:74:000000000-B4F4J:1:1111:14542:18181,M70265:74:000000000-B4F4J:1:1104:9119:13899,M70265:74:000000000-B4F4J:1:1114:10485:22262,M70265:74:000000000-B4F4J:1:2111:24444:4566,M70265:74:000000000-B4F4J:1:2113:24284:23644,M70265:74:000000000-B4F4J:1:2104:27141:8778,M70265:74:000000000-B4F4J:1:2117:13953:8079	YES_LDAS	38.9246	GT	1.4566	AG	1.7968	["HaasMedCancer","INTERCHROMOSOMAL[chr22--chrX]"]
SPECC1L-ADORA2A--NTRK2	22	6	ONLY_REF_SPLICE	SPECC1L-ADORA2A^ENSG00000258555.6	chr22:24334573:+	NTRK2^ENSG00000148053.16	chr9:84867243:+	M70265:74:000000000-B4F4J:1:1102:11424:4669,M70265:74:000000000-B4F4J:1:1104:6405:2923,M70265:74:000000000-B4F4J:1:1114:21462:20506,M70265:74:000000000-B4F4J:1:1108:9784:2773,M70265:74:000000000-B4F4J:1:2119:22961:20082,M70265:74:000000000-B4F4J:1:2116:27522:11441,M70265:74:000000000-B4F4J:1:1108:23971:7984,M70265:74:000000000-B4F4J:1:1113:23518:22366,M70265:74:000000000-B4F4J:1:1115:7837:20035,M70265:74:000000000-B4F4J:1:2104:21327:2495,M70265:74:000000000-B4F4J:1:1101:9462:22886,M70265:74:000000000-B4F4J:1:1108:5888:7340,M70265:74:000000000-B4F4J:1:2103:6756:2599,M70265:74:000000000-B4F4J:1:1114:10446:14926,M70265:74:000000000-B4F4J:1:1108:4192:20832,M70265:74:000000000-B4F4J:1:2110:22090:10259,M70265:74:000000000-B4F4J:1:1107:22426:22529,M70265:74:000000000-B4F4J:1:1102:20172:13788,M70265:74:000000000-B4F4J:1:1102:14526:14874,M70265:74:000000000-B4F4J:1:1108:28140:12698,M70265:74:000000000-B4F4J:1:1109:25864:11274,M70265:74:000000000-B4F4J:1:1109:21996:25012	M70265:74:000000000-B4F4J:1:1108:25374:7116,M70265:74:000000000-B4F4J:1:1102:12684:14754,M70265:74:000000000-B4F4J:1:1102:12917:5068,M70265:74:000000000-B4F4J:1:1101:28832:12042,M70265:74:000000000-B4F4J:1:1101:26215:18199,M70265:74:000000000-B4F4J:1:2110:18851:11245	YES_LDAS	8.3732	GT	1.9656	AG	1.8295	["INTERCHROMOSOMAL[chr22--chr9]"]
SPECC1L-ADORA2A--NTRK2	3	6	ONLY_REF_SPLICE	SPECC1L-ADORA2A^ENSG00000258555.6	chr22:24334573:+	NTRK2^ENSG00000148053.16	chr9:84861040:+	M70265:74:000000000-B4F4J:1:1108:15540:9098,M70265:74:000000000-B4F4J:1:2104:21069:2777,M70265:74:000000000-B4F4J:1:1116:20816:4555	M70265:74:000000000-B4F4J:1:1101:28832:12042,M70265:74:000000000-B4F4J:1:1102:12917:5068,M70265:74:000000000-B4F4J:1:1101:26215:18199,M70265:74:000000000-B4F4J:1:1108:25374:7116,M70265:74:000000000-B4F4J:1:1102:12684:14754,M70265:74:000000000-B4F4J:1:2110:18851:11245	YES_LDAS	2.6914	GT	1.9656	AG	1.7232	["INTERCHROMOSOMAL[chr22--chr9]"]'''
        with open(self.tmp_in, "w") as writer:
            writer.write(content)
        with open(self.tmp_in_abridged, "w") as writer:
            for line in content.split("\n"):
                fields = line.split("\t")
                abridged_line = "\t".join(fields[:8] + fields[10:])
                writer.write(abridged_line + "\n")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in, self.tmp_in_abridged, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_ReadWrite(self):
        # Read, convert to VCF, convert to dict and write
        with STARFusionIO(self.tmp_out, "w") as writer:
            with STARFusionIO(self.tmp_in) as reader:
                for first_bnd, second_bnd in reader:
                    writer.write((first_bnd, second_bnd))
        # Eval
        expected = None
        with open(self.tmp_out) as expec_reader:
            expected = expec_reader.readlines()
        observed = None
        with open(self.tmp_out) as obs_reader:
            observed = obs_reader.readlines()
        self.assertEqual(expected, observed)

    def test_ReadWriteAbridged(self):
        # Read, convert to VCF, convert to dict and write
        with STARFusionIO(self.tmp_out, "w") as writer:
            with STARFusionIO(self.tmp_in_abridged) as reader:
                for first_bnd, second_bnd in reader:
                    writer.write((first_bnd, second_bnd))
        # Eval
        expected = None
        with open(self.tmp_out) as expec_reader:
            expected = []
            for idx, line in enumerate(expec_reader.readlines()):
                if idx != 0:
                    fields = line.split("\t")
                    fields[8] = ""
                    fields[9] = ""
                    line = "\t".join(fields)
                expected.append(line)
        observed = None
        with open(self.tmp_out) as obs_reader:
            observed = obs_reader.readlines()
        self.assertEqual(expected, observed)


class ArribaIOTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        # Temporary files
        self.tmp_in_tsv = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out_tsv = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Create in TSV
        content = """#gene1	gene2	strand1(gene/fusion)	strand2(gene/fusion)	breakpoint1	breakpoint2	site1	site2	type	direction1	direction2	split_reads1	split_reads2	discordant_mates	coverage1	coverage2	confidence	closest_genomic_breakpoint1	closest_genomic_breakpoint2	filters	fusion_transcript	reading_frame	peptide_sequence	read_identifiers
BCR	ABL1	+/+	+/+	22:23289621	9:130854064	splice-site	splice-site	translocation	downstream	upstream	12	165	84	386	1134	high	.	.	duplicates(63)	CGGAGAGGCTGAAGAAGAAGCTGTCGGAGCAGGAGTCACTGCTGCTGCTTATGTCTCCCAGCATGGCCTTCAGGGTGCACAGCCGCAACGGCAAG___AGTTACACGTTCCTGATCTCCTCTGACTATGAGCGTGCAGAGTGGAGGGAGAACATCCGGGAGCAGCAGAAGAAGT___GTTTCAGAAGCTTCTCCCTGACATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAAcAAGGAAG|AAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAG___GTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTGTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAATGGCAGCTTCTTGGTGCGTGAGAGTGAGAG	in-frame	ERLKKKLSEQESLLLLMSPSMAFRVHSRNGKSYTFLISSDYERAEWRENIREQQKKCFRSFSLTSVELQMLTNSCVKLQTVHSIPLTINKE|eALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESE	.
ABL1	BCR	+/+	+/+	9:130714455	22:23292541	splice-site	splice-site	translocation	downstream	upstream	9	9	1	188	433	high	.	.	duplicates(1)	GAAACTTTTGGGTAACATTGCAATTACATGAAATTGATAACCGCGAAAATAATTGGAACTCCTGCTTGCAAGTGTCAACCTAAAAAAAGTGCTTCCTTTTGTTATGGAAGATGTCTTTCTGTGATTGACTTCAATTGCTGACTTGTGGAGATGCAGCGAATGTGAAATCCCACGTATATGCCATTTCCCTCTACGCTCGCTGACCGTTCTGGAAGATCTTGAACCCTCTTCTGGAAAGGGGTACCTATTATTACTTTATGGGGCAGCAGCCTGGAAAAGTACTTGGGGACCAAAGAAGGCCAAGCTTGCCTGCCCTGCATTTTATCAAAGGAGCAGGGAAGAAGGAATCATCGAGGCATGGGGGTCCACACTGCAATGTTTTTGTGGAACATG|ATCTGTACTGCACCCTGGAGGTGGATTCCTTTGGGTATTTTGTGAATAAAGCAAAGACGCGCGTCTACAGGGACACAGCTGAGCCAAACTGGAACGAG___GAATTTGAGATAGAGCTGGAGGGCTCCCAGACCCTGAGGATACTGTGCTATGAAAAGTGTTACAACAAGACGAAGATCCCCAAGGAGGACGGCGAGAGCACGG...___CTGGACCCGCAGGCCCTGCAGGACAGAGACTGGCAGCGCACCGTCATCGCCATGAATGGG	in-frame	MGQQPGKVLGDQRRPSLPALHFIKGAGKKESSRHGGPHCNVFVEH|dLYCTLEVDSFGYFVNKAKTRVYRDTAEPNWNEEFEIELEGSQTLRILCYEKCYNKTKIPKEDGEST	.
ABL1	BCR	+/+	+/+	9:130835525	22:23290339	splice-site	splice-site	translocation	downstream	upstream	1	4	1	27	328	high	.	.	.	CGGGCCTGAGCCGGGCCCGCGGACCGAGCTtGGAGAGGGGcTCCGGCCCCCGACGTGCTGGCGCGGGAAAATGTTGGAGATCTGCCTGAAGCTGGTGGGCTGCAAATCCAAGAAGGGGCTGTCCTCGTCCTCCAGCTGTTATCTGGAAG|ATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAA___ATCTGTACTGCACCCTGGAGGTGGATTCCTTTGGGTATTTTGTGAATAAAGCAAAGACGCGCGTCTACAGGGACACAGCTGAGCCAAACTGGAACGAG___GAATTTGAGATAGAGCTGGAGGGCTCCCAGACCCTGAGGATACTGTGCTATGAAAAGTGTTACAACAAGA	in-frame	MLEICLKLVGCKSKKGLSSSSSCYLE|dDESPGLYGFLNVIVHSATGFKQSSNLYCTLEVDSFGYFVNKAKTRVYRDTAEPNWNEEFEIELEGSQTLRILCYEKCYNK	.
ABL1	BCR	+/+	+/+	9:130714455	22:23290339	splice-site	splice-site	translocation	downstream	upstream	0	4	1	188	328	high	.	.	.	CCGTTCTGGAAGATCTTGAACCCTCTTCTGGAAAGGGGTACCTATTATTACTTTATGGGGCAGCAGCCTGGAAAAGTACTTGGGGACCAAAGAAGGCCAAGCTTGCCTGCCCTGCATTTTATCAAAGGAGCAGGGAAGAAGGAATCATCGAGGCATGGGGGTCCACACTGCAATGTTTTTGTGGAACATG|ATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAA___ATCTGTACTGCACCCTGGAGGTGGATTCCTTTGGGTATTTTGTGAATAAAGCAAAGACGCGCGTCTACAGGGACACAGCTGAGCCAAACTGGAACGAGG...___GAATTTGAGATAGAGC	in-frame	MGQQPGKVLGDQRRPSLPALHFIKGAGKKESSRHGGPHCNVFVEH|dDESPGLYGFLNVIVHSATGFKQSSNLYCTLEVDSFGYFVNKAKTRVYRDTAEPNWNE	.
BCR	ABL1	+/+	+/+	22:23289621	9:130854067	splice-site	splice-site	translocation	downstream	upstream	0	1	77	386	1134	medium	.	.	duplicates(14),inconsistently_clipped(1)	TATGTCTCCCAGCATGGCCTTCAGGGTGCACAGCCGCAACGGCAAG___AGTTACACGTTCCTGATCTCCTCTGACTATGAGCGTGCAGAGTGGAGGGAGAACATCCGGGAGCAGCAGAAGAAGT___GTTTCAGAAGCTTCTCCCTGACATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAAcAAGGAAG|CCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAG___GTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTGTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTC	in-frame	MSPSMAFRVHSRNGKSYTFLISSDYERAEWRENIREQQKKCFRSFSLTSVELQMLTNSCVKLQTVHSIPLTINKE|aLQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITP	.
ABL1	BCR	+/+	+/+	9:130835525	22:23292541	splice-site	splice-site	translocation	downstream	upstream	0	3	1	27	433	medium	.	.	.	GTTGGAGATCTGCCTGAAGtTGGTGGGCTGCAAATCCAAGAAGGGGCTGTCCTCGTCCTCCAGCTGTTATCTGGAAG|ATCTGTACTGCACCCTGGAGGTGGATTCCTTTGGGTATTTTGTGAATAAAGCAAAGACGCGCGTCTACAGGGACACAGCTGAGCCAAACTGGAACGAG___GAATTTGAGATAGAGCTGGAGGGCTCCCAGACCCTGAGGATACTGTGCTATGAAAAGTGTTACAACAAGACGAAGATCCCCA	in-frame	LEICLKLVGCKSKKGLSSSSSCYLE|dLYCTLEVDSFGYFVNKAKTRVYRDTAEPNWNEEFEIELEGSQTLRILCYEKCYNKTKIP	.
ABL1	BCR	+/+	+/+	9:130714455	22:23309424	splice-site	splice-site	translocation	downstream	upstream	0	1	0	188	197	low	.	.	.	GGAAGAAGGAATCATCGAGGCATGGGGGTCCACACTGCAATGTTTTTGTGGAACATG|CTGGACCCGCAGGCCCTGCAGGACAGAGACTGGCAGCGCACCGTCATCGCCATGAATGGG___ATCGAAGTAAAGCTCTCGGTCAAGTTCAACAGCAGGGAGTTCAGCTTGAAGAGGATGCCGTCCCGAAAACAGAC	out-of-frame	KKESSRHGGPHCNVFVEH|agpagpagqrlaahrhrhewdrskalgqvqqqgvqleedavpkt	.
FP671120.1	NAMPT	+/+	-/-	21:8216803	7:106284907	intron	5'UTR	translocation	downstream	downstream	0	1	18	6	8879	medium	.	.	duplicates(4),mismappers(1)	CCCGGGGAGAGTTCTCTTTTCTTTGTGAAGGGCAGGGCGCCCTGGAATGGGTTCGCCCCGAGAGAGGGGCCCGTGCCTTGtAAAGCGTCGCGGTTCCGGCGGCGTCCGGTGAGCTCTCGCTGGCC...AGGTAAGGGAAGTCGGCAAGCCGGATCCGTAACTTCGGGATAAGGATTGGCTCTAAGGGCTGGGTCGGTCGGGCTGGGGCGCGAAGCG|CGGCCCCTGTCCTCCGGCCCGAGATGAATCCTGCGGCAGAAGCCGAGTTCAACATCCTCCTGGCCACCGACTCCTACAAG___GTTACTCACTATAAACAATATCCACCCAACACAAGCAAAGTTTATTCCTACTTTGAATGCCGTGAAAAGAAGACAG	.	.	.
NAMPT	FP671120.1	-/-	+/+	7:106253077	21:8215509	CDS	intron	translocation	upstream	upstream	1	0	16	65535	81	medium	.	.	duplicates(1)	TTTCTTTATTACAAACGACTATCCTTCAAAACAGGCTAAAGTGTTATGTTTCTGAAAGCCAGATGTGGATTGGTTCAGTTATACCTATTAAGTGCTTAAAACTTAAATGTAAATTATGTATGTGC...AAGGGCCGATTATCTTTACATAGGACG|GCGAAAGACTAATCGAACCATCTAGTAGCTGGTTCCCTCCGAAGTTTCCCTCAGGATAGCTGGCGCTCTCGCAGACCCGACGCACCCCCGCCACGCAGTTTTATCCGGTAAAGCGAATGATTAGAGGTCTTGGGGCCGAAACGATCTCAACCTATTCTCAAACTTTAAATGGGTAAGAAGCCCGGCTCGCTGGCGTGGAGCCGGGCGTGGAATGCGAGTGCCTAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGCGGGATGAACCGAACGCCGGGTTAAGGCGCCCGATGCCGACGCTCATCAGACCCCA	out-of-frame	KGRLSLHRT|akd*	.
RUNX1	RUNX1T1	-/-	-/-	21:34859474	8:92017363	splice-site	splice-site	translocation	upstream	downstream	0	9	3	5513	24	medium	.	.	inconsistently_clipped(2)	GTTTGTCGGTCGAAGTGGAAGAG___GGAAAAGCTTCACTCTGACCATCACTGTCTTCACAAACCCACCGCAAGTCGCCACCTACCACAGAGCCATCAAAATCACAGTGGATGGGCCCCGAGAACCTCGAA|ATCGTACTGAGAAGCACTCCACAATGCCAGACTCACCTGTGGATGTGAAGACGCAATCTAGGCTGACTCCTCCAACAATGCCACCTCCCCCAACTACTCAAGGAGCTCCAAGAACCAGTTCATTTACACCGACAACGT___TAACTAATGGCACGAGCCATTCTCCTACAGCCTTGAATGGCGCCCCCTCACCAC	in-frame	FVGRSGRGKSFTLTITVFTNPPQVATYHRAIKITVDGPREPR|nRTEKHSTMPDSPVDVKTQSRLTPPTMPPPPTTQGAPRTSSFTPTTLTNGTSHSPTALNGAPSP	.
RUNX1	RUNX1T1	-/-	-/-	21:34859474	8:92017431	splice-site	splice-site	translocation	upstream	downstream	0	1	3	5513	8	medium	.	.	.	GTTTGTCGGTCGAAGTGGAAGAG___GGAAAAGCTTCACTCTGACCATCACTGTCTTCACAAACCCACCGCAAGTCGCCACCTACCACAGAGCCATCAAAATCACAGTGGATGGGCCCCGAGAACCTCGAA|ATAAACCCCACTTGAAAAACTGAGGTGCTTAAGGAGTAAAATAATATGTTCCTGGTGGCATCCTCCAGATCGTACTGAGAAGCACTCCACAATGCCAGACTCACCTGTGGATGTGAAGACGCAATCTAGGCTGACTCCTCCAACAATGCCACCTCCCCCAACTACTCAAGGAGCTCCAAGAACCAGTTCATTTACACCGACAACGT___TAACTAATGGCACGAGCCATTCTCCTACAGCCTTGAATGGCGCCCCC	out-of-frame	FVGRSGRGKSFTLTITVFTNPPQVATYHRAIKITVDGPREPR|nkphlkn*	.
FP671120.1	IKZF1	+/+	+/+	21:8216497	7:50319060	intron	5'UTR	translocation	downstream	upstream	1	0	6	19	1343	medium	.	.	.	CGCCGCAGGTGCAGATCTTGGTGGTAGTAGCAAATATTCAAACGAGAACTTTGAAGGCCGAAGTGGAGAAGGGTTCCATGTGAACAGCAGTTGAACATGGGTCAGTCGGTCCTGAGAGATGGGCG...GGAGTCGGGTTCAGATCCCCGAATCCGGAGTGGCGGAGATGGGCGCCGCGAGGCGTCCAGTGCGGTAACGCGACCGATCCCGGAGAAGCCGGCGGGAGCCCCGGGGA|CCATGGATGCTGATGAGGGTCAAGACATGTCCCAAGTTTCAG___GGAAGGAAAGCCCCCCTGTAAGCGATACTCCAGATGAGGGC	.	.	.
FP236383.3	NAMPT	-/+	-/-	21:8438438	7:106276455	exon	exon	translocation/3'-3'	downstream	downstream	0	1	6	59	39414	medium	.	.	.	CTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACG|ACTGGAAACTCCTTTTGAATTATGGTAGCAATCAGAATATTTTTATATTAGCCAGTTTTACCTTGAAGACCTATTTTTAAAAACTACCTGTGTCTCTGGACTTAGTTGCAAATGCATATTAAAACAAAAATCCCCCAATTTCTGTGCTTTCTTATTTGAAAGGCCATTTCTAGGGGGAAAACAGTTCCCAAAC	.	.	.
MN1	ETV6	-/-	+/+	22:27796763	12:11826402	splice-site	splice-site	translocation	upstream	upstream	0	3	0	12	2176	medium	.	.	.	CGgTGGCGCCCTGGGAGAAGGCCAAACCCCAGAACCCCAACAGCAAAGAAG|CACAGAATGAAAAGCACCTTTGAAGAATCGTAAGCCATGAATGATGACGGGCCGAAGGGAGTGAGGGCTCCTACCTCTTGTCACTCCTACCATCCTGCCTCCAGGACCCCTGCCACCTACTTGGGAGTCTGCCCTTAAAACTCAAGGGTCACCTGACTCATTCCTCCCCACTCATTCATCTGACTACCCAC	out-of-frame	vAPWEKAKPQNPNSKE|aqnekhl*	.
NACC2	NOTCH1	-/-	-/-	9:136049636	9:136519565	splice-site	splice-site	duplication	upstream	downstream	1	2	0	2	1575	medium	.	.	.	CGACAGCCCCACCTCGTACCACAACGAGGAGGACGAGGAGGACGACGAGGCCTACGACACCATGGTGGAGGAGCAGTACGGCCAGATGTACATCAAGGCCTCCGGCAGCTATGCAG|GCTTCACCGGCCAGAACTGTGAGGAAAATATCGACGATTGTCCAGGAAACAACTGCAAGAACGGGGGTGCCTGTGTGGACGGCGTGAACACCTACAACTGCCGCTGCCCGCCAGAGTGGACAG___GTCAGTACTGTACCGAGGATGTGGACGAGTGCCAGCTGATGC	in-frame	DSPTSYHNEEDEEDDEAYDTMVEEQYGQMYIKASGSYA|gFTGQNCEENIDDCPGNNCKNGGACVDGVNTYNCRCPPEWTGQYCTEDVDECQLM	.
HOOK3	KAT6A	+/+	-/-	8:42897188	8:42049302	splice-site	splice-site	inversion	downstream	downstream	1	1	0	1	5175	medium	.	.	.	GGATCCCCCGACAGAGCGGCGGCGGTGTCTGGCCAGGCGGTAGGCGCTGCCTGGCCGCGGCGGGGAAGATGTTCAGCGTAGAGTCGCTGGAGCGGGCGGAGCTGTGCGAGAGCCTCCTCACTTGG|GCATAAGATGCACATTTTTCTGCTCTGGAGCCGGGAATGAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTGGAAGACATTTCAAaCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGaGTGG	out-of-frame	MFSVESLERAELCESLLTW|a*	.
HOOK3	KAT6A	+/+	-/-	8:42897188	8:42049707	splice-site	splice-site	inversion	downstream	downstream	0	3	0	1	1615	low	.	.	.	GAGCGGCGGCGGTGTCTGGCCAGGCGGTAGGCGCTGCCTGGCCGCGGCGGGGAAGATGTTCAGCGTAGAGTCGCTGGAGCGGGCGGAGCTGTGCGAGAGCCTCCTCACTTGG|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGTCAAGGAAGAGTTGAAATTCCTAGACTG___GCATAAGATGCACATTT...TGAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTGGAAGACATTTCAATCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGGGT	out-of-frame	MFSVESLERAELCESLLTW|dsfy*	.
HOOK3	KAT6A	+/+	-/-	8:42906258	8:42049707	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	1615	low	.	.	.	CCCAGGTTCTTCAAAAGAT|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGTCAAGGAAGAGTTGAAATTCCTAGACTGGTGAGCAGAGATTGCGGAGTACAGAATAGCCATCACCCACTCTG___GCATAAGATGCACATTTTTCTGCTCTGGAGCCGGGAATGAAATATTCTTGAGTTCTTACAACTTTATGACGAG	out-of-frame	QVLQKm|ilstnpdtc*	.
EVI2B	FP671120.1	-/-	+/+	17:31304253	21:8215461	3'UTR	intron	translocation	upstream	upstream	0	0	17	45697	94	low	.	.	duplicates(3)	CAATTGTATCTCCTCTTCCAAATGATTCTACTAGTCTCCCTCCATCTCTGGACTGTCTCAATCAAGACTGTGGAGATCATAAATCTGAGATAATACAATCATTTCCACCGCTTGACTCACTTAACTTGCCCCTGCCACCAGTAGATTTTATGAAAAACCAAGAAGATTCCAACCTTGAGATCCAGTGTCAGGAGTTCTCTATTCCTCCCAACTCTGATCAAGATCTTAATGAATCCCTGCCACCTCCACCTGCAGAACTGTTATAAATATTACAAC...|...CCGTAGCGGTCCTGACGTGCAAATCGGTCGTCCGACCTGGGTATAGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTCCCTCCGAAGTTTCCCTCAGGATAGCTGGCGCTCTCGCAGACCCGACGCACCCCCGCCACGCAGTTTTATCCGGTAAAGCGAATGATTAGAGGTCTTGGGGCCGAAACGATCTCAACCTATTCTCAAACTTTAAATGGGTAAGAAGCCCGGCTCGCTGGCGTGGAGCCGGGCGTGGAATGCGAGTGCCTAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGCGGGATGAACCGAACGCCGG	.	.	.
PICALM	NAMPT	-/-	-/-	11:85958469	7:106275454	3'UTR	exon	translocation	upstream	downstream	0	0	9	5535	65535	low	.	.	.	TATCAAGAAAAGCCCTTATgTGGAAAGTGTGAAATTTGTATTTGGAAAAGCTGCCTGGAGAGAAGAACTGTGTCCTTTAaTGTATTTCAACAGGACTCTTTTGGGG...CACAAATACAGAAACAATAACTGAAATTAACTTTTCTTTTTTTAAAAAAAcTTATATTCAGTTTGCAGTAGACATTtCTTAAGTATTTGTATTTATTTATGATTATCAATTTTACATAACATTAA...|...GGAAATTACTGACTTACCAACCCTAAAGATGCCACATTCTTTACGAACGAAGTGAGACGTTTCATTTGGATGTAAGAAGCTTTTGAGTCATCTCAAAGTTGAAAAGTACCAAAATTAATCGAGAG...AAGTTAATATCATTTTTAAAAGCTGTGCAAGTAATAATTTTTTCTATTACACATCTGCTTTTACAATAAATAAACATTTATACTTTGAGTTATTAAAAGTACAAATTCCTTTTCTGACTCTGCTG	.	.	.
NAMPT	PICALM	-/-	-/-	7:106250498	11:86031505	3'UTR	CDS	translocation	upstream	downstream	0	0	7	46180	13362	low	.	.	.	AGAAAATACTTGAAAATAAATTACATTGCCTTTTGTCCATTAATCAGCAAATAAAACATGGCCTTAACAAAGTTGTTTGTGTTATTGTACAATTTGAAAATTATGTCGGGcCATACCCTATAGAA...|...CATTACAACTCATCATTTGATGGTGTATGGAAATGAG___CGTTTTATTCAaTATTTGGCTTCAAGAAACACGTTGTTTAACTTAAGCAATTTTTTGGATAAAAGTGGATTGCAAG___GATATGACATGTCTACATTTATTAGGC	.	.	.
EVI2B	NAMPT	-/-	-/-	17:31304449	7:106284844	CDS	CDS	translocation	upstream	downstream	0	0	9	65535	24459	low	.	.	mismappers(2)	CTTAGAAATTAAGTTGTTTGAATCAAGTGAAAACATTGAAGACTCCAACAACCCCAAAACAGAGAAAATAAAAGATCAAGTAAATGGTACATCAGAAGATAGTGCTGATGGTTCAACAGTTGGAA...GAAGGACAGGAAAGTAACCAATCTGACAAACCCACAATGACAATTGTATCTCCTCTTCCAAATGATTCTACTAGTCTCCCTCCATCTCTGGACTGTCTCAATCAAGACTGTGGAGATCAT...|...CCACCGACTCCTACAAG___GTTACTCACTATAAACAATATCCACCCAACACAAGCAAAGTTTATTCCTACTTTGAATGCCGTGAAAAGAAGACAGAAAACTCCAAATTAAGGAAGGTGAAATATGAG	.	.	.
FP671120.1	EP300	+/+	+/+	21:8215554	22:41131426	intron	CDS	translocation	downstream	upstream	0	0	7	81	4997	low	.	.	.	CCGAGGGCGCACCACCGGCCCGTCTCGCCCGCCGCGCCGGGGAGGTGGAGCACGAGCGCACGTGTTAGGACCCGAAAGATGGTGAACTATGCCTGGGCAGGGCGAAGCCAGAGGAAACTCTGGTGGAGGTCaGTAGCGGTCCTGACGTGCAAATCGGTCGTCCGACCTGGGTATAGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTCCCTCCGAAGTT...|...AGCTCTCTAGGGGTGGGTCAACAGTCTGCCCCCAACCTAAGCACTGTTAGTCAGATTGATCCCAGCTCCATAGAAAGAGCCTATGCAGCTCTTGGACTACCCTATCAAGTAAATCAGATGCCGAC	.	.	.
NAMPT	MSI2	-/-	+/+	7:106249653	17:57264377	3'UTR	exon	translocation	upstream	upstream	0	0	7	17396	39	low	.	.	mismappers(1)	AACTATAGCCCAAAGGCCAAATgTGGACTTCTTTTTATAAATGCAGATTGTCTATGGCTGCTTTCCtACTACTCCAGCCTAAGGTAAACAGCTGCAATAaAAGCCAAATGAGAATC...|...TTTATTTATTTATTTATTTTTGAGATGGGATCTCACTGTGTTGCCCAGGCTGGTCgTAAACTCCTGGGCCAAGCTATCCTCCCACCTCTGCCCCCCAGAGTGCTGGGATTATAGGCGTGAGCCAC	.	.	.
NAMPT	YPEL5	-/-	+/+	7:106250346	2:30158935	3'UTR	3'UTR	translocation	upstream	upstream	0	0	7	30632	14164	low	.	.	.	TAATCATTCTACATTAAAGAAAATAATGGTTCTTACTGGAATGTCTAGGCACTGTACAGTTATTATATATCTTGGTTGTTGTATTGTACCAGTGAAATGCCAAATTTGAAAGGCCTGT...|...ATTAATGAACTGCGGAACAAGAGGTTGTGAGAATCTAAGATGGAACCTTTCTTTCTTTCTTTCTTTTTTTTTAAATTTTGTATTTTCCATCCAACAGCAGTGTGTAGAGAGAATATTATGCAGA...TGTTTACATCTTGAGGCAtCAtAGTCTGTCTGCAGCTATGTGGTGAGCTATGTAAGGAAAAAAATCTGGGCTGTTAGAGTGAAAAAGTGTGTTTTATGTCAATTGTGAAAGGAAAATGTTAGGAG	.	.	.
NAMPT	RNF213	-/-	+/+	7:106248751	17:80339245	3'UTR	CDS	translocation	upstream	upstream	0	0	7	8923	1567	low	.	.	.	GATTAATTCTGACTTTCTGaCACTAAGCACTTGGTATTTGGCCATCTCCATTCTGAGCACCAAACGGTTAACACGAATGTCCACTAGAACTCTGaTGTGTGgCACCCTTAAATCAGTCTAAATCT...|...GCcCTCTGCTGGGgATATGgTGTTCAGGACGTGGATCGCCATGGCCTACTGCTCCCCCAAGCAGGGTGTGTCCCTCCAAATGGACTTTGGCTTGGACCTGGTGACGGAGCTTAAAGAAGGTGGAG	.	.	.
NAMPT	CDK6	-/-	-/-	7:106250497	7:92614642	3'UTR	3'UTR	deletion	upstream	downstream	0	0	7	46180	905	low	.	.	.	AAAATACTTGAAAATAAATTACATTGCCTTTTGTCCATTAATCAGCAAATAAAACATGGCCTTAACAAAGaTGTTTGTGTTATTGTACAATTTGAAAATTATGTCGGGACATACCCTATAGAAT...|...AATTCAATTTTAAAGACTCAAGGTGGTCAGTAAATAACAGGCATTTGTTCACTGAAGGTGATTCACCAAAATAGTCTTCTCAAATTAGAAAGTTAACCCCATGTCCTCAGCATTTCTTTTCTGGC	.	.	.
EP300	RNF213	+/+	+/+	22:41178663	17:80313125	CDS	CDS	translocation	downstream	upstream	0	0	6	6229	2426	low	.	.	.	CCCCATGAGaCCCCAGCAGCATATGCTCCCAAATCAGaCCCAGTCCCCACACCTACAAGGaCAGCAGATCCCTAATTCTCTCTCCAATCAAGTGCGCTCTCCCCAGCCTGTCCCTTCTCCACGGC...|...CACATCTTCAGGTGCCTCATTCACATACGTCAAGGAAATTGAG___GTCTGGAGGCGGCTGGTGGAAATCCAATTCCCCGCGGAGCATGGCTGGAAGGAGTCGTTGCTGGGAGACATGGAATGGAaGC	.	.	.
CREBBP	FP671120.4	-/-	-/+	16:3726895	21:8254637	3'UTR	exon	translocation/5'-5'	upstream	upstream	0	0	6	8953	53	low	.	.	.	AAAAAATAGGGATGATCACTCTTAGACCATGCTAATGTTACTAGAGAAGAAGCCTTCTTTTCTTTCTTCTATGTGAAACTTGAAATGAGGAAAAGCAATTCTAGTGTAAATCATGCAAGCGCTCTAATTCCTATAAATACGAAACTCGAGAAGATTC...|...TTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGcCCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAG	.	.	.
CSF3R	FP671120.4	-/-	-/+	1:36466050	21:8254626	3'UTR	exon	translocation/5'-5'	upstream	upstream	0	0	6	563	53	low	.	.	.	CTTTAAGGACCAGATCATGCTCCATCCAGCCCCACCCAATGGCCTTTTGTGCTTGTTTCCTATAACTTCAGTATTGTAAACTAGTTTTTGGTTTGCAGTTTTTGTTGTTGTTTATAGACACTCTT...|...GGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAAC	.	.	.
FP671120.1	KAT6A	+/+	-/-	21:8216699	8:42048929	intron	CDS	translocation	downstream	downstream	0	0	6	31	6417	low	.	.	.	CGTCCAGTGCGGTAACGCGACCGATCCCGGAGAAGCCGGCGGGAGCCCCGGGGAGAGTTCTCTTTTCTTTGTGAAGGGCAGGGCGCCCTGGAATGGGTTCGCCCCGAGAGAGGGGCCCGTGCCTT...CGTCGCGGTTCCGGCGGCGTCCGGTGAGCTCTCGCTGGCCCTTGAAAATCCGGGGGAGAGGGTGTAAATCTCGCGCCGGGCCGTACCCATATCCGCAGCAGGTCTCCAAGGTGAACAGCaTCTGG...|...ATCAAAAAAGTGAAAAAGCAGAAACAGCGTCCTTCAGAAGAAAGGATATGCAATGCTGTGTCTTCATCCCATGGCTTGGATCGTAAAACTGTTTTAGAACAATTGGAGTTGAGTGTTAAAGATGG...TCAATTCCTATAAAGATCCTGATAATCCTGGGCGAATAGCACTTCCTAAGCCTCGGAACCATGGAAAATTGGATAATAAACAAAATGTGGATTGGAATAAACTGATAAAGCGGGCAGTTGAGGG	.	.	.
SF3B1	FP671120.1	-/-	+/+	2:197396201	21:8215662	CDS	intron	translocation	upstream	upstream	0	0	6	11510	85	low	.	.	duplicates(3)	CTCTTAAAAGCCCACAAAAAGGCTATTCGTAGAGCCACAGTCAACACATTTGGTTATATTGCAAAGtaCATTGG___CCCTCATGATGTATTGGCTACACTTCTGAACAACCTCAAAGTTCAAGAAAGGCAGAACgGAGTTTGTACCACTGTAGCAATAGCTATTGTTGCAGAAACATGTTCACCCTTTACAGTACTaCCTGCCT...|...ATTCTCAAACTTTAAATGGGTAAGAAGCCCGGCTCGCTGGCGTGGAGCCGGGCGTGGAATGCGAGTGCCTAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGCGGGATGAACCGAACGCCGGG...GGCCATGGAAGTCGGAATCCGCTAAGGAGTGTGTAACAACTCACCTGCCGAATCAACTAGCCCTGAAAATGGATGGCGCTGGAGCGTCGGGCCCATACCCGGCCGTCGCCGGCAGTCGAGAGTGG	.	.	.
FP236383.1	RNF213	+/+	+/+	21:8399276	17:80298346	intron	CDS	translocation	downstream	upstream	0	0	5	19	1834	low	.	.	.	CCCGCGGACGCTACGCCGCGACGAGTAGGAGGGCCGCTGCGGTGAGCCTTGAAGCCTAGGGCGCGGGCCCGGGTGGAGCCGCCGCAGGTGCAGATCTcGGTGGTAGTcGCA...|...TGCCAAAGATGCATGGACACAAGGACGTACACCTGGCTGGGCGCCCTGCCTGTCCTGCACTGCTGTATGGAGCTGGCCCCGCGGCACAAGGATGCCTGGAGACAGCCTGAGGACACCTGGGCCGC	.	.	.
ZMYM2	CDK6	+/+	-/-	13:20086011	7:92612665	CDS	3'UTR	translocation	downstream	downstream	0	0	3	2391	487	low	.	.	.	TCTACAGATAGCCCTGTCTGGTATACGTCTACTTCACTGGACCGAAACACCTTGGAAAATATGCTTGTACGGGTTCTTCTAGTAAAAGATATTTATGATAAAGACAATTATGAACTGGATGAAGACACAGAC...|...GGCAGTGCTTATGAAGAAAAATTAAAGCACAAACATTCTGGCATTCAATCGTTGGCAGATTATCTTCTGATGACACAGAATGAAAGGGCATCTCAGCCTCTCTGAACTTTGTAAAAATCTGTCCC	.	.	.
PCM1	TRIP11	+/+	-/-	8:17957355	14:92005718	CDS	CDS	translocation	downstream	downstream	0	0	3	2629	2952	low	.	.	.	TTCTGAGCATGAAAATTCCGAGCCTGTTACTAACATTCG___AAATCCACAAGTAGCTTCCACTTGGAATGAAGTAAATAGTCATAGTAATGCACAGTGTGTTTCTAATAATAGAGATGGGCGAACAGTTAATT...|...ATACCTCTaCCTTACAGCTGGAACATGAGCATTTAATTAAACTCAATCAAAAGAAAGACATGGAAATAGCAGAACTCAAAAAGAATATTGAACAAATGGATACTGACCATAAAGAAACTAAGGAC	.	.	.
KDM5A	CDK6	-/-	-/-	12:292955	7:92614582	CDS	3'UTR	translocation	upstream	downstream	0	0	3	6245	1298	low	.	.	.	CTTTTTGGAGAAGGAAAACAGAAGTCCAAGGAGTTAAAGAAAATGGACAAACCTAGAAAGAAGAAATTAAAATTAGGTGCAGACAAATCAAAGGAGCTGAATAAACTGGCCAAGAAACTAGaAAA...|...ATTCACCAAAATAGTCTTCTCAAATTAGAAAGTTAACCCCATGTCCTCAGCATTTCTTTTCTGGCCAAAAGCAGTAAATTTGCTAGCAGTAAAAGATGAAGTTTTATACACACAGCAAAAAGGAG	.	.	.
NAMPT	NUP98	-/-	-/-	7:106261588	11:3713995	splice-site	splice-site	translocation	upstream	downstream	1	0	2	51206	3553	low	.	.	mismappers(1)	CACAGGCACCACTAATAATCAGACCTGATTCTGGAAACCCTCTTGACACTGTGTTAAAG___GTTTTGGAGATTTTAGGTAAGAAGTTTCCTGTTACTGAGAACTCAAAGGGTTACAAGTTGCTGCCACCTTATCTTAGAGTTATTCAAGGGGATGGAGTAGATATTAATACCTTACAAGAG|GAAGGCTGAAGTTACATTGGATGGAGTTTGGCCAACAGATAAAACATCTCGTTGTTTAATAAAGAGCC	out-of-frame	QAPLIIRPDSGNPLDTVLKVLEILGKKFPVTENSKGYKLLPPYLRVIQGDGVDINTLQE|eg*	.
TTL	MSI2	+/+	+/+	2:112529260	17:57257517	3'UTR	CDS	translocation	downstream	upstream	0	0	2	252	1126	low	.	.	.	CTGCTCTGGACTATGGATTGGACGTCAGAGCATATTGGAGGTTGCCTGTGTGTTCCCCACCCcTCCCTTCGGTAACACTCTGCCACACTAAGCTCTGT...|...TCATGAGAGATCCCACTACGAAACGCTCCAG___AGGCTTCGGTTTCGTCACGTTCGCAGACCCAGCAAGTGTAGATAAAGTATTAGGTCAGCCCCACCATGAGTTAGATTCCAAGACG___ATTGACCCC	.	.	.
LYN	KAT6A	+/+	-/-	8:55999549	8:42049707	splice-site	splice-site	inversion	downstream	downstream	0	1	1	0	1615	low	.	.	.	CCTATGGGAAAATTCCCTACCCAG|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGTCAAGGAAGAGTTGAAATTCCTAGACTGGTGAGCAG...AAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTGGAAGACATTTCAATCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGGG	out-of-frame	YGKIPYP|gfflliqilvevltsflgkdvkeelkfldw*	.
IKBKB	KAT6A	+/+	-/-	8:42272205	8:42049707	splice-site	splice-site	inversion	downstream	downstream	0	2	0	0	1615	low	.	.	.	GCTGGTCACCTTCCCTGACAACGCAGACATGTGGGGCCTGGGAAATGAAAGAGCGCCTTGGGACAGGGGGATTTGGAAATGTCATCCGATGGCACAATCAG|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGTCA...GAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTGGAAGACATTTCAATCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGGGTGGCATTCTGTTTTGTGAAAGGAG	out-of-frame	MSSDGTIr|ilstnpdtc*	.
NAMPT	IKZF1	-/-	+/+	7:106284828	7:50319048	splice-site	splice-site	inversion	upstream	upstream	0	1	1	16050	1343	low	.	.	.	CGCAGCCGCGCCAGGGAGCTCGCGGCGCGCGGCCCCTGTaCTCCGGCCCGAGATGAATCCgGCGGCAGAAGCCGAGTTCAACATCCTCCTGGCCACCGACTCCTACAAG|ATAACCTGAGGACCATGGATGCTGATGAGGGTCAAGACATGTCCCAAGTTTCAG___GGAAGGAAAGCCCCCCTGTAAGCGATACTCCAGATGAGGGCGA...CCACCACCTCGGGAGGACAGCAAAGCTCCAAGAGTGACAGAGTCGTGG___CCAGTAATGTTAAAGTAGAGACTCAGAGTGATGAAGAGAATGGGCGTGCCTGTGAAATGAATGGGGAAGAATGgGCG	out-of-frame	MNPAAEAEFNILLATDSYK|it*	.
WNK1	ERC1	+/+	+/+	12:830160	12:1408149	splice-site	splice-site	deletion	downstream	upstream	0	1	0	0	1095	low	.	.	.	CCGTCGCGTGACCAGT|ACGCAAAATCGAATGAAGCTAATGGCCGACAACTACGAGGATGACCACTTCAAATCCTCCCATTCCAATCAAACAAATCACAAGCCCTCCCCAGACCAG	in-frame	RRVTS|TQNRMKLMADNYEDDHFKSSHSNQTNHKPSPD	.
WNK1	ERC1	+/+	+/+	12:754324	12:1027750	splice-site	splice-site	deletion/read-through	downstream	upstream	1	0	0	1	1074	low	.	.	.	GCGAGAAGTGGCAGCGGCGGCGGCAGCGCCAAGGAGCCACAGGAGGAACGGAGCCAGCAGCAGGATGATATCGAAGAGCTGGAGACCAAGGCCGTGGGAATGTCTAACGATGGCCGCTTTCTCAA...CTGAAACCACCGTGGAAGTCGCCTGGTGTGAACTGCAG|ATGGTGTAAGATACTTCTTCAAGATTGGACAGCTGGGGACCTTCTTCTGATTAACCTTAAACCAACTTGTAGCCATAGAGACACC	out-of-frame	ETTVEVAWCELQ|mv*	.
STAU2	KAT6A	-/-	-/-	8:73746783	8:42049302	splice-site	splice-site	deletion	upstream	downstream	0	1	0	0	5175	low	.	.	.	GTTTGCCAATGTTGGAGCCGTCTGCAAAGTGTCCCCGGCAAGAAG|GCATAAGATGCACATTTTTCTGCTCTGGAGCCGGGAATGAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGgGAAATTCATTGGGAAGTTGGAAGACATTTCAAaCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGGGTG	.	.	.
STAU2	KAT6A	-/-	-/-	8:73746783	8:42049707	splice-site	splice-site	deletion	upstream	downstream	0	1	0	0	1615	low	.	.	.	GCCGTCTGCAAAGTGTCCCCGGCAAGAAG|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGTCAAGGAAGAGTTGAAATTCCTAGACTGGTGAGCAGAGATTGCGGAGTACAGAATAGCC	.	.	.
NUP98	NSD1	-/-	+/+	11:3744509	5:177235821	splice-site	splice-site	translocation	upstream	upstream	1	0	0	851	982	low	.	.	.	CCTGTTTGGCAATAACAAGCTTACTACATTTGGAAGCAGCACAACCAGTGCACCTTCATTTGGTACAACCAGTGGCGGGCTCTTTG___GTTTTGGCACAAATACCAGTGGGAATAGTATTTTTGGAA...TTGGGGCCCCTGGATTTAATACTACGgCAGCCAgTTTGGGCTTTGGAGCCCCCCAGGCCCCAGTAG|CTGTGCGGTCAGAGAAGAAACGCCTTAGGAAGCCAAGCAAGTGGCTTTTGGAATATAC	in-frame	GAPGFNTTaAsLGFGAPQAPV|aVRSEKKRLRKPSKWLLEY	.
JAK2	RFX3	+/+	-/-	9:5126446	9:3293258	splice-site	splice-site	inversion	downstream	downstream	0	1	0	2449	1	low	.	.	.	ATTTGATAGAACTTTTGAAGAATAATGGAAGATTACCAAGACCAGATGGATGCCCAGATGAG|CTCCAGTGGCTGTTGGACAATTATGAGACAGCAGAAGGAGTGAGCCTTCCCAGAAGCACTCTGTACAACCACTACCTTCGACACTGTCAGGAACACAAACTGGACCCAGTCAATGCTGCCTCTTTTGGAAAATTAATAAGATCAATTTTTATGGGGCTACGAACCAGGAGATTGGGCA	in-frame	LIELLKNNGRLPRPDGCPDE|LQWLLDNYETAEGVSLPRSTLYNHYLRHCQEHKLDPVNAASFGKLIRSIFMGLRTRRLG	.
NCOA2	KAT6A	-/-	-/-	8:70296744	8:42049302	splice-site	splice-site	deletion	upstream	downstream	0	1	0	1087	5175	low	.	.	.	TATAGATCTTCTGCACTGTTTACAGGCACAG|GCATAAGATGCACATTTTTCTGCTCTGGAGCCGGGAATGAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTaGAAGACATTTCAAaCAACAGGT	.	.	.
NCOA2	KAT6A	-/-	-/-	8:70403700	8:42049707	splice-site	splice-site	deletion	upstream	downstream	0	1	0	701	1615	low	.	.	.	GCTTCGGCGGCGAAGGTCAGCGCCGACGGCAGCCGGCACCTGACGGCGTGACCGACCCGAGCCG|GATTCTTTCTACTAATCCAGATACTTGTTGAAGTGCTGACTAGTTTCTTGGGGAAAGATGT...GAAATATTCTTGAGTTCTTACAACTTTATGACGAGACCCATGTGTGGTGCTATTGAGAAATTCATTGGGAAGTTGGAAGACATTTCAAaCAACAGGTTGTTTTGGTTTCTATAGTACAATTGGGG	.	.	.
LINC00921	CREBBP	+/+	-/-	16:3263883	16:3851009	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	4923	low	.	.	.	AGCCACGAAGTGGTGGGAAATACGAGATTCAAATATCTGGTGGCCGTGTGAC|ATTTTGGATCATTGTTTGACTTGGAAAATGATCTTCCTGATGAGCTGATACCCAATGGAGGAGAATTAGGCCTTTTAAACAGTGGGAACCTTGTTCCAGATGCTGCTTCCAAACATAAACAACTGTCGGAGC	.	.	.
LINC00921	CREBBP	+/+	-/-	16:3263883	16:3731937	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	10946	low	.	.	.	ACGAGATTCAAATATCTGGTGGCCGTGTGAC|GGCAGTCAGGGCGACAGCAAGAATGCCAAGAAGAAGAACAACAAGAAAACCAACAAGAACAAAAGCAGCATCAGCCGCGCCAACAAGAAGAAGCCCAGCATGCCCAACGTGTCCAATGACCTGTCCCAGAAGCTGTATGCCACCATGGAGAAGCACAAGGAG___GTCTTCTTCGTGATCCACCTGCACGCTGGGCCTGTCATCAACACCCTGCC	.	.	.
ZNF263	CREBBP	+/+	-/-	16:3263883	16:3731937	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	10946	low	.	.	.	ACGAGATTCAAATATCTGGTGGCCGTGTGAC|GGCAGTCAGGGCGACAGCAAGAATGCCAAGAAGAAGAACAACAAGAAAACCAACAAGAACAAAAGCAGCATCAGCCGCGCCAACAAGAAGAAGCCCAGCATGCCCAACGTGTCCAATGACCTGTCCCAGAAGCTGTATGCCACCATGGAGAAGCACAAGGAG___GTCTTCTTCGTGATCCACCTGCACGCTGGGCCTGTCATCAACACCCTGCC	.	.	.
ZNF263	CREBBP	+/+	-/-	16:3263883	16:3851009	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	4923	low	.	.	.	AGCCACGAAGTGGTGGGAAATACGAGATTCAAATATCTGGTGGCCGTGTGAC|ATTTTGGATCATTGTTTGACTTGGAAAATGATCTTCCTGATGAGCTGATACCCAATGGAGGAGAATTAGGCCTTTTAAACAGTGGGAACCTTGTTCCAGATGCTGCTTCCAAACATAAACAACTGTCGGAGC	.	.	.
CREBBP	CORO7-PAM16	-/-	-/-	16:3736036	16:4365545	splice-site	splice-site	duplication	upstream	downstream	1	0	0	9417	0	low	.	.	.	TACAAAAAGATGCTGGtCAAGGCGTTTGCAGAGCGGATCATCCATGACTACAAG___GATATTTTCAAACAAGCAACTGAAGACAGGCTCACCAGTGCCAAGGAACTGCCCTATTTTGAAGGTGATTT...GAAGAGAGCACTGCAGCCAGTGAAACCACTGAG|GTGTCTCGTGCCTCTGCTGGACCCTGACTCTGGGCTCCTGGTCCTGGCAGGAAAGGTGAGTGAGAAGCTGGGGCTCCCACCCACATCCACA	out-of-frame	EESTAASETTE|vsrasagp*	.
CREBBP	CORO7-PAM16	-/-	-/-	16:3767720	16:4405567	splice-site	splice-site	duplication	upstream	downstream	0	1	0	10192	1	low	.	.	.	AAGAAACCTGAAGTGAAAGTAGAAGTTAAAGAGGAAGAAGAGAGTAGCAGTAACGGCACAGCCTCTCAGTCAACATCTCCTTCGCAGCCGCGCAAAAAAA|AGCTGGCAGCCCATGGGGACCTGGTGCAGAGCGCCGTCTGGAGCCGAGATGGAGCCCTGGTGGGCACGGCGTGCAAG___GACAAGCAGCTGCaGATCTTTGACCCCAGAACAAAGCCGCGGGCC	in-frame	KKPEVKVEVKEEEESSSNGTASQSTSPSQPRKK|kLAAHGDLVQSAVWSRDGALVGTACKDKQLqIFDPRTKPR	.
CREBBP	CORO7	-/-	-/-	16:3736036	16:4365545	splice-site	splice-site	duplication	upstream	downstream	1	0	0	9417	0	low	.	.	.	TACAAAAAGATGCTGGtCAAGGCGTTTGCAGAGCGGATCATCCATGACTACAAG___GATATTTTCAAACAAGCAACTGAAGACAGGCTCACCAGTGCCAAGGAACTGCCCTATTTTGAAGGTGATTT...GAAGAGAGCACTGCAGCCAGTGAAACCACTGAG|GTGTCTCGTGCCTCTGCTGGACCCTGACTCTGGGCTCCTGGTCCTGGCAGGAAAGGTGAGTGAGAAGCTGGGGCTCCCACCCACATCCACA	out-of-frame	EESTAASETTE|vsrasagp*	.
CREBBP	CORO7	-/-	-/-	16:3767720	16:4405567	splice-site	splice-site	duplication	upstream	downstream	0	1	0	10192	1	low	.	.	.	AAGAAACCTGAAGTGAAAGTAGAAGTTAAAGAGGAAGAAGAGAGTAGCAGTAACGGCACAGCCTCTCAGTCAACATCTCCTTCGCAGCCGCGCAAAAAAA|AGCTGGCAGCCCATGGGGACCTGGTGCAGAGCGCCGTCTGGAGCCGAGATGGAGCCCTGGTGGGCACGGCGTGCAAG___GACAAGCAGCTGCaGATCTTTGACCCCAGAACAAAGCCGCGGGCC	in-frame	KKPEVKVEVKEEEESSSNGTASQSTSPSQPRKK|kLAAHGDLVQSAVWSRDGALVGTACKDKQLqIFDPRTKPR	.
PRKCB	CREBBP	+/+	-/-	16:23837406	16:3851009	splice-site	splice-site	inversion	downstream	downstream	0	1	0	0	4923	low	.	.	.	GGGGCTTCGGGAAGCAGGGATTCCAGTGCCAAG|ATTTTGGATCATTGTTTGACTTGGAAAATGATCTTCCTGATGAGCTGATACCCAATGGAGGAGAATTAGGCCTTTTAAACAGTGGGAACCTTGTTCCAGATGCTGCTTCCAAACATAAACAACTGTCGGAGC	.	.	.
CREBBP	AC005520.1	-/-	+/+	16:3736036	14:73892027	splice-site	splice-site	translocation	upstream	upstream	1	0	0	9417	0	low	.	.	.	ATTTTCAAACAAGCAACTGAAGACAGGCTCACCAGTGCCAAGGAACTGCCCTATTTTGAAGGTGATTTCTGGCCCAATGTGTTAGAAGAGAGCATTAAGGAACTAGAACAAGAAGAAGAGGAGAGGAAAAAGGAAGAGAGCACTGCAGCCAGTGAAACCACTGAG|GTTACATTGATTACCCACCTAGTACAACATCTTACGGGAAG	out-of-frame	IFKQATEDRLTSAKELPYFEGDFWPNVLEESIKELEQEEEERKKEESTAASETTE|vtlithlvqhltg	.
CREBBP	ZNF410	-/-	+/+	16:3736036	14:73892027	splice-site	splice-site	translocation	upstream	upstream	1	0	0	9417	0	low	.	.	.	ATTTTCAAACAAGCAACTGAAGACAGGCTCACCAGTGCCAAGGAACTGCCCTATTTTGAAGGTGATTTCTGGCCCAATGTGTTAGAAGAGAGCATTAAGGAACTAGAACAAGAAGAAGAGGAGAGGAAAAAGGAAGAGAGCACTGCAGCCAGTGAAACCACTGAG|GTTACATTGATTACCCACCTAGTACAACATCTTACGGGAAG	out-of-frame	IFKQATEDRLTSAKELPYFEGDFWPNVLEESIKELEQEEEERKKEESTAASETTE|vtlithlvqhltg	.
BCL6	EIF4A2	-/-	+/+	3:187736090	3:186785883	splice-site	splice-site	inversion	upstream	upstream	0	1	0	3265	4526	low	.	.	.	CACTAACCTTGGAGCCGATGGGATTGAGTGACTGGCACTTGGGACCACAGAGAAATGTCAGAGTGTTTGGTTACAGACTCAAGGAAACCTCTCATTTTAGAGTGCTCATTTG|ATCCAAAAGGTAATTCTGGCACTTGGAGACTATATGGGAGCCACTTGTCATGCCTGCATTGGTGGAACAAATGTTCGAAATGAAATGCAAAAACTGCAGGCTGAAGCACCACATATTGTTGTTGGTACACCCGGG	.	.	.
NAMPT	MCM3AP	-/-	-/-	7:106284828	21:46285394	splice-site	splice-site	translocation	upstream	downstream	0	1	0	16050	920	low	.	.	.	CTGCGGCAGAAGCCGAGTTCAACATCCTCCTGGCCACCGACTCCTACAAG|CTTGTTCTAAAGACGAAGGGGAACTCCCTAGTGCTACAGGCTGTTTCAGGAAGCCAATTAATTTTGTACACATAATACTTAATTACCTTCTAATAATTGGAGCAGAAGATGAACCCAACTAATCCTTTCAGTGGGCAGCAGCCTAGTGCTTTTTCGGCGTCTTCTAGTAATGTAGGAACACTTCCATCTAAGCCGCC	out-of-frame	AAEAEFNILLATDSYK|lvlktkgnslvlqavsgsqlilyt*	.
PDE8A	DMD	+/.	-/.	15:85126150	X:31161503	intron	intron	translocation	upstream	upstream	5	1	0	494	83	medium	.	.	duplicates(4)	.	.	.	."""
        with open(self.tmp_in_tsv, "w") as writer:
            writer.write(content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_tsv, self.tmp_out_tsv]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_ReadWrite(self):
        # Read, convert to VCF, convert to dict and write
        with ArribaIO(self.tmp_out_tsv, "w") as writer:
            with ArribaIO(self.tmp_in_tsv) as reader:
                for first_bnd, second_bnd in reader:
                    writer.write((first_bnd, second_bnd))
        # Eval
        expected = None
        with open(self.tmp_out_tsv) as expec_reader:
            expected = expec_reader.readlines()
        observed = None
        with open(self.tmp_out_tsv) as obs_reader:
            observed = obs_reader.readlines()
        self.assertEqual(expected, observed)


class BreakendVCFIOTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp = os.path.join(tmp_folder, unique_id + ".vcf")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGetSub(self):
        content = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
4	888888	.	T	G	.	PASS	.
13	123456	bnd_U	C	C[2 : 321682[,C[17:198983[	6	PASS	SVTYPE=BND;MATEID=bnd_V,bnd_Z
17	198983	bnd_Z	A	]13:123456]A	6	PASS	SVTYPE=BND;MATEID=bnd_U'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        try:
            pysam.tabix_compress(self.tmp, self.tmp + ".gz")
            pysam.tabix_index(self.tmp + ".gz", preset="vcf")
            with BreakendVCFIO(self.tmp + ".gz", "i") as reader:
                # Missing chromosome
                observed = [elt.getName() for elt in reader.getSub("1", 1437, 11437)]
                self.assertEqual(observed, [])
                # No variants
                observed = [elt.getName() for elt in reader.getSub("2", 1234580, 1234590)]
                self.assertEqual(observed, [])
                # Fusion 2 mates
                observed = [elt for elt in reader.getSub("13", 123450, 123460)]
                for idx, curr in enumerate(observed):
                    observed[idx] = (
                        curr[0].getName(),
                        [elt.getName() for elt in curr[1]]
                    )
                self.assertEqual(
                    observed,
                    [(
                        "13:123456=C/C[2 : 321682[/C[17:198983[",
                        ["2:321682=T/]13 : 123456]T", "17:198983=A/]13:123456]A"]
                    )]
                )
        finally:
            for curr_file in [self.tmp + ".gz", self.tmp + ".gz.tbi"]:
                if os.path.exists(curr_file):
                    os.remove(curr_file)

    def test_read(self):
        # Create input tmp
        content = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
4	888888	.	T	G	.	PASS	.
13	123456	bnd_U	C	C[2 : 321682[,C[17:198983[	6	PASS	SVTYPE=BND;MATEID=bnd_V,bnd_Z
17	198983	bnd_Z	A	]13:123456]A	6	PASS	SVTYPE=BND;MATEID=bnd_U'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        # Expected
        expected = [
            "bnd_V\tbnd_U_0\t]13 : 123456]T\t@@\tbnd_U_0\tbnd_V\tC[2 : 321682[",
            "bnd_U_1\tbnd_Z\tC[17:198983[\t@@\tbnd_Z\tbnd_U_1\t]13:123456]A"
        ]
        # Observed
        observed = []
        with BreakendVCFIO(self.tmp) as reader:
            for first, second in reader:
                first_str = "{}\t{}\t{}".format(first.id, first.info["MATEID"][0], first.alt[0])
                second_str = "{}\t{}\t{}".format(second.id, second.info["MATEID"][0], second.alt[0])
                observed.append("{}\t@@\t{}".format(first_str, second_str))
        # Eval
        self.assertEqual(expected, observed)

    def test_read2(self):
        # Create input tmp
        content = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
4	888888	.	T	G	.	PASS	.
13	123456	bnd_U	C	C[2 : 321682[	6	PASS	MATEID=bnd_V;SVTYPE=BND
18	888888	.	T	G	.	PASS	.
20	888888	.	T	G	.	PASS	.'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        # Expected
        expected = [
            "bnd_V\tbnd_U\t]13 : 123456]T\t@@\tbnd_U\tbnd_V\tC[2 : 321682["
        ]
        # Observed
        observed = []
        with BreakendVCFIO(self.tmp) as reader:
            for first, second in reader:
                first_str = "{}\t{}\t{}".format(first.id, first.info["MATEID"][0], first.alt[0])
                second_str = "{}\t{}\t{}".format(second.id, second.info["MATEID"][0], second.alt[0])
                observed.append("{}\t@@\t{}".format(first_str, second_str))
        # Eval
        self.assertEqual(expected, observed)

    def test_isValid(self):
        # Create input tmp
        content = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
4	888888	.	T	G	.	PASS	.
13	123456	bnd_U	C	C[2 : 321682[	6	PASS	MATEID=bnd_V;SVTYPE=BND
18	888888	.	T	G	.	PASS	.
20	888888	.	T	G	.	PASS	.'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        self.assertTrue(BreakendVCFIO.isValid(self.tmp))
        # Create input tmp
        content = '''##fileformat=VCFv4.3
##FILTER=<ID=HLA,Description="One breakend is located on HLA.">
##FILTER=<ID=IG,Description="One breakend is located on immunoglobulin.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
##FILTER=<ID=Readthrough,Description="The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= 100000).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
20	888888	.	T	G	.	PASS	.'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        self.assertTrue(not BreakendVCFIO.isValid(self.tmp))

    def test_write(self):
        data = [(
            VCFRecord("2", 321682, "bnd_V", "T", ["]13 : 123456]T"], "6", ["PASS"], {"SVTYPE": "BND", "MATEID": ["bnd_U"]}),
            VCFRecord("13", 123456, "bnd_U", "C", ["C[2:321682["], "6", ["PASS"], {"SVTYPE": "BND", "MATEID": ["bnd_V"]})
        )]
        with BreakendVCFIO(self.tmp, "w") as writer:
            writer.info = {
                "MATEID": HeaderInfoAttr(id="MATEID", number="A", type="String", description="ID of mate breakend."),
                "SVTYPE": HeaderInfoAttr(id="SVTYPE", number="1", type="String", description="Type of structural variant."),
            }
            writer.filter = {"MATEID": HeaderFilterAttr(id="Inner", description="The two breakends are located in the same gene.")}
            writer.writeHeader()
            for first, second in data:
                writer.write(first, second)
        # Expected
        expected = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	MATEID=bnd_U;SVTYPE=BND
13	123456	bnd_U	C	C[2:321682[	6	PASS	MATEID=bnd_V;SVTYPE=BND'''
        # Observed
        observed = None
        with open(self.tmp) as reader:
            observed = "".join(reader.readlines()).strip()
        # Eval
        self.assertEqual(expected, observed)

    def test_readWithAnnot(self):
        # Create input tmp
        content = '''##fileformat=VCFv4.3
##INFO=<ID=TESTANN,Number=.,Type=String,Description="Annotations. Format: SYMBOL|Gene|Feature|Feature_type|STRAND">
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_V	T	]13 : 123456]T	6	PASS	TESTANN=SPECC1L|ENSG00000258555|ENST00000358654|Transcript|+;SVTYPE=BND;MATEID=bnd_U
4	888888	.	T	G	.	PASS	.
13	123456	bnd_U	C	C[2 : 321682[	6	PASS	TESTANN=NTRK2|ENSG00000148053|ENST00000376213|Transcript|-;MATEID=bnd_V;SVTYPE=BND
18	888888	.	T	G	.	PASS	.
20	888888	.	T	G	.	PASS	.'''
        with open(self.tmp, "w") as writer:
            writer.write(content)
        # Expected
        expected = [
            "SPECC1L_NTRK2\t+_-"
        ]
        # Observed
        observed = []
        with BreakendVCFIO(self.tmp, "r", "TESTANN") as reader:
            for first, second in reader:
                observed.append(
                    "{}_{}\t{}_{}".format(
                        first.info["TESTANN"][0]["SYMBOL"],
                        second.info["TESTANN"][0]["SYMBOL"],
                        first.info["TESTANN"][0]["STRAND"],
                        second.info["TESTANN"][0]["STRAND"]
                    )
                )
        # Eval
        self.assertEqual(expected, observed)

    def test_writeWithAnnot(self):
        data = [
            (
                VCFRecord(
                    "2", 321682, "bnd_A", "T", ["]13 : 123456]T"], "6", ["PASS"],
                    {
                        "MATEID": ["bnd_C"],
                        "SVTYPE": "BND",
                        "TESTANN": [{
                            "SYMBOL": "SPECC1L",
                            "STRAND": "+"
                        }]
                    }
                ),
                VCFRecord(
                    "13", 123456, "bnd_C", "C", ["C[2 : 321682["], "6", ["PASS"],
                    {
                        "MATEID": ["bnd_A"],
                        "SVTYPE": "BND",
                        "TESTANN": [{
                            "SYMBOL": "NTRK2",
                            "STRAND": "-"
                        }]
                    }
                )
            ),
            (
                VCFRecord(
                    "8", 323698, "bnd_B", "T", ["]15:654321]T"], "30", ["PASS"],
                    {
                        "MATEID": ["bnd_D"],
                        "SVTYPE": "BND"
                    }
                ),
                VCFRecord(
                    "15", 654321, "bnd_D", "C", ["C[8:323698["], "30", ["PASS"],
                    {
                        "MATEID": ["bnd_B"],
                        "SVTYPE": "BND",
                        "TESTANN": [
                            {"SYMBOL": "GENF", "STRAND": "-"},
                            {"SYMBOL": "GENC", "STRAND": "+"}
                        ]
                    }
                )
            )
        ]
        # Expected
        expected = '''##fileformat=VCFv4.3
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=TESTANN,Number=.,Type=String,Description="Annotations. Format: SYMBOL|STRAND">
##FILTER=<ID=Inner,Description="The two breakends are located in the same gene.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	321682	bnd_A	T	]13 : 123456]T	6	PASS	MATEID=bnd_C;SVTYPE=BND;TESTANN=SPECC1L|+
13	123456	bnd_C	C	C[2 : 321682[	6	PASS	MATEID=bnd_A;SVTYPE=BND;TESTANN=NTRK2|-
8	323698	bnd_B	T	]15:654321]T	30	PASS	MATEID=bnd_D;SVTYPE=BND
15	654321	bnd_D	C	C[8:323698[	30	PASS	MATEID=bnd_B;SVTYPE=BND;TESTANN=GENF|-,GENC|+'''
        # Observed
        with BreakendVCFIO(self.tmp, "w", "TESTANN") as writer:
            writer.info = {
                "MATEID": HeaderInfoAttr(id="MATEID", number="A", type="String", description="ID of mate breakend."),
                "SVTYPE": HeaderInfoAttr(id="SVTYPE", number="1", type="String", description="Type of structural variant."),
                "TESTANN": HeaderInfoAttr(id="TESTANN", number=".", type="String", description="Annotations. Format: SYMBOL|STRAND")
            }
            writer.filter = {"MATEID": HeaderFilterAttr(id="Inner", description="The two breakends are located in the same gene.")}
            writer.ANN_titles = ["SYMBOL", "STRAND"]
            writer.writeHeader()
            for first, second in data:
                writer.write(first, second)
        observed = None
        with open(self.tmp) as reader:
            observed = "".join(reader.readlines()).strip()
        # Eval
        self.assertEqual(expected, observed)


class UtilsTest(unittest.TestCase):
    def test_getCoordDictFromCoordStr(self):
        # Ok
        data = [
            {
                "input": "chr1:10:+",
                "expected": {
                    "chrom": "chr1",
                    "pos": 10,
                    "strand": "+"
                }
            },
            {
                "input": "1:222526778:-",
                "expected": {
                    "chrom": "1",
                    "pos": 222526778,
                    "strand": "-"
                }
            }
        ]
        expected = []
        observed = []
        for curr in data:
            expected.append(curr["expected"])
            observed.append(getCoordDictFromCoordStr(curr["input"]))
        self.assertEqual(expected, observed)
        # Invalid field
        with self.assertRaises(Exception) as context:
            getCoordDictFromCoordStr("chr1:a:+")
        with self.assertRaises(Exception) as context:
            getCoordDictFromCoordStr("chr1:20:a")

    def test_getAltFromCoord(self):
        data = [
            # First shard is + and second is -
            {
                "first_coord": {
                    "chrom": "chr2",
                    "pos": 321681,
                    "strand": "+"
                },
                "second_coord": {
                    "chrom": "chr17",
                    "pos": 198982,
                    "strand": "-"
                },
                "expected": ("N]chr17:198982]", "N]chr2:321681]")
            },
            # First shard is + and second is +
            {
                "first_coord": {
                    "chrom": "chr13",
                    "pos": 123456,
                    "strand": "+"
                },
                "second_coord": {
                    "chrom": "chr2",
                    "pos": 321682,
                    "strand": "+"
                },
                "expected": ("N[chr2:321682[", "]chr13:123456]N")
            },
            # First shard is - and second is -
            {
                "first_coord": {
                    "chrom": "chr9",
                    "pos": 88969186,
                    "strand": "-"
                },
                "second_coord": {
                    "chrom": "chr9",
                    "pos": 102784508,
                    "strand": "-"
                },
                "expected": ("]chr9:102784508]N", "N[chr9:88969186[")
            },
            # First shard is - and second is +
            {
                "first_coord": {
                    "chrom": "chr13",
                    "pos": 123457,
                    "strand": "-"
                },
                "second_coord": {
                    "chrom": "chr17",
                    "pos": 198983,
                    "strand": "+"
                },
                "expected": ("[chr17:198983[N", "[chr13:123457[N")
            }
        ]
        expected = []
        observed = []
        for curr in data:
            expected.append(curr["expected"])
            observed.append(getAltFromCoord(curr["first_coord"], curr["second_coord"]))
        self.assertEqual(expected, observed)

    def test_getCoordStr(self):
        data = [
            # First shard is + and second is -
            {
                "in": VCFRecord(
                    region="chr2",
                    position=321681,
                    knownSNPId="1",
                    refAllele="N",  # ref
                    altAlleles=["N]chr17:198982]"],  # alt
                    info={"RNA_FIRST": True, "MATEID": "2"}
                ),
                "expected": {
                    "chrom": "chr2",
                    "pos": 321681,
                    "strand": "+"
                },
            },
            {
                "in": VCFRecord(
                    region="chr17",
                    position=198982,
                    knownSNPId="1",
                    refAllele="N",  # ref
                    altAlleles=["N]chr2:321681]"],  # alt
                    info={"MATEID": "1"}
                ),
                "expected": {
                    "chrom": "chr17",
                    "pos": 198982,
                    "strand": "-"
                }
            },
            # First shard is + and second is +
            {
                "in": VCFRecord(
                    region="chr13",
                    position=123456,
                    knownSNPId="3",
                    refAllele="N",  # ref
                    altAlleles=["N[chr2:321682["],  # alt
                    info={"RNA_FIRST": True, "MATEID": "4"}
                ),
                "expected": {
                    "chrom": "chr13",
                    "pos": 123456,
                    "strand": "+"
                }
            },
            {
                "in": VCFRecord(
                    region="chr2",
                    position=321682,
                    knownSNPId="4",
                    refAllele="N",  # ref
                    altAlleles=["]chr13:123456]N"],  # alt
                    info={"MATEID": "3"}
                ),
                "expected": {
                    "chrom": "chr2",
                    "pos": 321682,
                    "strand": "+"
                }
            },
            # First shard is - and second is -
            {
                "in": VCFRecord(
                    region="chr9",
                    position=88969186,
                    knownSNPId="5",
                    refAllele="N",  # ref
                    altAlleles=["]chr9:102784508]N"],  # alt
                    info={"RNA_FIRST": True, "MATEID": "6"}
                ),
                "expected": {
                    "chrom": "chr9",
                    "pos": 88969186,
                    "strand": "-"
                }
            },
            {
                "in": VCFRecord(
                    region="chr9",
                    position=102784508,
                    knownSNPId="6",
                    refAllele="N",  # ref
                    altAlleles=["N[chr9:88969186["],  # alt
                    info={"MATEID": "5"}
                ),
                "expected": {
                    "chrom": "chr9",
                    "pos": 102784508,
                    "strand": "-"
                }
            },
            # First shard is - and second is +
            {
                "in": VCFRecord(
                    region="chr13",
                    position=123457,
                    knownSNPId="7",
                    refAllele="N",  # ref
                    altAlleles=["[chr17:198983[N"],  # alt
                    info={"RNA_FIRST": True, "MATEID": "8"}
                ),
                "expected": {
                    "chrom": "chr13",
                    "pos": 123457,
                    "strand": "-"
                }
            },
            {
                "in": VCFRecord(
                    region="chr17",
                    position=198983,
                    knownSNPId="8",
                    refAllele="N",  # ref
                    altAlleles=["[chr13:123457[N"],  # alt
                    info={"MATEID": "7"}
                ),
                "expected": {
                    "chrom": "chr17",
                    "pos": 198983,
                    "strand": "+"
                }
            }
        ]
        expected = []
        observed = []
        for curr in data:
            expected.append(curr["expected"])
            observed.append(getCoordStr(curr["in"]))
        self.assertEqual(expected, observed)


#####################################################################
#
# MAIN
#
#####################################################################
if __name__ == '__main__' :
    unittest.main()
