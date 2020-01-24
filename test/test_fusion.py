#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.vcf import VCFRecord
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
                "expected": STARFusionIO,
                "content": '#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots'
            },
            {
                "expected": ArribaIO,
                "content": 'gene1	gene2	strand1	strand2	breakpoint1	breakpoint2	site1	site2	type	direction1	direction2	split_reads1	split_reads2	discordant_mates	confidence'
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
            with open(self.tmp_in, "w") as writer:
                writer.write(curr_test["content"])
            with FusionFileReader.factory(self.tmp_in, "r", "splA", annot_field="AnnotField") as reader:
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
        self.tmp_in_tsv = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out_tsv = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Create in TSV
        content = '''#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots
MN1--BEND2	90	39	ONLY_REF_SPLICE	MN1^ENSG00000169184.6	chr22:27796763:-	BEND2^ENSG00000177324.14	chrX:18195442:-	M70265:74:000000000-B4F4J:1:2116:12294:7575,M70265:74:000000000-B4F4J:1:2115:7184:4292,M70265:74:000000000-B4F4J:1:1104:5050:21345,M70265:74:000000000-B4F4J:1:1119:9539:18627,M70265:74:000000000-B4F4J:1:2115:8678:9709,M70265:74:000000000-B4F4J:1:2119:19014:22293,M70265:74:000000000-B4F4J:1:2105:16442:9076,M70265:74:000000000-B4F4J:1:1117:12801:7721,M70265:74:000000000-B4F4J:1:1111:14401:18304,M70265:74:000000000-B4F4J:1:1112:22240:14651,M70265:74:000000000-B4F4J:1:1108:20469:22260,M70265:74:000000000-B4F4J:1:2102:11105:8918,M70265:74:000000000-B4F4J:1:2108:19409:16466,M70265:74:000000000-B4F4J:1:1113:20196:5091,M70265:74:000000000-B4F4J:1:1105:10766:13643,M70265:74:000000000-B4F4J:1:1103:19101:13070,M70265:74:000000000-B4F4J:1:2107:24460:10546,M70265:74:000000000-B4F4J:1:1106:19866:14418,M70265:74:000000000-B4F4J:1:1106:5053:13766,M70265:74:000000000-B4F4J:1:2113:20643:17000,M70265:74:000000000-B4F4J:1:2102:25631:11285,M70265:74:000000000-B4F4J:1:1119:27699:19981,M70265:74:000000000-B4F4J:1:2113:4962:18403,M70265:74:000000000-B4F4J:1:2106:25372:20356,M70265:74:000000000-B4F4J:1:2108:16489:16372,M70265:74:000000000-B4F4J:1:1105:14999:8683,M70265:74:000000000-B4F4J:1:1107:13181:5229,M70265:74:000000000-B4F4J:1:1103:9694:18262,M70265:74:000000000-B4F4J:1:1106:7893:19397,M70265:74:000000000-B4F4J:1:1119:28664:16235,M70265:74:000000000-B4F4J:1:1119:25900:13031,M70265:74:000000000-B4F4J:1:2117:5341:12350,M70265:74:000000000-B4F4J:1:2108:21016:10605,M70265:74:000000000-B4F4J:1:1111:6184:20952,M70265:74:000000000-B4F4J:1:1111:7115:5525,M70265:74:000000000-B4F4J:1:1110:27378:16651,M70265:74:000000000-B4F4J:1:1108:14118:21000,M70265:74:000000000-B4F4J:1:2110:8046:7787,M70265:74:000000000-B4F4J:1:2112:29162:12358,M70265:74:000000000-B4F4J:1:1114:12762:7330,M70265:74:000000000-B4F4J:1:1112:6478:12309,M70265:74:000000000-B4F4J:1:2114:8353:4989,M70265:74:000000000-B4F4J:1:1110:12165:19870,M70265:74:000000000-B4F4J:1:2119:26398:14585,M70265:74:000000000-B4F4J:1:1104:19297:14091,M70265:74:000000000-B4F4J:1:1114:12488:10493,M70265:74:000000000-B4F4J:1:2110:19477:14905,M70265:74:000000000-B4F4J:1:1111:29396:14060,M70265:74:000000000-B4F4J:1:1112:5997:8955,M70265:74:000000000-B4F4J:1:1114:21233:24893,M70265:74:000000000-B4F4J:1:1117:16546:12881,M70265:74:000000000-B4F4J:1:2102:15683:16355,M70265:74:000000000-B4F4J:1:1102:7568:7922,M70265:74:000000000-B4F4J:1:1104:24497:10850,M70265:74:000000000-B4F4J:1:1112:12330:15917,M70265:74:000000000-B4F4J:1:2104:23488:20289,M70265:74:000000000-B4F4J:1:2115:18539:5639,M70265:74:000000000-B4F4J:1:2108:27004:11103,M70265:74:000000000-B4F4J:1:1108:7280:6580,M70265:74:000000000-B4F4J:1:1111:13093:9434,M70265:74:000000000-B4F4J:1:1105:22026:13215,M70265:74:000000000-B4F4J:1:1103:12233:24724,M70265:74:000000000-B4F4J:1:1114:21833:4508,M70265:74:000000000-B4F4J:1:1102:5182:21387,M70265:74:000000000-B4F4J:1:2102:20234:8226,M70265:74:000000000-B4F4J:1:1104:5911:7048,M70265:74:000000000-B4F4J:1:1104:19911:7039,M70265:74:000000000-B4F4J:1:2101:23116:15940,M70265:74:000000000-B4F4J:1:1102:28272:19048,M70265:74:000000000-B4F4J:1:2107:3513:10128,M70265:74:000000000-B4F4J:1:2101:3654:8461,M70265:74:000000000-B4F4J:1:1119:19110:18013,M70265:74:000000000-B4F4J:1:1107:18400:1636,M70265:74:000000000-B4F4J:1:2115:15123:16288,M70265:74:000000000-B4F4J:1:2110:13483:19384,M70265:74:000000000-B4F4J:1:2101:9347:11364,M70265:74:000000000-B4F4J:1:1103:14624:5192,M70265:74:000000000-B4F4J:1:2103:25588:18903,M70265:74:000000000-B4F4J:1:1117:15805:5600,M70265:74:000000000-B4F4J:1:1119:13495:19891,M70265:74:000000000-B4F4J:1:1108:7984:5492,M70265:74:000000000-B4F4J:1:1117:9253:13094,M70265:74:000000000-B4F4J:1:1112:8428:11232,M70265:74:000000000-B4F4J:1:2103:22447:6340,M70265:74:000000000-B4F4J:1:1111:26219:5015,M70265:74:000000000-B4F4J:1:1112:14939:4411,M70265:74:000000000-B4F4J:1:2109:18888:3219,M70265:74:000000000-B4F4J:1:2110:22338:16190,M70265:74:000000000-B4F4J:1:2114:11325:2713,M70265:74:000000000-B4F4J:1:2114:6567:13905	M70265:74:000000000-B4F4J:1:1111:7154:15877,M70265:74:000000000-B4F4J:1:2117:22762:16468,M70265:74:000000000-B4F4J:1:1112:3545:11649,M70265:74:000000000-B4F4J:1:1102:13385:24356,M70265:74:000000000-B4F4J:1:1103:21865:16107,M70265:74:000000000-B4F4J:1:2109:19125:3230,M70265:74:000000000-B4F4J:1:1113:11062:16439,M70265:74:000000000-B4F4J:1:2115:11207:19666,M70265:74:000000000-B4F4J:1:2102:1922:14639,M70265:74:000000000-B4F4J:1:2105:29283:14469,M70265:74:000000000-B4F4J:1:1109:15243:9627,M70265:74:000000000-B4F4J:1:1108:22583:13263,M70265:74:000000000-B4F4J:1:1103:23453:4405,M70265:74:000000000-B4F4J:1:1113:9694:20456,M70265:74:000000000-B4F4J:1:2103:19632:8012,M70265:74:000000000-B4F4J:1:1102:18799:15205,M70265:74:000000000-B4F4J:1:1102:19075:11592,M70265:74:000000000-B4F4J:1:2114:28568:18618,M70265:74:000000000-B4F4J:1:1112:21951:11679,M70265:74:000000000-B4F4J:1:1111:22886:13062,M70265:74:000000000-B4F4J:1:1115:17822:21329,M70265:74:000000000-B4F4J:1:1116:21890:11280,M70265:74:000000000-B4F4J:1:1104:10237:9204,M70265:74:000000000-B4F4J:1:1112:3494:8607,M70265:74:000000000-B4F4J:1:2103:5353:4657,M70265:74:000000000-B4F4J:1:2101:18273:12696,M70265:74:000000000-B4F4J:1:2109:13674:17280,M70265:74:000000000-B4F4J:1:1103:18734:3105,M70265:74:000000000-B4F4J:1:2102:17275:5848,M70265:74:000000000-B4F4J:1:1112:24048:22285,M70265:74:000000000-B4F4J:1:1115:21982:9046,M70265:74:000000000-B4F4J:1:1109:27605:16863,M70265:74:000000000-B4F4J:1:1111:14542:18181,M70265:74:000000000-B4F4J:1:1104:9119:13899,M70265:74:000000000-B4F4J:1:1114:10485:22262,M70265:74:000000000-B4F4J:1:2111:24444:4566,M70265:74:000000000-B4F4J:1:2113:24284:23644,M70265:74:000000000-B4F4J:1:2104:27141:8778,M70265:74:000000000-B4F4J:1:2117:13953:8079	YES_LDAS	38.9246	GT	1.4566	AG	1.7968	["HaasMedCancer","INTERCHROMOSOMAL[chr22--chrX]"]
SPECC1L-ADORA2A--NTRK2	22	6	ONLY_REF_SPLICE	SPECC1L-ADORA2A^ENSG00000258555.6	chr22:24334573:+	NTRK2^ENSG00000148053.16	chr9:84867243:+	M70265:74:000000000-B4F4J:1:1102:11424:4669,M70265:74:000000000-B4F4J:1:1104:6405:2923,M70265:74:000000000-B4F4J:1:1114:21462:20506,M70265:74:000000000-B4F4J:1:1108:9784:2773,M70265:74:000000000-B4F4J:1:2119:22961:20082,M70265:74:000000000-B4F4J:1:2116:27522:11441,M70265:74:000000000-B4F4J:1:1108:23971:7984,M70265:74:000000000-B4F4J:1:1113:23518:22366,M70265:74:000000000-B4F4J:1:1115:7837:20035,M70265:74:000000000-B4F4J:1:2104:21327:2495,M70265:74:000000000-B4F4J:1:1101:9462:22886,M70265:74:000000000-B4F4J:1:1108:5888:7340,M70265:74:000000000-B4F4J:1:2103:6756:2599,M70265:74:000000000-B4F4J:1:1114:10446:14926,M70265:74:000000000-B4F4J:1:1108:4192:20832,M70265:74:000000000-B4F4J:1:2110:22090:10259,M70265:74:000000000-B4F4J:1:1107:22426:22529,M70265:74:000000000-B4F4J:1:1102:20172:13788,M70265:74:000000000-B4F4J:1:1102:14526:14874,M70265:74:000000000-B4F4J:1:1108:28140:12698,M70265:74:000000000-B4F4J:1:1109:25864:11274,M70265:74:000000000-B4F4J:1:1109:21996:25012	M70265:74:000000000-B4F4J:1:1108:25374:7116,M70265:74:000000000-B4F4J:1:1102:12684:14754,M70265:74:000000000-B4F4J:1:1102:12917:5068,M70265:74:000000000-B4F4J:1:1101:28832:12042,M70265:74:000000000-B4F4J:1:1101:26215:18199,M70265:74:000000000-B4F4J:1:2110:18851:11245	YES_LDAS	8.3732	GT	1.9656	AG	1.8295	["INTERCHROMOSOMAL[chr22--chr9]"]
SPECC1L-ADORA2A--NTRK2	3	6	ONLY_REF_SPLICE	SPECC1L-ADORA2A^ENSG00000258555.6	chr22:24334573:+	NTRK2^ENSG00000148053.16	chr9:84861040:+	M70265:74:000000000-B4F4J:1:1108:15540:9098,M70265:74:000000000-B4F4J:1:2104:21069:2777,M70265:74:000000000-B4F4J:1:1116:20816:4555	M70265:74:000000000-B4F4J:1:1101:28832:12042,M70265:74:000000000-B4F4J:1:1102:12917:5068,M70265:74:000000000-B4F4J:1:1101:26215:18199,M70265:74:000000000-B4F4J:1:1108:25374:7116,M70265:74:000000000-B4F4J:1:1102:12684:14754,M70265:74:000000000-B4F4J:1:2110:18851:11245	YES_LDAS	2.6914	GT	1.9656	AG	1.7232	["INTERCHROMOSOMAL[chr22--chr9]"]'''
        with open(self.tmp_in_tsv, "w") as writer:
            writer.write(content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_tsv, self.tmp_out_tsv]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_ReadWrite(self):
        # Read, convert to VCF, convert to dict and write
        with STARFusionIO(self.tmp_out_tsv, "w") as writer:
            with STARFusionIO(self.tmp_in_tsv) as reader:
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


class ArribaIOTest(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        # Temporary files
        self.tmp_in_tsv = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out_tsv = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Create in TSV
        content = """gene1	gene2	strand1	strand2	breakpoint1	breakpoint2	site1	site2	type	direction1	direction2	split_reads1	split_reads2	discordant_mates	confidence
FLI1	EWSR1	+/+	+/+	11:128675327	22:29688126	splice-site	splice-site	translocation	downstream	upstream	0	4	0	high
EWSR1	FLI1	+/+	+/+	22:29684776	11:128679052	splice-site	splice-site	translocation	downstream	upstream	1	0	2	high
S100A6	EEF1A1	-/-	-/-	1:153507149	6:74229778	3'UTR	5'UTR	translocation	upstream	downstream	0	1	11	high
EEF1A1	EEF2	-/-	-/-	6:74227963	19:3980028	exon	exon	translocation	upstream	downstream	2	0	10	high
AHNAK	HLA-B	-/-	-/-	11:62286163	6:31322973	exon	exon	translocation	upstream	downstream	1	0	10	high
B2M	HLA-B	+/+	-/-	15:45010186	6:31324076	3'UTR	exon	translocation	downstream	downstream	0	0	10	medium
KRT6A	LDHA	-/-	+/+	12:52881705	11:18418365	exon	5'UTR	translocation	upstream	upstream	0	0	10	medium
FLNA	KRT6A	-/-	-/-	X:153578196	12:52886822	5'UTR	exon	translocation	upstream	downstream	0	0	10	medium
GAPDH	MUC5B	+/+	+/+	12:6647133	11:1264669	exon	exon	translocation	downstream	upstream	0	0	10	medium
MUC5B	GAPDH	+/+	+/+	11:1282959	12:6646126	3'UTR	exon	translocation	downstream	upstream	0	0	9	low
GAPDH	ALDOA	+/+	+/+	12:6647376	16:30077176	3'UTR	5'UTR	translocation	downstream	upstream	0	0	10	medium
S100A6	GAPDH	-/-	+/+	1:153507153	12:6645674	3'UTR	UTR	translocation	upstream	upstream	0	0	10	medium
BCAS4	BCAS3	+/+	+/+	20:49411711	17:59445688	splice-site	splice-site	translocation	downstream	upstream	140	129	300	high
BCAS4	BCAS3	+/+	+/+	20:49411711	17:59430949	splice-site	splice-site	translocation	downstream	upstream	10	3	299	high
BCAS4	BCAS3	+/+	+/+	20:49411711	17:59431723	splice-site	splice-site	translocation	downstream	upstream	0	1	299	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22029432	splice-site	splice-site	deletion/read-through	downstream	upstream	3	3	20	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:21995047	splice-site	splice-site	deletion/read-through	downstream	upstream	1	0	22	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22092307	splice-site	splice-site	deletion/read-through	downstream	upstream	1	3	14	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22046750	splice-site	splice-site	deletion/read-through	downstream	upstream	0	1	16	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22097257	splice-site	splice-site	deletion/read-through	downstream	upstream	4	1	10	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22096371	splice-site	splice-site	deletion/read-through	downstream	upstream	2	1	12	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22112319	splice-site	splice-site	deletion/read-through	downstream	upstream	4	3	7	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22120503	splice-site	splice-site	deletion/read-through	downstream	upstream	5	1	4	high
MTAP	CDKN2B-AS1	+/+	+/+	9:21818202	9:22120199	splice-site	splice-site	deletion/read-through	downstream	upstream	1	0	4	high
MTAP	UBA52P6	+/+	+/+	9:21818202	9:22012127	splice-site	5'UTR	deletion/read-through	downstream	upstream	1	2	1	high
MTAP	RP11-408N14.1	+/+	-/+	9:21818202	9:22210663	splice-site	intron	deletion/read-through/5'-5'	downstream	upstream	3	4	3	medium"""
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
