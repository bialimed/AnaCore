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
import re
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FusionCatcherToVCF(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_input = os.path.join(tmp_folder, unique_id + "_fusions.tsv")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_fusions.vcf")

        # Exec command
        self.cmd = [
            "fusionCatcherToVCF.py",
            "--annotation-field", "TESTANN",
            "--sample-name", "testSample",
            "--input-fusions", self.tmp_input,
            "--output-fusions", self.tmp_output
        ]

        # Create input file
        content = """Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect
ANKRD10	MPHOSPH6	exon-exon	0	3	6	27	BOWTIE;BOWTIE+BLAT;BOWTIE+STAR	13:111558380:-	16:82197799:-	ENSG00000088448	ENSG00000135698	ENSE00003468677	ENSE00003484650	TAGAATGCATCAGTGCCCTTGTGGCGAATGGGGCTCACGTCGA*TTTATGCAAAGGGGACTGGACTCAGAAACCAAGAAACAACTAG	out-of-frame
HTN3	MSANTD3	known,m23,exon-exon	0	44	27	30	BOWTIE;BOWTIE+BLAT;BOWTIE+STAR	4:70894233:+	9:103204188:+	ENSG00000205649	ENSG00000066697	ENSE00002491917	ENSE00001465026	ACTTCAGCTTCACTGACTTCTGGATTCTCCTCTTGAGTAAAAG*GATAGCTAGCGGCCAGGAGAAATACAGTGGAAAATGCAAAACA	UTR/UTR
AC018630.6	AL139099.5	m2	0	336	6	42	BOWTIE+BLAT;BOWTIE+STAR	12:11034887:-	14:50053459:+	ENSG00000275778	ENSG00000283029			CCATCAGCAAGGTCCTCCCCCACCTCCTCCTGGAAAGCCCCAGGGACCAC*CAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGT	CDS(truncated)/exonic(no-known-CDS)
AL139099.5	AC018630.6	m2	0	197	6	23	BOWTIE+STAR	14:50053551:+	12:11294449:-	ENSG00000283029	ENSG00000275778			AGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAAT*AAGCCACTGCACTCCAGCCTGGGTGACAGAGCAAGACCCTGTCTCCCTCC	exonic(no-known-CDS)/intronic
KMT2C	EGFR	oncogene,cancer,tumor,exon-exon	0	9	14	29	BOWTIE	7:152132711:-	7:55209979:+	ENSG00000055609	ENSG00000146648	ENSE00001642231	ENSE00003541288	CTCGCAAAGATGGCGCTTCCCCTTTCCAGAGAGCCAGAAAGAA*TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTT	out-of-frame
MN1	BEND2	oncogene,cancer,tumor,exon-exon	0	40	19	30	BOWTIE;BOWTIE+BLAT;BOWTIE+STAR	22:28192751:-	X:18213562:-	ENSG00000169184	ENSG00000177324	ENSE00001150385	ENSE00001275482	CCCTGGGAGAAGGCCAAACCCCAGAACCCCAACAGCAAAGAAG*CTGAGAAAATTATACTTACTGAAATGCCAGGAACAACAGAAAC	in-frame"""
        with open(self.tmp_input, "w") as FH_out:
            FH_out.write(content)


    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_input, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """##fileformat=VCFv4.1
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend.">
##INFO=<ID=RNA_CONTIG,Number=1,Type=String,Description="The sequence of the breakend spanning contig.">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this break-end is 5' in the fusion transcript.">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Aligning method used for mapping the reads and finding the fusion genes. Here are two methods used which are: (i) BOWTIE = only Bowtie aligner is used for mapping the reads on the genome and exon-exon fusion junctions, (ii) BOWTIE+BLAT = Bowtie aligner is used for mapping reads on the genome and BLAT is used for mapping reads for finding the fusion junction, (iii) BOWTIE+STAR = Bowtie aligner is used for mapping reads on the genome and STAR is used for mapping reads for finding the fusion junction, (iv) BOWTIE+BOWTIE2 = Bowtie aligner is used for mapping reads on the genome and Bowtie2 is used for mapping reads for finding the fusion junction.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=TESTANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|Gene|EXON|Effect.">
##FORMAT=<ID=CM,Number=1,Type=Integer,Description="Count of reads mapping simultaneously on both genes which form the fusion gene. This is an indication how similar are the DNA/RNA sequences of the genes forming the fusion gene (i.e. what is their homology because highly homologous genes tend to appear show as candidate fusion genes). In case of completely different sequences of the genes involved in forming a fusion gene then here it is expected to have the value zero.">
##FORMAT=<ID=LA,Number=1,Type=Integer,Description="Longest anchor (hangover) found among the unique reads mapping on the fusion junction.">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Count of pairs of reads supporting the fusion (including also the multimapping reads).">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Count of unique reads (i.e. unique mapping positions) mapping on the fusion junction. Shortly, here are counted all the reads which map on fusion junction minus the PCR duplicated reads.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	testSample
13	111558380	e131255a-a545-45d7-99ec-76bbfe1559c1	N	]16:82197799]N	.	exon-exon	MATEID=68638193-2763-43ef-a615-950893d685d7;RNA_FIRST;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=ANKRD10|ENSG00000088448|ENSE00003468677|out-of-frame	CM:LA:PR:SR	0:27:3:6
16	82197799	68638193-2763-43ef-a615-950893d685d7	N	N[13:111558380[	.	exon-exon	MATEID=e131255a-a545-45d7-99ec-76bbfe1559c1;RNA_CONTIG=TAGAATGCATCAGTGCCCTTGTGGCGAATGGGGCTCACGTCGA*TTTATGCAAAGGGGACTGGACTCAGAAACCAAGAAACAACTAG;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=MPHOSPH6|ENSG00000135698|ENSE00003484650|out-of-frame	CM:LA:PR:SR	0:27:3:6
4	70894233	91b19b81-9191-43a9-8e0a-c56c90d37d4a	N	N[9:103204188[	.	known;m23;exon-exon	MATEID=c8f6c747-5c60-4b70-819e-8b9e4f5d5825;RNA_FIRST;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=HTN3|ENSG00000205649|ENSE00002491917|UTR/UTR	CM:LA:PR:SR	0:30:44:27
9	103204188	c8f6c747-5c60-4b70-819e-8b9e4f5d5825	N	]4:70894233]N	.	known;m23;exon-exon	MATEID=91b19b81-9191-43a9-8e0a-c56c90d37d4a;RNA_CONTIG=ACTTCAGCTTCACTGACTTCTGGATTCTCCTCTTGAGTAAAAG*GATAGCTAGCGGCCAGGAGAAATACAGTGGAAAATGCAAAACA;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=MSANTD3|ENSG00000066697|ENSE00001465026|UTR/UTR	CM:LA:PR:SR	0:30:44:27
12	11034887	c6800f9e-b5ac-4d82-8154-776fca7bd70e	N	[14:50053459[N	.	m2	MATEID=8c153d50-920c-42d7-82e2-7d3ee6c371c6;RNA_FIRST;SOURCES=BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=AC018630.6|ENSG00000275778||CDS(truncated)/exonic(no-known-CDS)	CM:LA:PR:SR	0:42:336:6
14	50053459	8c153d50-920c-42d7-82e2-7d3ee6c371c6	N	[12:11034887[N	.	m2	MATEID=c6800f9e-b5ac-4d82-8154-776fca7bd70e;RNA_CONTIG=CCATCAGCAAGGTCCTCCCCCACCTCCTCCTGGAAAGCCCCAGGGACCAC*CAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGT;SOURCES=BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=AL139099.5|ENSG00000283029||CDS(truncated)/exonic(no-known-CDS)	CM:LA:PR:SR	0:42:336:6
14	50053551	c33b318c-73a0-4b64-98b2-5346678cad99	N	N]12:11294449]	.	m2	MATEID=7cbb584d-5f88-4fc1-aa9d-54b9ae2d4b99;RNA_FIRST;SOURCES=BOWTIE+STAR;SVTYPE=BND;TESTANN=AL139099.5|ENSG00000283029||exonic(no-known-CDS)/intronic	CM:LA:PR:SR	0:23:197:6
12	11294449	7cbb584d-5f88-4fc1-aa9d-54b9ae2d4b99	N	N]14:50053551]	.	m2	MATEID=c33b318c-73a0-4b64-98b2-5346678cad99;RNA_CONTIG=AGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAAT*AAGCCACTGCACTCCAGCCTGGGTGACAGAGCAAGACCCTGTCTCCCTCC;SOURCES=BOWTIE+STAR;SVTYPE=BND;TESTANN=AC018630.6|ENSG00000275778||exonic(no-known-CDS)/intronic	CM:LA:PR:SR	0:23:197:6
7	152132711	1c7f1b75-010c-45a3-93a9-5ca5a432988a	N	[7:55209979[N	.	oncogene;cancer;tumor;exon-exon	MATEID=43304d3f-135a-4979-ba1d-8b76b06f3ab5;RNA_FIRST;SOURCES=BOWTIE;SVTYPE=BND;TESTANN=KMT2C|ENSG00000055609|ENSE00001642231|out-of-frame	CM:LA:PR:SR	0:29:9:14
7	55209979	43304d3f-135a-4979-ba1d-8b76b06f3ab5	N	[7:152132711[N	.	oncogene;cancer;tumor;exon-exon	MATEID=1c7f1b75-010c-45a3-93a9-5ca5a432988a;RNA_CONTIG=CTCGCAAAGATGGCGCTTCCCCTTTCCAGAGAGCCAGAAAGAA*TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTT;SOURCES=BOWTIE;SVTYPE=BND;TESTANN=EGFR|ENSG00000146648|ENSE00003541288|out-of-frame	CM:LA:PR:SR	0:29:9:14
22	28192751	03eb6c26-bb30-43ba-86c4-4db5f3619d85	N	]X:18213562]N	.	oncogene;cancer;tumor;exon-exon	MATEID=a7ee3d94-1a0c-4fc7-9ca4-e3c5ac4be1df;RNA_FIRST;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=MN1|ENSG00000169184|ENSE00001150385|in-frame	CM:LA:PR:SR	0:30:40:19
X	18213562	a7ee3d94-1a0c-4fc7-9ca4-e3c5ac4be1df	N	N[22:28192751[	.	oncogene;cancer;tumor;exon-exon	MATEID=03eb6c26-bb30-43ba-86c4-4db5f3619d85;RNA_CONTIG=CCCTGGGAGAAGGCCAAACCCCAGAACCCCAACAGCAAAGAAG*CTGAGAAAATTATACTTACTGAAATGCCAGGAACAACAGAAAC;SOURCES=BOWTIE,BOWTIE+BLAT,BOWTIE+STAR;SVTYPE=BND;TESTANN=BEND2|ENSG00000177324|ENSE00001275482|in-frame	CM:LA:PR:SR	0:30:40:19"""
        observed = None
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        # Remove unique element (IDs)
        expected = [re.sub("MATEID=[^;]+;", "", elt) for elt in expected.split("\n")]
        expected = [re.sub("\s[^\s]+-[^\s]+-[^\s]+-[^\s]+-[^\s]+\s", " id ", elt) for elt in expected]
        observed = [re.sub("MATEID=[^;]+;", "", elt) for elt in observed.split("\n")]
        observed = [re.sub("\s[^\s]+-[^\s]+-[^\s]+-[^\s]+-[^\s]+\s", " id ", elt) for elt in observed]
        # Compare
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
