#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.13.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sequenceIO import IdxFastaIO
import os
import pysam
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.vcf import VCFRecord, VCFIO, getHeaderAttr, HeaderInfoAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestVCFHeader(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_variants = os.path.join(tmp_folder, unique_id + "_in.vcf")
        self.tmp_out_variants = os.path.join(tmp_folder, unique_id + "_out.vcf")
        content = """##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##SAMPLE=<ID=NA00001,Assay="WholeGenome",Description="Patient germline genome from unaffected",Disease="None",DOI="url",Ethnicity="AFR">
##SAMPLE=<ID=NA00002,Assay="Exome",Description="European patient exome from breast cancer",Disease="Cancer",Ethnicity="CEU",Tissue="Breast">
##PEDIGREE=<ID=TumourSample,Original=GermlineID>
##PEDIGREE=<ID=SomaticNonTumour,Original=GermlineID>
##PEDIGREE=<ID=ChildID,Father=FatherID,Mother=MotherID>
##PEDIGREE=<ID=SampleID,Name_1=Ancestor_1,Name_N=Ancestor_N>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",Source="myImputationProgram",Version="3.1">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SRC,Number=.,Type=String,Description="Caller source">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership",Version="129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype \\"test\\"">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
"""
        with open(self.tmp_in_variants, "w") as FH_variants:
            FH_variants.write(content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_variants, self.tmp_out_variants]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testReadWriteHeader(self):
        with VCFIO(self.tmp_in_variants) as FH_in:
            with VCFIO(self.tmp_out_variants, "w") as FH_out:
                FH_out.copyHeader(FH_in)
                FH_out.writeHeader()
        expected = None
        with open(self.tmp_in_variants) as FH_exp:
            expected = FH_exp.readlines()
        observed = None
        with open(self.tmp_out_variants) as FH_obs:
            observed = FH_obs.readlines()
        self.assertTrue(len(expected) != 0)
        self.assertEqual(sorted(expected), sorted(observed))

    def testHeaderAttrDel(self):
        data = HeaderInfoAttr(
            id="NS",
            number="1",
            type="Integer",
            description="Number of Samples With Data",
            source="myImputationProgram",
            version="3.1"
        )
        del(data.source)
        expected = {"id": "NS", "number": "1", "type": "Integer", "description": "Number of Samples With Data", "version": "3.1"}
        self.assertEqual(expected, data.datastore)

    def testHeaderAttrKeys(self):
        data = HeaderInfoAttr(
            id="NS",
            number="1",
            type="Integer",
            description="Number of Samples With Data",
            source="myImputationProgram",
            version="3.1"
        )
        expected = sorted(["id", "number", "type", "description", "source", "version"])
        self.assertEqual(expected, sorted(data.keys()))

    def testHeaderAttrStr(self):
        data = HeaderInfoAttr(
            id="NS",
            number="1",
            type="Integer",
            description="Number of Samples With Data",
            source="myImputationProgram",
            version="3.1"
        )
        expected = '<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",Source="myImputationProgram",Version="3.1">'
        self.assertEqual(expected, str(data))

    def testGetHeaderAttr(self):
        data = [
            {
                "text": '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",Source="myImputationProgram",Version="3.1">',
                "expected_dict": {"id": "NS", "number": "1", "type": "Integer", "description": "Number of Samples With Data", "source": "myImputationProgram", "version": "3.1"},
                "expected_cls": "HeaderInfoAttr"
            },
            {
                "text": '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                "expected_dict": {"id": "DP", "number": "1", "type": "Integer", "description": "Total Depth"},
                "expected_cls": "HeaderInfoAttr"
            },
            {
                "text": '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                "expected_dict": {"id": "AF", "number": "A", "type": "Float", "description": "Allele Frequency"},
                "expected_cls": "HeaderInfoAttr"
            },
            {
                "text": '##INFO=<ID=SRC,Number=.,Type=String,Description="Caller source">',
                "expected_dict": {"id": "SRC", "number": ".", "type": "String", "description": "Caller source"},
                "expected_cls": "HeaderInfoAttr"
            },
            {
                "text": '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, 2018-01",Version="129">',
                "expected_dict": {"id": "DB", "number": "0", "type": "Flag", "description": "dbSNP membership, 2018-01", "version": "129"},
                "expected_cls": "HeaderInfoAttr"
            },
            {
                "text": '##FILTER=<ID=q10,Description="Quality below 10",Source="myImputationProgram",Version=10>',
                "expected_dict": {"id": "q10", "description": "Quality below 10", "source": "myImputationProgram", "version": "10"},
                "expected_cls": "HeaderFilterAttr"
            },
            {
                "text": '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
                "expected_dict": {"id": "s50", "description": "Less than 50% of samples have data"},
                "expected_cls": "HeaderFilterAttr"
            },
            {
                "text": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "expected_dict": {"id": "GT", "number": "1", "type": "String", "description": "Genotype"},
                "expected_cls": "HeaderFormatAttr"
            },
            {
                "text": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype \\"test\\"">',
                "expected_dict": {"id": "GT", "number": "1", "type": "String", "description": 'Genotype "test"'},
                "expected_cls": "HeaderFormatAttr"
            }
        ]
        for curr in data:
            attr = getHeaderAttr(curr["text"])
            self.assertEqual(curr["expected_dict"], attr.datastore)
            self.assertEqual(curr["expected_cls"], attr.__class__.__name__)


class TestNoneInfo(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.vcf")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.vcf")
        with open(self.tmp_in, "w") as FH_variants:
            FH_variants.write("""##fileformat=VCFv4.3
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	.	G	A	29.0	PASS	DP=.;AF=.	GT	0|0	1|0	1/1
20	14371	.	T	A	29.0	PASS	.	GT	0|0	1|0	1/1
20	14372	.	C	A	29.0	PASS	DP=300;AF=.	GT	0|0	1|0	1/1
20	14373	.	C	A	29.0	PASS	DP=300	GT	0|0	1|0	1/1
20	14374	.	G	A	29.0	PASS	DP=.;AF=0.5	GT	0|0	1|0	1/1
20	14375	.	T	A	29.0	PASS	AF=0.5	GT	0|0	1|0	1/1
20	14376	.	C	A	29.0	PASS	DP=300;AF=0.5	GT	0|0	1|0	1/1""")

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testRecord(self):
        with VCFIO(self.tmp_in) as reader:
            for record in reader:
                if record.pos in {14374, 14375, 14376}:
                    self.assertTrue("AF" in record.info and record.info["AF"] == [0.5])
                else:
                    self.assertTrue("AF" not in record.info)
                if record.pos in {14372, 14373, 14376}:
                    self.assertTrue("DP" in record.info and record.info["DP"] == 300)
                else:
                    self.assertTrue("DP" not in record.info)

    def testWrite(self):
        expected_content = """##fileformat=VCFv4.3
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	.	G	A	29.0	PASS	.	GT	0|0	1|0	1/1
20	14371	.	T	A	29.0	PASS	.	GT	0|0	1|0	1/1
20	14372	.	C	A	29.0	PASS	DP=300	GT	0|0	1|0	1/1
20	14373	.	C	A	29.0	PASS	DP=300	GT	0|0	1|0	1/1
20	14374	.	G	A	29.0	PASS	AF=0.5	GT	0|0	1|0	1/1
20	14375	.	T	A	29.0	PASS	AF=0.5	GT	0|0	1|0	1/1
20	14376	.	C	A	29.0	PASS	AF=0.5;DP=300	GT	0|0	1|0	1/1"""
        # Read and write VCF
        with VCFIO(self.tmp_in) as reader:
            with VCFIO(self.tmp_out, "w") as writer:
                writer.copyHeader(reader)
                writer.writeHeader()
                for record in reader:
                    writer.write(record)
        # Compare input content with output
        observed_content = None
        with open(self.tmp_out) as reader:
            observed_content = "".join(reader.readlines()).strip()
        self.assertEqual(observed_content, expected_content)


class TestVCFIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_variants = os.path.join(tmp_folder, unique_id + "_in.vcf")
        self.tmp_in_variants_with_spl_info = os.path.join(tmp_folder, unique_id + "_in_detailed.vcf")
        self.tmp_out_variants = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Create VCF
        spec_example = """##fileformat=VCFv4.3
##fileDate=20090805"
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4
20	1230237	.	T	-	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3"""
        with open(self.tmp_in_variants, "w") as FH_variants:
            FH_variants.write(spec_example)

        # Create VCF with spl info
        spec_example_with_spl_info = spec_example.replace(
            "#CHROM",
            """##SAMPLE=<ID=NA00001,Assay=WholeGenome,Ethnicity=AFR,Disease=None,Description="Patient germline genome from unaffected">
##SAMPLE=<ID=NA00003,Assay=Exome,Ethnicity=CEU,Description="European patient exome from breast cancer",Disease=Cancer,Tissue=Breast>
#CHROM""")
        with open(self.tmp_in_variants_with_spl_info, "w") as FH_variants:
            FH_variants.write(spec_example_with_spl_info)

    def testIter(self):
        with VCFIO(self.tmp_in_variants) as FH_vcf:
            # Header
            self.assertEqual(FH_vcf.samples, ["NA00001", "NA00002", "NA00003"])
            self.assertEqual(sorted(list(FH_vcf.format.keys())), sorted(["GT", "GQ", "DP", "HQ"]))
            self.assertEqual(sorted(list(FH_vcf.info.keys())), sorted(["NS", "DP", "AF", "AA", "DB", "H2"]))
            # Records
            expected_records = ["20_14370_G_A", "20_17330_T_A", "20_1110696_A_G,T", "20_1230237_T_-", "20_1234567_GTC_G,GTCT"]
            readed_records = list()
            for variant in FH_vcf:
                readed_records.append(
                    "_".join([variant.chrom, str(variant.pos), variant.ref, ",".join(variant.alt)])
                )
            self.assertEqual(expected_records, readed_records)

    def testGetSub(self):
        try:
            pysam.tabix_compress(self.tmp_in_variants, self.tmp_in_variants + ".gz")
            pysam.tabix_index(self.tmp_in_variants + ".gz", preset="vcf")
            with VCFIO(self.tmp_in_variants + ".gz", "i") as reader:
                # Missing chromosome
                observed = [elt.getName() for elt in reader.getSub("1", 1437, 11437)]
                self.assertEqual(observed, [])
                # No variants
                observed = [elt.getName() for elt in reader.getSub("20", 1234580, 1234590)]
                self.assertEqual(observed, [])
                # At substit
                observed = [elt.getName() for elt in reader.getSub("20", 17330, 17330)]
                self.assertEqual(observed, ["20:17330=T/A"])
                # Contains several variants
                observed = [elt.getName() for elt in reader.getSub("20", 14330, 17350)]
                self.assertEqual(observed, ["20:14370=G/A", "20:17330=T/A"])
                # Overlap deletion start
                observed = [elt.getName() for elt in reader.getSub("20", 1234550, 1234567)]
                self.assertEqual(observed, ["20:1234567=GTC/G/GTCT"])
                # In deletion
                observed = [elt.getName() for elt in reader.getSub("20", 1234569, 1234569)]
                self.assertEqual(observed, ["20:1234567=GTC/G/GTCT"])
        finally:
            for curr_file in [self.tmp_in_variants + ".gz", self.tmp_in_variants + ".gz.tbi"]:
                if os.path.exists(curr_file):
                    os.remove(curr_file)

    def testParseHeader(self):
        with VCFIO(self.tmp_in_variants_with_spl_info) as FH_vcf:
            # FORMAT
            self.assertEqual(
                {key: elt.datastore for key, elt in FH_vcf.format.items()},
                {
                    'GT': {'id': 'GT', 'description': 'Genotype', 'type': 'String', 'number': '1'},
                    'GQ': {'id': 'GQ', 'description': 'Genotype Quality', 'type': 'Integer', 'number': '1'},
                    'DP': {'id': 'DP', 'description': 'Read Depth', 'type': 'Integer', 'number': '1'},
                    'HQ': {'id': 'HQ', 'description': 'Haplotype Quality', 'type': 'Integer', 'number': '2'}
                }
            )
            # INFO
            self.assertEqual(
                {key: elt.datastore for key, elt in FH_vcf.info.items()},
                {
                    'NS': {'id': 'NS', 'description': 'Number of Samples With Data', 'type': 'Integer', 'number': '1'},
                    'DP': {'id': 'DP', 'description': 'Total Depth', 'type': 'Integer', 'number': '1'},
                    'AF': {'id': 'AF', 'description': 'Allele Frequency', 'type': 'Float', 'number': 'A'},
                    'AA': {'id': 'AA', 'description': 'Ancestral Allele', 'type': 'String', 'number': '1'},
                    'DB': {'id': 'DB', 'description': 'dbSNP membership, build 129', 'type': 'Flag', 'number': '0'},
                    'H2': {'id': 'H2', 'description': 'HapMap2 membership', 'type': 'Flag', 'number': '0'}}
            )
            # Samples list
            self.assertEqual(FH_vcf.samples, ["NA00001", "NA00002", "NA00003"])
            # SAMPLE
            self.assertDictEqual(
                {key: elt.datastore for key, elt in FH_vcf.sample_info.items()},
                {
                    "NA00001": {
                        "id": "NA00001",
                        "assay": "WholeGenome",
                        "ethnicity": "AFR",
                        "disease": "None",
                        "description": "Patient germline genome from unaffected"
                    },
                    "NA00003": {
                        "id": "NA00003",
                        "assay": "Exome",
                        "ethnicity": "CEU",
                        "disease": "Cancer",
                        "description": "European patient exome from breast cancer",
                        "tissue": "Breast"
                    }
                }
            )
            # EXTRA
            self.assertEqual(
                FH_vcf.extra_header,
                [
                    '##fileDate=20090805"',
                    '##source=myImputationProgramV3.1',
                    '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta',
                    '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
                    '##phasing=partial'
                ]
            )

    def testReadTypeFlag(self):
        expected = [
            {"DB": True, "H2": True},
            {"DB": False, "H2": False},
            {"DB": True, "H2": False},
            {"DB": False, "H2": False},
            {"DB": False, "H2": False}
        ]
        observed = []
        with VCFIO(self.tmp_in_variants) as FH_vcf:
            for record in FH_vcf:
                observed.append({
                    "DB": "DB" in record.info,
                    "H2": "H2" in record.info
                })
        self.assertEqual(expected, observed)

    def testWrite(self):
        expected_content = """##fileformat=VCFv4.3
##fileDate=20090805"
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	A	29.0	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3.0	q10	AF=0.017;DP=11;NS=3	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.
20	1110696	rs6040355	A	G,T	67.0	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.
20	1230237	.	T	-	47.0	PASS	AA=T;DP=13;NS=3	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
20	1234567	microsat1	GTC	G,GTCT	50.0	PASS	AA=G;DP=9;NS=3	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3"""
        # Read and write VCF
        with VCFIO(self.tmp_in_variants) as reader:
            with VCFIO(self.tmp_out_variants, "w") as writer:
                writer.copyHeader(reader)
                writer.writeHeader()
                for record in reader:
                    writer.write(record)
        # Compare input content with output
        observed_content = None
        with open(self.tmp_out_variants) as reader:
            observed_content = "".join(reader.readlines()).strip()
        self.assertEqual(expected_content, observed_content)

    def testWriteWithSplInfo(self):
        expected_content = """##fileformat=VCFv4.3
##fileDate=20090805"
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##SAMPLE=<ID=NA00001,Assay="WholeGenome",Description="Patient germline genome from unaffected",Disease="None",Ethnicity="AFR">
##SAMPLE=<ID=NA00003,Assay="Exome",Description="European patient exome from breast cancer",Disease="Cancer",Ethnicity="CEU",Tissue="Breast">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	A	29.0	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3.0	q10	AF=0.017;DP=11;NS=3	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.
20	1110696	rs6040355	A	G,T	67.0	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.
20	1230237	.	T	-	47.0	PASS	AA=T;DP=13;NS=3	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
20	1234567	microsat1	GTC	G,GTCT	50.0	PASS	AA=G;DP=9;NS=3	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3"""
        # Read and write VCF
        with VCFIO(self.tmp_in_variants_with_spl_info) as reader:
            with VCFIO(self.tmp_out_variants, "w") as writer:
                writer.copyHeader(reader)
                writer.writeHeader()
                for record in reader:
                    writer.write(record)
        # Compare input content with output
        observed_content = None
        with open(self.tmp_out_variants) as reader:
            observed_content = "".join(reader.readlines()).strip()
        self.assertEqual(expected_content, observed_content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_variants, self.tmp_in_variants_with_spl_info, self.tmp_out_variants]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


class TestVCFRecord(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fa = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + ".fa.fai")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")

        # AD, AF and DP evaluation
        self.freq_data = dict()
        self.freq_data["wout_spl_rAD_aAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
11	1	wout_spl_rAD_aAF_1	A	T	29	PASS	expModel=model_0;AD=90,10;DP=100	.
11	2	wout_spl_rAD_aAF_2	A	T	29	PASS	expModel=model_0;AF=0.1;AD=90,10	.
11	3	wout_spl_rAD_aAF_3	A	T	29	PASS	expModel=model_0;AF=0.1;DP=100	.
11	1	wout_spl_rAD_aAF_4	A	T,G	29	PASS	expModel=model_1;AD=85,10,5;DP=100	.
11	2	wout_spl_rAD_aAF_5	A	T,G	29	PASS	expModel=model_1;AF=0.1,0.05;AD=85,10,5	.
11	3	wout_spl_rAD_aAF_6	A	T,G	29	PASS	expModel=model_1;AF=0.1,0.05;DP=100	."""
        self.freq_data["wout_spl_aAD_rAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
11	1	wout_spl_aAD_rAF_1	A	T	29	PASS	expModel=model_0;AD=10;DP=100	.
11	2	wout_spl_aAD_rAF_2	A	T	29	PASS	expModel=model_0;AF=0.9,0.1;AD=10	.
11	3	wout_spl_aAD_rAF_3	A	T	29	PASS	expModel=model_0;AF=0.9,0.1;DP=100	.
11	1	wout_spl_aAD_rAF_4	A	T,G	29	PASS	expModel=model_1;AD=10,5;DP=100	.
11	2	wout_spl_aAD_rAF_5	A	T,G	29	PASS	expModel=model_1;AF=0.85,0.1,0.05;AD=10,5	.
11	3	wout_spl_aAD_rAF_6	A	T,G	29	PASS	expModel=model_1;AF=0.85,0.1,0.05;DP=100	."""
        self.freq_data["one_spl_rAD_aAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
##INFO=<ID=skipAD,Number=1,Type=String,Description="Variant must be skipped by AD calculation">
##INFO=<ID=skipAF,Number=1,Type=String,Description="Variant must be skipped by AF calculation">
##INFO=<ID=skipDP,Number=1,Type=String,Description="Variant must be skipped by DP calculation">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
11	1	one_spl_rAD_aAF_1	A	T	29	PASS	expModel=model_2;skipAF=True;AD=90,10;DP=100	.	.
11	2	one_spl_rAD_aAF_2	A	T	29	PASS	expModel=model_2;skipDP=True;AF=0.1;AD=90,10	.	.
11	3	one_spl_rAD_aAF_3	A	T	29	PASS	expModel=model_2;skipAD=True;AF=0.1;DP=100	.	.
11	1	one_spl_rAD_aAF_4	A	T	29	PASS	expModel=model_2	AD:DP	90,10:100
11	2	one_spl_rAD_aAF_5	A	T	29	PASS	expModel=model_2	AF:AD	0.1:90,10
11	3	one_spl_rAD_aAF_6	A	T	29	PASS	expModel=model_2	AF:DP	0.1:100
11	1	one_spl_rAD_aAF_7	A	T,G	29	PASS	expModel=model_3;skipAF=True;AD=85,10,5;DP=100	.	.
11	2	one_spl_rAD_aAF_8	A	T,G	29	PASS	expModel=model_3;skipDP=True;AF=0.1,0.05;AD=85,10,5	.	.
11	3	one_spl_rAD_aAF_9	A	T,G	29	PASS	expModel=model_3;skipAD=True;AF=0.1,0.05;DP=100	.	.
11	1	one_spl_rAD_aAF_10	A	T,G	29	PASS	expModel=model_3	AD:DP	85,10,5:100
11	2	one_spl_rAD_aAF_11	A	T,G	29	PASS	expModel=model_3	AF:AD	0.1,0.05:85,10,5
11	3	one_spl_rAD_aAF_12	A	T,G	29	PASS	expModel=model_3	AF:DP	0.1,0.05:100"""
        self.freq_data["one_spl_aAD_rAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
##INFO=<ID=skipAD,Number=1,Type=String,Description="Variant must be skipped by AD calculation">
##INFO=<ID=skipAF,Number=1,Type=String,Description="Variant must be skipped by AF calculation">
##INFO=<ID=skipDP,Number=1,Type=String,Description="Variant must be skipped by DP calculation">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
11	1	one_spl_aAD_rAF_1	A	T	29	PASS	expModel=model_2;skipAF=True;AD=10;DP=100	.	.
11	2	one_spl_aAD_rAF_2	A	T	29	PASS	expModel=model_2;skipDP=True;AF=0.9,0.1;AD=10	.	.
11	3	one_spl_aAD_rAF_3	A	T	29	PASS	expModel=model_2;skipAD=True;AF=0.9,0.1;DP=100	.	.
11	1	one_spl_aAD_rAF_4	A	T	29	PASS	expModel=model_2	AD:DP	10:100
11	2	one_spl_aAD_rAF_5	A	T	29	PASS	expModel=model_2	AF:AD	0.9,0.1:10
11	3	one_spl_aAD_rAF_6	A	T	29	PASS	expModel=model_2	AF:DP	0.9,0.1:100
11	1	one_spl_rAD_aAF_7	A	T,G	29	PASS	expModel=model_3;skipAF=True;AD=10,5;DP=100	.	.
11	2	one_spl_rAD_aAF_8	A	T,G	29	PASS	expModel=model_3;skipDP=True;AF=0.85,0.1,0.05;AD=10,5	.	.
11	3	one_spl_rAD_aAF_9	A	T,G	29	PASS	expModel=model_3;skipAD=True;AF=0.85,0.1,0.05;DP=100	.	.
11	1	one_spl_rAD_aAF_10	A	T,G	29	PASS	expModel=model_3	AD:DP	10,5:100
11	2	one_spl_rAD_aAF_11	A	T,G	29	PASS	expModel=model_3	AF:AD	0.85,0.1,0.05:10,5
11	3	one_spl_rAD_aAF_12	A	T,G	29	PASS	expModel=model_3	AF:DP	0.85,0.1,0.05:100"""
        self.freq_data["two_spl_rAD_aAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
##INFO=<ID=skipAD,Number=1,Type=String,Description="Variant must be skipped by AD calculation">
##INFO=<ID=skipAF,Number=1,Type=String,Description="Variant must be skipped by AF calculation">
##INFO=<ID=skipDP,Number=1,Type=String,Description="Variant must be skipped by DP calculation">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
11	1	two_spl_rAD_aAF_1	A	T	29	PASS	expModel=model_4	AD:DP	48,2:50	42,8:50
11	2	two_spl_rAD_aAF_2	A	T	29	PASS	expModel=model_4	AF:AD	0.04:48,2	0.16:42,8
11	3	two_spl_rAD_aAF_3	A	T	29	PASS	expModel=model_4	AF:DP	0.04:50	0.16:50
11	1	two_spl_rAD_aAF_4	A	T,G	29	PASS	expModel=model_5	AD:DP	43,2,5:50	41,8,1:50
11	2	two_spl_rAD_aAF_5	A	T,G	29	PASS	expModel=model_5	AF:AD	0.04,0.1:43,2,5	0.16,0.02:41,8,1
11	3	two_spl_rAD_aAF_6	A	T,G	29	PASS	expModel=model_5	AF:DP	0.04,0.1:50	0.16,0.02:50"""
        self.freq_data["two_spl_aAD_rAF"] = """##fileformat=VCFv4.0
##INFO=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=expModel,Number=1,Type=String,Description="Expected model for variant">
##INFO=<ID=skipAD,Number=1,Type=String,Description="Variant must be skipped by AD calculation">
##INFO=<ID=skipAF,Number=1,Type=String,Description="Variant must be skipped by AF calculation">
##INFO=<ID=skipDP,Number=1,Type=String,Description="Variant must be skipped by DP calculation">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=AF,Number=R,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
11	1	two_spl_rAD_aAF_1	A	T	29	PASS	expModel=model_4	AD:DP	2:50	8:50
11	2	two_spl_rAD_aAF_2	A	T	29	PASS	expModel=model_4	AF:AD	0.96,0.04:2	0.84,0.16:8
11	3	two_spl_rAD_aAF_3	A	T	29	PASS	expModel=model_4	AF:DP	0.96,0.04:50	0.84,0.16:50
11	1	two_spl_rAD_aAF_4	A	T,G	29	PASS	expModel=model_5	AD:DP	2,5:50	8,1:50
11	2	two_spl_rAD_aAF_5	A	T,G	29	PASS	expModel=model_5	AF:AD	0.96,0.04,0.1:2,5	0.84,0.16,0.02:8,1
11	3	two_spl_rAD_aAF_6	A	T,G	29	PASS	expModel=model_5	AF:DP	0.96,0.04,0.1:50	0.84,0.16,0.02:50"""
        self.freq_expected = {
            "model_0": {
                "AF": {"pop": [0.1]},
                "AD": {"pop": [10]},
                "DP": {"pop": 100}
            },
            "model_1": {
                "AF": {"pop": [0.1, 0.05]},
                "AD": {"pop": [10, 5]},
                "DP": {"pop": 100}
            },
            "model_2": {
                "AF": {"splA": [0.1], "pop": [0.1]},
                "AD": {"splA": [10], "pop": [10]},
                "DP": {"splA": 100, "pop": 100}
            },
            "model_3": {
                "AF": {"splA": [0.1, 0.05], "pop": [0.1, 0.05]},
                "AD": {"splA": [10, 5], "pop": [10, 5]},
                "DP": {"splA": 100, "pop": 100}
            },
            "model_4": {
                "AF": {"splA": [0.04], "splB": [0.16], "pop": [0.1]},
                "AD": {"splA": [2], "splB": [8], "pop": [10]},
                "DP": {"splA": 50, "splB": 50, "pop": 100}
            },
            "model_5": {
                "AF": {"splA": [0.04, 0.1], "splB": [0.16, 0.02], "pop": [0.1, 0.06]},
                "AD": {"splA": [2, 5], "splB": [8, 1], "pop": [10, 6]},
                "DP": {"splA": 50, "splB": 50, "pop": 100}
            }
        }

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fa, self.tmp_fai, self.tmp_variants]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGetAltAF(self):
        for curr_dataset in sorted(self.freq_data):
            if not curr_dataset.startswith("wout_spl"):
                # Write test data
                content = self.freq_data[curr_dataset]
                with open(self.tmp_variants, "w") as FH_variants:
                    FH_variants.write(content)
                # Parse
                expected_records = list()
                observed_records = list()
                with VCFIO(self.tmp_variants) as FH_vcf:
                    for variant in FH_vcf:
                        if "skipAF" not in variant.info:
                            try:
                                for idx_spl, curr_spl in enumerate(FH_vcf.samples):
                                    expected_AF = self.freq_expected[variant.info["expModel"]]["AF"][curr_spl]
                                    expected_records.append("{} {}".format(variant.id, expected_AF))
                                    observed_records.append("{} {}".format(variant.id, variant.getAltAF(curr_spl)))
                            except Exception:
                                raise Exception('Error in TestVCFRecord.testGetAltAF() for variant "' + variant.id + '".')
                # Assert
                self.assertEqual(expected_records, observed_records)

    def testGetAltAD(self):
        for curr_dataset in sorted(self.freq_data):
            if not curr_dataset.startswith("wout_spl"):
                # Write test data
                content = self.freq_data[curr_dataset]
                with open(self.tmp_variants, "w") as FH_variants:
                    FH_variants.write(content)
                # Parse
                expected_records = list()
                observed_records = list()
                with VCFIO(self.tmp_variants) as FH_vcf:
                    for variant in FH_vcf:
                        if "skipAD" not in variant.info:
                            try:
                                for idx_spl, curr_spl in enumerate(FH_vcf.samples):
                                    expected_AD = self.freq_expected[variant.info["expModel"]]["AD"][curr_spl]
                                    expected_records.append("{} {}".format(variant.id, expected_AD))
                                    observed_records.append("{} {}".format(variant.id, variant.getAltAD(curr_spl)))
                            except Exception:
                                raise Exception('Error in TestVCFRecord.testGetAltAD() for variant "' + variant.id + '".')
                # Assert
                self.assertEqual(expected_records, observed_records)

    def testGetDP(self):
        for curr_dataset in sorted(self.freq_data):
            if not curr_dataset.startswith("wout_spl"):
                # Write test data
                content = self.freq_data[curr_dataset]
                with open(self.tmp_variants, "w") as FH_variants:
                    FH_variants.write(content)
                # Parse
                expected_records = list()
                observed_records = list()
                with VCFIO(self.tmp_variants) as FH_vcf:
                    for variant in FH_vcf:
                        if "skipDP" not in variant.info:
                            try:
                                for curr_spl in FH_vcf.samples:
                                    expected_DP = self.freq_expected[variant.info["expModel"]]["DP"][curr_spl]
                                    expected_records.append("{} {}".format(variant.id, expected_DP))
                                    observed_records.append("{} {}".format(variant.id, variant.getDP(curr_spl)))
                            except Exception:
                                raise Exception('Error in TestVCFRecord.testGetDP() for variant "' + variant.id + '".')
                # Assert
                self.assertEqual(expected_records, observed_records)

    def testGetPopAltAD(self):
        for curr_dataset in sorted(self.freq_data):
            # Write test data
            content = self.freq_data[curr_dataset]
            with open(self.tmp_variants, "w") as FH_variants:
                FH_variants.write(content)
            # Parse
            expected_records = list()
            observed_records = list()
            with VCFIO(self.tmp_variants) as FH_vcf:
                for variant in FH_vcf:
                    try:
                        expected_AD = self.freq_expected[variant.info["expModel"]]["AD"]["pop"]
                        expected_records.append("{} {}".format(variant.id, expected_AD))
                        observed_records.append("{} {}".format(variant.id, variant.getPopAltAD()))
                    except Exception:
                        raise Exception('Error in TestVCFRecord.testGetPopAltAD() for variant "' + variant.id + '".')
            # Assert
            self.assertEqual(expected_records, observed_records)

    def testGetPopAltAF(self):
        for curr_dataset in sorted(self.freq_data):
            # Write test data
            content = self.freq_data[curr_dataset]
            with open(self.tmp_variants, "w") as FH_variants:
                FH_variants.write(content)
            # Parse
            expected_records = list()
            observed_records = list()
            with VCFIO(self.tmp_variants) as FH_vcf:
                for variant in FH_vcf:
                    try:
                        expected_AF = self.freq_expected[variant.info["expModel"]]["AF"]["pop"]
                        expected_records.append("{} {}".format(variant.id, expected_AF))
                        observed_records.append("{} {}".format(variant.id, variant.getPopAltAF()))
                    except Exception:
                        raise Exception('Error in TestVCFRecord.testGetPopAltAF() for variant "' + variant.id + '".')
            # Assert
            self.assertEqual(expected_records, observed_records)

    def testGetPopDP(self):
        for curr_dataset in sorted(self.freq_data):
            # Write test data
            content = self.freq_data[curr_dataset]
            with open(self.tmp_variants, "w") as FH_variants:
                FH_variants.write(content)
            # Parse
            expected_records = list()
            observed_records = list()
            with VCFIO(self.tmp_variants) as FH_vcf:
                for variant in FH_vcf:
                    try:
                        expected_DP = self.freq_expected[variant.info["expModel"]]["DP"]["pop"]
                        expected_records.append("{} {}".format(variant.id, expected_DP))
                        observed_records.append("{} {}".format(variant.id, variant.getPopDP()))
                    except Exception:
                        raise Exception('Error in TestVCFRecord.testGetPopDP() for variant "' + variant.id + '".')
            # Assert
            self.assertEqual(expected_records, observed_records)

    def testNormalizeSingleAllele(self):
        # Test substitution one nt
        substitution = VCFRecord("artificial_1", 18, None, "TC", ["TA"], 230)
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "C" and substitution.alt[0] == "A" and substitution.pos == 19)
        # Test substitution multi nt
        substitution = VCFRecord("artificial_1", 18, None, "TCtgA", ["TaGCc"], 230)
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "CTGA" and substitution.alt[0] == "AGCC" and substitution.pos == 19)
        # Test substitution multi nt with possible remove
        substitution = VCFRecord("artificial_1", 18, None, "TCtgATT", ["TaGGctT"], 230)
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "CTGA" and substitution.alt[0] == "AGGC" and substitution.pos == 19)

        # Test insertion
        insertion = VCFRecord("artificial_1", 18, None, "T", ["TA"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "-" and insertion.alt[0] == "A" and insertion.pos == 19)
        # Test insertion multi nt and remove upstream and downstream
        insertion = VCFRecord("artificial_1", 18, None, "TGAT", ["TCGAGAT"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "-" and insertion.alt[0] == "CGA" and insertion.pos == 19)
        # Test insertion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord("artificial_1", 18, None, "TCtgATTAGC", ["TaGGctTATGCGC"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "CTGATTA" and insertion.alt[0] == "AGGCTTATGC" and insertion.pos == 19)

        # Test deletion
        insertion = VCFRecord("artificial_1", 18, None, "TA", ["T"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "A" and insertion.alt[0] == "-" and insertion.pos == 19)
        # Test insertion multi nt and remove upstream and downstream
        insertion = VCFRecord("artificial_1", 18, None, "TCGAGAT", ["TGAT"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "CGA" and insertion.alt[0] == "-" and insertion.pos == 19)
        # Test insertion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord("artificial_1", 18, None, "TaGGctTATGCGC", ["TCtgATTAGC"], 230)
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "AGGCTTATGC" and insertion.alt[0] == "CTGATTA" and insertion.pos == 19)

    def testGetMostUpstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  32
        # Test fix deletion
        deletion = VCFRecord("artificial_1", 18, None, "TtTaAGC", ["T"], 230)
        upstream = deletion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 19 and upstream.ref == "TTAAGC" and upstream.alt[0] == "-")
        # Test homopolymer deletion
        deletion = VCFRecord("artificial_1", 18, None, "TTT", ["T"], 230)
        upstream = deletion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 17 and upstream.ref == "TT" and upstream.alt[0] == "-")
        # Test complex deletion
        deletion = VCFRecord("artificial_1", 25, None, "CGAgCC", ["C"], 230)
        upstream = deletion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 22 and upstream.ref == "AGCCG" and upstream.alt[0] == "-")
        # Test deletion on start
        deletion = VCFRecord("artificial_1", 1, None, ref[0], ["-"], 230)
        upstream = deletion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 1 and upstream.ref == "N" and upstream.alt[0] == "-")
        # Test deletion on end
        deletion = VCFRecord("artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
        upstream = deletion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "-")
        # Test fix insertion
        insertion = VCFRecord("artificial_1", 4, None, "N", ["NGGTT"], 230)
        upstream = insertion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 5 and upstream.ref == "-" and upstream.alt[0] == "GGTT")
        # Test homopolymer insertion
        insertion = VCFRecord("artificial_1", 18, None, "T", ["TT"], 230)
        upstream = insertion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 17 and upstream.ref == "-" and upstream.alt[0] == "T")
        # Test complex insertion
        insertion = VCFRecord("artificial_1", 25, None, "C", ["CGAgCC"], 230)
        upstream = insertion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 22 and upstream.ref == "-" and upstream.alt[0] == "AGCCG")
        # Test insertion on start
        insertion = VCFRecord("artificial_1", 1, None, "-", ["A"], 230)
        upstream = insertion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 1 and upstream.ref == "-" and upstream.alt[0] == "A")
        # Test insertion on end
        insertion = VCFRecord("artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230)
        upstream = insertion.getMostUpstream(ref)
        self.assertTrue(upstream.pos == 34 and upstream.ref == "-" and upstream.alt[0] == "G")

    def testGetMostDownstream(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  32
        # Test fix deletion
        deletion = VCFRecord("artificial_1", 18, None, "TtTaAGC", ["T"], 230)
        downstream = deletion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 19 and downstream.ref == "TTAAGC" and downstream.alt[0] == "-")
        # Test homopolymer deletion
        deletion = VCFRecord("artificial_1", 18, None, "TTT", ["T"], 230)
        downstream = deletion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 19 and downstream.ref == "TT" and downstream.alt[0] == "-")
        # Test complex deletion
        deletion = VCFRecord("artificial_1", 25, None, "CGAgCC", ["C"], 230)
        downstream = deletion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 28 and downstream.ref == "GCCGA" and downstream.alt[0] == "-")
        # Test deletion on start
        deletion = VCFRecord("artificial_1", 1, None, ref[0], ["-"], 230)
        downstream = deletion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "-")
        # Test deletion on end
        deletion = VCFRecord("artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
        downstream = deletion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "-")
        # Test fix insertion
        insertion = VCFRecord("artificial_1", 4, None, "N", ["NGGTT"], 230)
        downstream = insertion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 5 and downstream.ref == "-" and downstream.alt[0] == "GGTT")
        # Test homopolymer insertion
        insertion = VCFRecord("artificial_1", 18, None, "T", ["TT"], 230)
        downstream = insertion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 21 and downstream.ref == "-" and downstream.alt[0] == "T")
        # Test complex insertion
        insertion = VCFRecord("artificial_1", 25, None, "C", ["CGAgCC"], 230)
        downstream = insertion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 33 and downstream.ref == "-" and downstream.alt[0] == "GCCGA")
        # Test insertion on start
        insertion = VCFRecord("artificial_1", 1, None, "-", ["A"], 230)
        downstream = insertion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 1 and downstream.ref == "-" and downstream.alt[0] == "A")
        # Test insertion on end
        insertion = VCFRecord("artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230)
        downstream = insertion.getMostDownstream(ref)
        self.assertTrue(downstream.pos == 34 and downstream.ref == "-" and downstream.alt[0] == "G")

    def testFastDownstreamed(self):
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        with open(self.tmp_fa, "w") as writer:
            writer.write(">artificial_1\n")
            writer.write("nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT")
        with open(self.tmp_fai, "w") as writer:
            writer.write("artificial_1\t33\t14\t33\t34")
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test fix deletion
            deletion = VCFRecord("artificial_1", 18, None, "TTTaAGC", ["T"], 230)
            downstream = deletion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTTAAGC" and downstream.alt[0] == "T")
            # Test homopolymer deletion
            deletion = VCFRecord("artificial_1", 16, None, "GTT", ["G"], 230)
            downstream = deletion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTT" and downstream.alt[0] == "T")
            # Test complex deletion
            deletion = VCFRecord("artificial_1", 25, None, "CGAgCC", ["C"], 230)
            downstream = deletion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 27 and downstream.ref == "AGCCGA" and downstream.alt[0] == "A")
            # Test deletion on start
            deletion = VCFRecord("artificial_1", 1, None, ref[0], ["-"], 230)
            downstream = deletion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 3 and downstream.ref == "NN" and downstream.alt[0] == "N")
            # Test deletion on end
            deletion = VCFRecord("artificial_1", len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
            downstream = deletion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 32 and downstream.ref == "AT" and downstream.alt[0] == "A")
            # Test fix insertion
            insertion = VCFRecord("artificial_1", 4, None, "N", ["NGGTT"], 230)
            downstream = insertion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "NGGTT")
            # Test homopolymer insertion
            insertion = VCFRecord("artificial_1", 18, None, "T", ["TT"], 230)
            downstream = insertion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 20 and downstream.ref == "T" and downstream.alt[0] == "TT")
            # Test complex insertion
            insertion = VCFRecord("artificial_1", 25, None, "C", ["CGAgCC"], 230)
            downstream = insertion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 32 and downstream.ref == "A" and downstream.alt[0] == "AGCCGA")
            # Test insertion on end
            insertion = VCFRecord("artificial_1", len(ref), None, ref[-1], [ref[-1] + "G"], 230)
            downstream = insertion.fastDownstreamed(seq_handler, 10)
            self.assertTrue(downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "TG")

    def testRefStart(self):
        data = [
            # Substitutions
            ["artificial_1", 10, "substit_1", "A", ["T"], None, None, {"expected_start": 10, "expected_end": 10}],
            ["artificial_1", 10, "substit_2", "AA", ["TT"], None, None, {"expected_start": 10, "expected_end": 11}],
            # Deletions
            ["artificial_1", 12, "del_1", "A", ["-"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "del_2", "AAA", ["-"], None, None, {"expected_start": 10, "expected_end": 12}],
            ["artificial_1", 10, "del_3", "AAA", ["A"], None, None, {"expected_start": 11, "expected_end": 12}],
            ["artificial_1", 10, "del_4", "AGC", ["A"], None, None, {"expected_start": 11, "expected_end": 12}],
            ["artificial_1", 10, "del_5", "AAAT", ["TG"], None, None, {"expected_start": 10, "expected_end": 13}],
            # Insertions
            ["artificial_1", 12, "ins_1", "A", ["TG"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 12, "ins_2", "A", ["AGT"], None, None, {"expected_start": 12.5, "expected_end": 12.5}],
            ["artificial_1", 11, "ins_3", "AA", ["AGT"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "ins_4", "AAA", ["AAGT"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "ins_5", "AATC", ["AACGT"], None, None, {"expected_start": 12, "expected_end": 13}],
        ]
        for curr_data in data:
            record = VCFRecord(*curr_data)
            self.assertTrue(record.refStart() == record.info["expected_start"])

    def testRefEnd(self):
        data = [
            # Substitutions
            ["artificial_1", 10, "substit_1", "A", ["T"], None, None, {"expected_start": 10, "expected_end": 10}],
            ["artificial_1", 10, "substit_2", "AA", ["TT"], None, None, {"expected_start": 10, "expected_end": 11}],
            # Deletions
            ["artificial_1", 12, "del_1", "A", ["-"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "del_2", "AAA", ["-"], None, None, {"expected_start": 10, "expected_end": 12}],
            ["artificial_1", 10, "del_3", "AAA", ["A"], None, None, {"expected_start": 11, "expected_end": 12}],
            ["artificial_1", 10, "del_4", "AGC", ["A"], None, None, {"expected_start": 11, "expected_end": 12}],
            ["artificial_1", 10, "del_5", "AAAT", ["TG"], None, None, {"expected_start": 10, "expected_end": 13}],
            # Insertions
            ["artificial_1", 12, "ins_1", "A", ["TG"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 12, "ins_2", "A", ["AGT"], None, None, {"expected_start": 12.5, "expected_end": 12.5}],
            ["artificial_1", 11, "ins_3", "AA", ["AGT"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "ins_4", "AAA", ["AAGT"], None, None, {"expected_start": 12, "expected_end": 12}],
            ["artificial_1", 10, "ins_5", "AATC", ["AACGT"], None, None, {"expected_start": 12, "expected_end": 13}],
        ]
        for curr_data in data:
            record = VCFRecord(*curr_data)
            self.assertTrue(record.refEnd() == record.info["expected_end"])

    def testIsDeletion(self):
        data = [
            {"variant": VCFRecord("chr1", 2, None, "A", ["T"]), "expected": False},
            {"variant": VCFRecord("chr1", 2, None, "A", ["AT"]), "expected": False},
            {"variant": VCFRecord("chr1", 1, None, "CA", ["C"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GT"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GTT"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CAC", ["GT"]), "expected": True}
        ]
        for curr in data:
            self.assertEqual(curr["variant"].isDeletion(), curr["expected"])

    def testIsInsAndDel(self):
        data = [
            {"variant": VCFRecord("chr1", 2, None, "A", ["T"]), "expected": False},
            {"variant": VCFRecord("chr1", 2, None, "A", ["AT"]), "expected": False},
            {"variant": VCFRecord("chr1", 1, None, "CA", ["C"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GT"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GTT"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CAC", ["GT"]), "expected": True}
        ]
        for curr in data:
            self.assertEqual(curr["variant"].isInsAndDel(), curr["expected"])

    def testIsInsertion(self):
        data = [
            {"variant": VCFRecord("chr1", 2, None, "A", ["T"]), "expected": False},
            {"variant": VCFRecord("chr1", 2, None, "A", ["AT"]), "expected": True},
            {"variant": VCFRecord("chr1", 1, None, "CA", ["C"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GT"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GTT"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CAC", ["GT"]), "expected": False}
        ]
        for curr in data:
            self.assertEqual(curr["variant"].isInsertion(), curr["expected"])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
