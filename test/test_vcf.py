#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.16.0'
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

from anacore.vcf import _refDownIter, _refUpIter, getHeaderAttr, HeaderInfoAttr, VCFIO, VCFRecord, VCFSymbAltRecord


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
11	2	two_spl_rAD_aAF_5	A	T,G	29	PASS	expModel=model_5	AF:AD	0.86,0.04,0.1:2,5	0.82,0.16,0.02:8,1
11	3	two_spl_rAD_aAF_6	A	T,G	29	PASS	expModel=model_5	AF:DP	0.86,0.04,0.1:50	0.82,0.16,0.02:50"""
        self.freq_expected = {
            "model_0": {
                "AF": {"pop": [0.1], "pop_ref": 0.9},
                "AD": {"pop": [10], "pop_ref": 90},
                "DP": {"pop": 100}
            },
            "model_1": {
                "AF": {"pop": [0.1, 0.05], "pop_ref": 0.85},
                "AD": {"pop": [10, 5], "pop_ref": 85},
                "DP": {"pop": 100}
            },
            "model_2": {
                "AF": {"splA": [0.1], "pop": [0.1], "pop_ref": 0.9},
                "AD": {"splA": [10], "pop": [10], "pop_ref": 90},
                "DP": {"splA": 100, "pop": 100}
            },
            "model_3": {
                "AF": {"splA": [0.1, 0.05], "pop": [0.1, 0.05], "pop_ref": 0.85},
                "AD": {"splA": [10, 5], "pop": [10, 5], "pop_ref": 85},
                "DP": {"splA": 100, "pop": 100}
            },
            "model_4": {
                "AF": {"splA": [0.04], "splB": [0.16], "pop": [0.1], "pop_ref": 0.9},
                "AD": {"splA": [2], "splB": [8], "pop": [10], "pop_ref": 90},
                "DP": {"splA": 50, "splB": 50, "pop": 100}
            },
            "model_5": {
                "AF": {"splA": [0.04, 0.1], "splB": [0.16, 0.02], "pop": [0.1, 0.06], "pop_ref": 0.84},
                "AD": {"splA": [2, 5], "splB": [8, 1], "pop": [10, 6], "pop_ref": 84},
                "DP": {"splA": 50, "splB": 50, "pop": 100}
            }
        }

    def testContainsIndel(self):
        data = [
            {"variant": VCFRecord("chr1", 2, None, "A", ["T"]), "expected": False},
            {"variant": VCFRecord("chr1", 2, None, "A", ["AT"]), "expected": True},
            {"variant": VCFRecord("chr1", 1, None, "CA", ["C"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GT"]), "expected": False},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GTT"]), "expected": True},
            {"variant": VCFRecord("chr1", 15, None, "CAC", ["GT"]), "expected": True}
        ]
        for curr in data:
            self.assertEqual(curr["variant"].containsIndel(), curr["expected"])

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

    def testGetPopRefAD(self):
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
                        expected_AD = self.freq_expected[variant.info["expModel"]]["AD"]["pop_ref"]
                        expected_records.append("{} {}".format(variant.id, expected_AD))
                        observed_records.append("{} {}".format(variant.id, variant.getPopRefAD()))
                    except Exception:
                        raise Exception('Error in TestVCFRecord.testGetPopRefAD() for variant "{}".'.format(variant.id))
            # Assert
            self.assertEqual(expected_records, observed_records)

    def testGetPopRefAF(self):
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
                        expected_AF = self.freq_expected[variant.info["expModel"]]["AF"]["pop_ref"]
                        expected_records.append("{} {}".format(variant.id, expected_AF))
                        observed_records.append("{} {}".format(variant.id, variant.getPopRefAF()))
                    except Exception:
                        raise Exception('Error in TestVCFRecord.testGetPopRefAF() for variant "{}".'.format(variant.id))
            # Assert
            self.assertEqual(expected_records, observed_records)

    def testNormalizeSingleAllele(self):
        chrom = "artificial_1"
        # Test substitution one nt already normalized
        substitution = VCFRecord(chrom, 19, None, "C", ["A"])
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "C" and substitution.alt[0] == "A" and substitution.pos == 19)
        # Test substitution one nt
        substitution = VCFRecord(chrom, 18, None, "TC", ["TA"])
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "C" and substitution.alt[0] == "A" and substitution.pos == 19)
        # Test substitution multi nt
        substitution = VCFRecord(chrom, 18, None, "TCtgA", ["TaGCc"])
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "CTGA" and substitution.alt[0] == "AGCC" and substitution.pos == 19)
        # Test substitution multi nt with possible remove
        substitution = VCFRecord(chrom, 18, None, "TCtgATT", ["TaGGctT"])
        substitution.normalizeSingleAllele()
        self.assertTrue(substitution.ref == "CTGA" and substitution.alt[0] == "AGGC" and substitution.pos == 19)
        # Test insertion
        insertion = VCFRecord(chrom, 18, None, "T", ["TA"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "-" and insertion.alt[0] == "A" and insertion.pos == 19)
        # Test insertion multi nt and remove upstream and downstream
        insertion = VCFRecord(chrom, 18, None, "TGAT", ["TCGAGAT"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "-" and insertion.alt[0] == "CGA" and insertion.pos == 19)
        # Test insertion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord(chrom, 18, None, "TCtgATTAGC", ["TaGGctTATGCGC"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "CTGATTA" and insertion.alt[0] == "AGGCTTATGC" and insertion.pos == 19)
        # Test deletion
        insertion = VCFRecord(chrom, 18, None, "TA", ["T"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "A" and insertion.alt[0] == "-" and insertion.pos == 19)
        # Test deletion multi nt and remove upstream and downstream
        insertion = VCFRecord(chrom, 18, None, "TCGAGAT", ["TGAT"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "CGA" and insertion.alt[0] == "-" and insertion.pos == 19)
        # Test deletion multi nt with possible remove and downstream and without complete standardization
        insertion = VCFRecord(chrom, 18, None, "TaGGctTATGCGC", ["TCtgATTAGC"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.ref == "AGGCTTATGC" and insertion.alt[0] == "CTGATTA" and insertion.pos == 19)
        # Test insertion on start already normalized
        insertion = VCFRecord(chrom, 1, None, "-", ["A"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.pos == 1 and insertion.ref == "-" and insertion.alt[0] == "A")
        # Test insertion on start
        insertion = VCFRecord(chrom, 1, None, "G", ["AG"])
        insertion.normalizeSingleAllele()
        self.assertTrue(insertion.pos == 1 and insertion.ref == "-" and insertion.alt[0] == "A")
        # Test deletion on start already normalized
        deletion = VCFRecord(chrom, 1, None, "G", ["-"])
        deletion.normalizeSingleAllele()
        self.assertTrue(deletion.pos == 1 and deletion.ref == "G" and deletion.alt[0] == "-")
        # Test deletion on start
        deletion = VCFRecord(chrom, 1, None, "GN", ["N"])
        deletion.normalizeSingleAllele()
        self.assertTrue(deletion.pos == 1 and deletion.ref == "G" and deletion.alt[0] == "-")

    def testUpstreamed(self):
        chrom = "artificial_1"
        ref = "gnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test substit
            substit = VCFRecord(chrom, 18, None, "T", ["A"])
            upstream = substit.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "T" and upstream.alt[0] == "A")
            # Test large substit
            substit = VCFRecord(chrom, 18, None, "TTT", ["AGC"])
            upstream = substit.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "TTT" and upstream.alt[0] == "AGC")
            # Test fix deletion
            deletion = VCFRecord(chrom, 18, None, "TtTaAGC", ["T"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "TTTAAGC" and upstream.alt[0] == "T")
            deletion = VCFRecord(chrom, 19, None, "tTaAGC", ["-"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "TTTAAGC" and upstream.alt[0] == "T")
            # Test homopolymer deletion
            deletion = VCFRecord(chrom, 18, None, "TTT", ["T"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 16 and upstream.ref == "GTT" and upstream.alt[0] == "G")
            deletion = VCFRecord(chrom, 19, None, "TT", ["-"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 16 and upstream.ref == "GTT" and upstream.alt[0] == "G")
            # Test complex deletion
            deletion = VCFRecord(chrom, 25, None, "CGAgCC", ["C"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 21 and upstream.ref == "AAGCCG" and upstream.alt[0] == "A")
            deletion = VCFRecord(chrom, 26, None, "GAgCC", ["-"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 21 and upstream.ref == "AAGCCG" and upstream.alt[0] == "A")
            # Test deletion on start
            deletion = VCFRecord(chrom, 1, None, ref[0], ["-"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "GN" and upstream.alt[0] == "N")
            deletion = VCFRecord(chrom, 1, None, "GN", ["N"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "GN" and upstream.alt[0] == "N")
            # Test deletion on end
            deletion = VCFRecord(chrom, len(ref) - 1, None, ref[-2:], [ref[-2]])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 32 and upstream.ref == "AT" and upstream.alt[0] == "A")
            deletion = VCFRecord(chrom, len(ref), None, ref[-1:], ["-"])
            upstream = deletion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 32 and upstream.ref == "AT" and upstream.alt[0] == "A")
            # Test fix insertion
            insertion = VCFRecord(chrom, 4, None, "N", ["NGGTT"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 4 and upstream.ref == "N" and upstream.alt[0] == "NGGTT")
            insertion = VCFRecord(chrom, 5, None, "-", ["GGTT"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 4 and upstream.ref == "N" and upstream.alt[0] == "NGGTT")
            # Test homopolymer insertion
            insertion = VCFRecord(chrom, 18, None, "T", ["TT"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 16 and upstream.ref == "G" and upstream.alt[0] == "GT")
            insertion = VCFRecord(chrom, 19, None, "-", ["T"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 16 and upstream.ref == "G" and upstream.alt[0] == "GT")
            # Test complex insertion
            insertion = VCFRecord(chrom, 25, None, "C", ["CGAgCC"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 21 and upstream.ref == "A" and upstream.alt[0] == "AAGCCG")
            insertion = VCFRecord(chrom, 26, None, "-", ["GAgCC"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 21 and upstream.ref == "A" and upstream.alt[0] == "AAGCCG")
            # Test insertion on start
            insertion = VCFRecord(chrom, 1, None, "-", ["A"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "G" and upstream.alt[0] == "AG")
            insertion = VCFRecord(chrom, 1, None, "G", ["AG"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "G" and upstream.alt[0] == "AG")
            # Test insertion on end
            insertion = VCFRecord(chrom, len(ref), None, ref[-1], [ref[-1] + "G"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "TG")
            insertion = VCFRecord(chrom, len(ref) + 1, None, "-", ["G"])
            upstream = insertion.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "TG")
        chrom = "artificial_2"
        ref = "CCGTGATGCAT"  # length: 11
        #      | | | | | |
        #      1 3 5 7 9 11
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test insertion
            substit = VCFRecord(chrom, 8, None, "G", ["GATG"])
            upstream = substit.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 3 and upstream.ref == "G" and upstream.alt[0] == "GTGA")
            # Test deletion
            substit = VCFRecord(chrom, 5, None, "GATG", ["G"])
            upstream = substit.upstreamed(seq_handler)
            self.assertTrue(upstream.pos == 3 and upstream.ref == "GTGA" and upstream.alt[0] == "G")

    def testGetMostUpstream(self):
        chrom = "artificial_1"
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test substit
            substit = VCFRecord(chrom, 18, None, "T", ["A"])
            upstream = substit.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "T" and upstream.alt[0] == "A")
            # Test large substit
            substit = VCFRecord(chrom, 18, None, "TTT", ["AGC"])
            upstream = substit.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 18 and upstream.ref == "TTT" and upstream.alt[0] == "AGC")
            # Test fix deletion
            deletion = VCFRecord(chrom, 18, None, "TtTaAGC", ["T"], 230)
            upstream = deletion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 19 and upstream.ref == "TTAAGC" and upstream.alt[0] == "-")
            # Test homopolymer deletion
            deletion = VCFRecord(chrom, 18, None, "TTT", ["T"], 230)
            upstream = deletion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 17 and upstream.ref == "TT" and upstream.alt[0] == "-")
            # Test complex deletion
            deletion = VCFRecord(chrom, 25, None, "CGAgCC", ["C"], 230)
            upstream = deletion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 22 and upstream.ref == "AGCCG" and upstream.alt[0] == "-")
            # Test deletion on start
            deletion = VCFRecord(chrom, 1, None, ref[0], ["-"], 230)
            upstream = deletion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "N" and upstream.alt[0] == "-")
            # Test deletion on end
            deletion = VCFRecord(chrom, len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
            upstream = deletion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 33 and upstream.ref == "T" and upstream.alt[0] == "-")
            # Test fix insertion
            insertion = VCFRecord(chrom, 4, None, "N", ["NGGTT"], 230)
            upstream = insertion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 5 and upstream.ref == "-" and upstream.alt[0] == "GGTT")
            # Test homopolymer insertion
            insertion = VCFRecord(chrom, 18, None, "T", ["TT"], 230)
            upstream = insertion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 17 and upstream.ref == "-" and upstream.alt[0] == "T")
            # Test complex insertion
            insertion = VCFRecord(chrom, 25, None, "C", ["CGAgCC"], 230)
            upstream = insertion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 22 and upstream.ref == "-" and upstream.alt[0] == "AGCCG")
            # Test insertion on start
            insertion = VCFRecord(chrom, 1, None, "-", ["A"], 230)
            upstream = insertion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 1 and upstream.ref == "-" and upstream.alt[0] == "A")
            # Test insertion on end
            insertion = VCFRecord(chrom, len(ref), None, ref[-1], [ref[-1] + "G"], 230)
            upstream = insertion.getMostUpstream(seq_handler)
            self.assertTrue(upstream.pos == 34 and upstream.ref == "-" and upstream.alt[0] == "G")

    def testGetMostDownstream(self):
        chrom = "artificial_1"
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test substit
            substit = VCFRecord(chrom, 18, None, "T", ["A"])
            downstream = substit.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "T" and downstream.alt[0] == "A")
            # Test large substit
            substit = VCFRecord(chrom, 18, None, "TTT", ["AGC"])
            downstream = substit.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTT" and downstream.alt[0] == "AGC")
            # Test fix deletion
            deletion = VCFRecord(chrom, 18, None, "TtTaAGC", ["T"], 230)
            downstream = deletion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 19 and downstream.ref == "TTAAGC" and downstream.alt[0] == "-")
            # Test homopolymer deletion
            deletion = VCFRecord(chrom, 18, None, "TTT", ["T"], 230)
            downstream = deletion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 19 and downstream.ref == "TT" and downstream.alt[0] == "-")
            # Test complex deletion
            deletion = VCFRecord(chrom, 25, None, "CGAgCC", ["C"], 230)
            downstream = deletion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 28 and downstream.ref == "GCCGA" and downstream.alt[0] == "-")
            # Test deletion on start
            deletion = VCFRecord(chrom, 1, None, ref[0], ["-"], 230)
            downstream = deletion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "-")
            # Test deletion on end
            deletion = VCFRecord(chrom, len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
            downstream = deletion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "-")
            # Test fix insertion
            insertion = VCFRecord(chrom, 4, None, "N", ["NGGTT"], 230)
            downstream = insertion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 5 and downstream.ref == "-" and downstream.alt[0] == "GGTT")
            # Test homopolymer insertion
            insertion = VCFRecord(chrom, 18, None, "T", ["TT"], 230)
            downstream = insertion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 21 and downstream.ref == "-" and downstream.alt[0] == "T")
            # Test complex insertion
            insertion = VCFRecord(chrom, 25, None, "C", ["CGAgCC"], 230)
            downstream = insertion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 33 and downstream.ref == "-" and downstream.alt[0] == "GCCGA")
            # Test insertion on start
            insertion = VCFRecord(chrom, 1, None, "-", ["A"], 230)
            downstream = insertion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 1 and downstream.ref == "-" and downstream.alt[0] == "A")
            # Test insertion on end
            insertion = VCFRecord(chrom, len(ref), None, ref[-1], [ref[-1] + "G"], 230)
            downstream = insertion.getMostDownstream(seq_handler)
            self.assertTrue(downstream.pos == 34 and downstream.ref == "-" and downstream.alt[0] == "G")

    def testDownstreamed(self):
        chrom = "artificial_1"
        ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test substit
            substit = VCFRecord(chrom, 18, None, "T", ["A"])
            downstream = substit.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "T" and downstream.alt[0] == "A")
            # Test large substit
            substit = VCFRecord(chrom, 18, None, "TTT", ["AGC"])
            downstream = substit.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTT" and downstream.alt[0] == "AGC")
            # Test fix deletion
            deletion = VCFRecord(chrom, 18, None, "TTTaAGC", ["T"], 230)
            downstream = deletion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTTAAGC" and downstream.alt[0] == "T")
            # Test homopolymer deletion
            deletion = VCFRecord(chrom, 16, None, "GTT", ["G"], 230)
            downstream = deletion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 18 and downstream.ref == "TTT" and downstream.alt[0] == "T")
            # Test complex deletion
            deletion = VCFRecord(chrom, 25, None, "CGAgCC", ["C"], 230)
            downstream = deletion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 27 and downstream.ref == "AGCCGA" and downstream.alt[0] == "A")
            # Test deletion on start
            deletion = VCFRecord(chrom, 1, None, ref[0], ["-"], 230)
            downstream = deletion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 3 and downstream.ref == "NN" and downstream.alt[0] == "N")
            # Test deletion on end
            deletion = VCFRecord(chrom, len(ref) - 1, None, ref[-2:], [ref[-2]], 230)
            downstream = deletion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 32 and downstream.ref == "AT" and downstream.alt[0] == "A")
            # Test fix insertion
            insertion = VCFRecord(chrom, 4, None, "N", ["NGGTT"], 230)
            downstream = insertion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 4 and downstream.ref == "N" and downstream.alt[0] == "NGGTT")
            # Test homopolymer insertion
            insertion = VCFRecord(chrom, 18, None, "T", ["TT"], 230)
            downstream = insertion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 20 and downstream.ref == "T" and downstream.alt[0] == "TT")
            # Test complex insertion
            insertion = VCFRecord(chrom, 25, None, "C", ["CGAgCC"], 230)
            downstream = insertion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 32 and downstream.ref == "A" and downstream.alt[0] == "AGCCGA")
            # Test insertion on end
            insertion = VCFRecord(chrom, len(ref), None, ref[-1], [ref[-1] + "G"], 230)
            downstream = insertion.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 33 and downstream.ref == "T" and downstream.alt[0] == "TG")
        chrom = "artificial_2"
        ref = "CCGTGATGCAT"  # length: 11
        #      | | | | | |
        #      1 3 5 7 9 11
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test insertion
            substit = VCFRecord(chrom, 3, None, "G", ["GTGA"])
            downstream = substit.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 8 and downstream.ref == "G" and downstream.alt[0] == "GATG")
            # Test deletion
            substit = VCFRecord(chrom, 3, None, "GTGA", ["G"])
            downstream = substit.downstreamed(seq_handler)
            self.assertTrue(downstream.pos == 5 and downstream.ref == "GATG" and downstream.alt[0] == "G")

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

    def testStandardize(self):
        chrom = "artificial_1"
        ref = "gnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #      | | | | | | | | | | | | | | | | |
        #      1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                  13  17  21  25  29  33
        fastIdxRef(chrom, ref, self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Test substit
            substit = VCFRecord(chrom, 18, None, "T", ["A"])
            substit.standardize(seq_handler)
            self.assertTrue(substit.pos == 18 and substit.ref == "T" and substit.alt[0] == "A")
            # Test substit with useless prev
            substitution = VCFRecord(chrom, 20, None, "TA", ["TC"])
            substitution.normalizeSingleAllele()
            self.assertTrue(substitution.pos == 21 and substitution.ref == "A" and substitution.alt[0] == "C")
            # Test substit with useless prev and next
            substitution = VCFRecord(chrom, 20, None, "TAA", ["TCA"])
            substitution.normalizeSingleAllele()
            self.assertTrue(substitution.pos == 21 and substitution.ref == "A" and substitution.alt[0] == "C")
            # Test large substit
            substit = VCFRecord(chrom, 18, None, "TTT", ["AGC"])
            substit.standardize(seq_handler)
            self.assertTrue(substit.pos == 18 and substit.ref == "TTT" and substit.alt[0] == "AGC")
            # Test fix deletion
            deletion = VCFRecord(chrom, 18, None, "TtTaAGC", ["T"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 18 and deletion.ref == "TTTAAGC" and deletion.alt[0] == "T")
            deletion = VCFRecord(chrom, 19, None, "tTaAGC", ["-"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 18 and deletion.ref == "TTTAAGC" and deletion.alt[0] == "T")
            # Test homopolymer deletion
            deletion = VCFRecord(chrom, 18, None, "TTT", ["T"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 16 and deletion.ref == "GTT" and deletion.alt[0] == "G")
            deletion = VCFRecord(chrom, 19, None, "TT", ["-"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 16 and deletion.ref == "GTT" and deletion.alt[0] == "G")
            # Test complex deletion
            deletion = VCFRecord(chrom, 25, None, "CGAgCC", ["C"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 21 and deletion.ref == "AAGCCG" and deletion.alt[0] == "A")
            deletion = VCFRecord(chrom, 26, None, "GAgCC", ["-"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 21 and deletion.ref == "AAGCCG" and deletion.alt[0] == "A")
            # Test deletion on start
            deletion = VCFRecord(chrom, 1, None, ref[0], ["-"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 1 and deletion.ref == "GN" and deletion.alt[0] == "N")
            deletion = VCFRecord(chrom, 1, None, "GN", ["N"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 1 and deletion.ref == "GN" and deletion.alt[0] == "N")
            # Test deletion on end
            deletion = VCFRecord(chrom, len(ref) - 1, None, ref[-2:], [ref[-2]])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 32 and deletion.ref == "AT" and deletion.alt[0] == "A")
            deletion = VCFRecord(chrom, len(ref), None, ref[-1:], ["-"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 32 and deletion.ref == "AT" and deletion.alt[0] == "A")
            # Test deletion multi nt with useless up and down and without complete normalization
            deletion = VCFRecord(chrom, 18, None, "TtTaAGCCGAG", ["TTGGAG"])
            deletion.standardize(seq_handler)
            self.assertTrue(deletion.pos == 20 and deletion.ref == "TAAGCC" and deletion.alt[0] == "G")
            # Test fix insertion
            insertion = VCFRecord(chrom, 4, None, "N", ["NGGTT"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 4 and insertion.ref == "N" and insertion.alt[0] == "NGGTT")
            insertion = VCFRecord(chrom, 5, None, "-", ["GGTT"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 4 and insertion.ref == "N" and insertion.alt[0] == "NGGTT")
            # Test homopolymer insertion
            insertion = VCFRecord(chrom, 18, None, "T", ["TT"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 16 and insertion.ref == "G" and insertion.alt[0] == "GT")
            insertion = VCFRecord(chrom, 19, None, "-", ["T"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 16 and insertion.ref == "G" and insertion.alt[0] == "GT")
            # Test complex insertion
            insertion = VCFRecord(chrom, 25, None, "C", ["CGAgCC"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 21 and insertion.ref == "A" and insertion.alt[0] == "AAGCCG")
            insertion = VCFRecord(chrom, 26, None, "-", ["GAgCC"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 21 and insertion.ref == "A" and insertion.alt[0] == "AAGCCG")
            # Test insertion multi nt and remove upstream and downstream
            insertion = VCFRecord(chrom, 11, None, "ATGAT", ["ATCGAGAT"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 12 and insertion.ref == "T" and insertion.alt[0] == "TCGA")
            # Test insertion on start
            insertion = VCFRecord(chrom, 1, None, "-", ["A"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 1 and insertion.ref == "G" and insertion.alt[0] == "AG")
            insertion = VCFRecord(chrom, 1, None, "G", ["AG"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 1 and insertion.ref == "G" and insertion.alt[0] == "AG")
            # Test insertion on end
            insertion = VCFRecord(chrom, len(ref), None, ref[-1], [ref[-1] + "G"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 33 and insertion.ref == "T" and insertion.alt[0] == "TG")
            insertion = VCFRecord(chrom, len(ref) + 1, None, "-", ["G"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 33 and insertion.ref == "T" and insertion.alt[0] == "TG")
            # Test insertion multi nt with useless up and down and without complete normalization
            insertion = VCFRecord(chrom, 11, None, "aTgATGTT", ["aTTTCAGTT"])
            insertion.standardize(seq_handler)
            self.assertTrue(insertion.pos == 13 and insertion.ref == "GAT" and insertion.alt[0] == "TTCA")

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

    def testType(self):
        data = [
            {"variant": VCFRecord("chr1", 2, None, "A", ["T"]), "expected": "snp"},
            {"variant": VCFRecord("chr1", 2, None, "A", ["AT"]), "expected": "indel"},
            {"variant": VCFRecord("chr1", 1, None, "CA", ["C"]), "expected": "indel"},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GT"]), "expected": "variation"},
            {"variant": VCFRecord("chr1", 15, None, "CA", ["GTT"]), "expected": "indel"},
            {"variant": VCFRecord("chr1", 15, None, "CAC", ["GT"]), "expected": "indel"}
        ]
        for curr in data:
            self.assertEqual(curr["variant"].type(), curr["expected"])


class TestVCFSymbAltRecord(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_fa = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + ".fa.fai")
        with open(self.tmp_variants, "w") as writer:
            writer.write("""##fileformat=VCFv4.3
##fileDate=20100501
##reference=1000GenomesPilot-NCBI36
##assembly=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta
##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001
1	100	cnv_notation	T	<CNV:TR>	.	.	END=130;SVLEN=30	GT	1/2
1	2827694	rs2376870	CGTGGATGCGGGGAC	C	.	PASS	SVTYPE=DEL;END=2827708;HOMLEN=1;HOMSEQ=G;SVLEN=-14	GT:GQ	1/1:14
2	321682	.	T	<DEL>	6	PASS	SVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62	GT:GQ	0/1:12
2	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV;END=421681	.	.
2	14477084	.	C	<DEL:ME:ALU>	12	PASS	SVTYPE=DEL;END=14477381;SVLEN=-297;CIPOS=-22,18;CIEND=-12,32	GT:GQ	0/1:12
3	9425916	.	C	<INS:ME:L1>	23	PASS	SVTYPE=INS;END=9425916;SVLEN=6027;CIPOS=-16,22	GT:GQ	1/1:15
3	12665100	.	A	<DUP>	14	PASS	SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-500,500;CIEND=-500,500	GT:GQ:CN:CNQ	./.:0:3:16.2
4	18665128	.	T	<DUP:TANDEM>	11	PASS	SVTYPE=DUP;END=18665204;SVLEN=76;CIPOS=-10,10;CIEND=-10,10	GT:GQ:CN:CNQ	./.:0:5:8.3
17	41242722	.	N	<DEL>	.	.	END=41246879	GT	./.""")

    def tearDown(self):
        for curr_file in [self.tmp_variants, self.tmp_fa, self.tmp_fai]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testDownstreamed(self):
        chrom = "chrA"
        fastIdxRef(chrom, "ATGCGAAAAAAATGTCCGTGATGCAT", self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # IMPRECISE
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "IMPRECISE": True})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (var.pos, var.ref, var.alt[0])
            )
            # DEL with CIPOS
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DEL>"], info={"SVLEN": [2], "CIPOS": (0, 5)})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"]),
                (10, "A", "<DEL>", (-5, 0))
            )
            # DEL without CIPOS
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DEL>"], info={"SVLEN": [2]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (10, "A", "<DEL>")
            )
            # DEL with CIPOS and END
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DEL>"], info={"SVLEN": [2], "CIPOS": (0, 5), "END": 7})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"], down_var.info["END"]),
                (10, "A", "<DEL>", (-5, 0), 12)
            )
            # DEL without CIPOS and END
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DEL>"], info={"SVLEN": [2], "END": 7})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["END"]),
                (10, "A", "<DEL>", 12)
            )
            # DEL complex
            var = VCFSymbAltRecord(chrom, 18, None, "G", ["<DEL>"], info={"SVLEN": [3]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (20, "G", "<DEL>")
            )
            # DUP with CIPOS
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-1, 4)})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"]),
                (10, "A", "<DUP>", (-5, 0))
            )
            # DUP without CIPOS
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<DUP>"], info={"SVLEN": [2]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (10, "A", "<DUP>")
            )
            # DUP without CIPOS and with convert to INS
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<DUP>"], info={"SVLEN": [2]})
            down_var = var.downstreamed(seq_handler, convert_dup=True)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (12, "A", "<INS>")
            )
            # DUP with CIPOS and END
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-1, 4), "END": 8})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"], down_var.info["END"]),
                (10, "A", "<DUP>", (-5, 0), 12)
            )
            # DUP with CIPOS and END and convert to INS
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-1, 4), "END": 8})
            down_var = var.downstreamed(seq_handler, convert_dup=True)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"], down_var.info["END"]),
                (12, "A", "<INS>", (-7, 0), 12)
            )
            # DUP without CIPOS and END
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DUP>"], info={"SVLEN": [2], "END": 7})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["END"]),
                (10, "A", "<DUP>", 12)
            )
            # DUP complex move
            var = VCFSymbAltRecord(chrom, 18, None, "G", ["<DUP>"], info={"SVLEN": [3]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (20, "G", "<DUP>")
            )
            # DUP complex move to INS
            var = VCFSymbAltRecord(chrom, 18, None, "G", ["<DUP>"], info={"SVLEN": [3]})
            down_var = var.downstreamed(seq_handler, convert_dup=True)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (23, "G", "<INS>")
            )
            # INV with END
            var = VCFSymbAltRecord(chrom, 3, None, "G", ["<INV>"], info={"END": 11})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["END"]),
                (3, "G", "<INV>", 11)
            )
            # INS with CIPOS
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<INS>"], info={"CIPOS": [-1,4]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"]),
                (10, "A", "<INS>", (-5, 0))
            )
            # INS without CIPOS => cannot be moved
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={})
            down_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (9, "A", "<INS>")
            )
            # INS with CIPOS and END
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<INS>"], info={"CIPOS": (-1, 4), "END": 6})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"], down_var.info["END"]),
                (10, "A", "<INS>", (-5, 0), 10)
            )
            # INS without CIPOS and with END => cannot be moved
            var = VCFSymbAltRecord(chrom, 6, None, "A", ["<INS>"], info={"END": 6})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["END"]),
                (6, "A", "<INS>", 6)
            )
            # CNV with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1)})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"]),
                (10, "A", "<CNV:TR>", (-5, 0))
            )
            # CNV without CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2]})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (10, "A", "<CNV:TR>")
            )
            # CNV with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1), "END": 11})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["CIPOS"], down_var.info["END"]),
                (10, "A", "<CNV:TR>", (-5, 0), 12)
            )
            # CNV without CIPOS and with END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "END": 11})
            down_var = var.downstreamed(seq_handler)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0], down_var.info["END"]),
                (10, "A", "<CNV:TR>", 12)
            )
            # CNV complex move
            var = VCFSymbAltRecord(chrom, 18, None, "G", ["<CNV>"], info={"SVLEN": [3]})
            down_var = var.downstreamed(seq_handler, convert_dup=True)
            self.assertEqual(
                (down_var.pos, down_var.ref, down_var.alt[0]),
                (20, "G", "<CNV>")
            )

    def testEnd(self):
        # With END
        expected = [130, 321887, 421681, 14477381, 9425916, 12686200, 18665204, 41246879]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.__class__ is VCFSymbAltRecord:
                    observed.append(record.end)
        self.assertEqual(observed, expected)
        # Without END
        expected = [130, 321887, None, 14477381, 9425916, 12686200, 18665204, None]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.__class__ is VCFSymbAltRecord:
                    if "END" in record.info:
                        del record.info["END"]
                    observed.append(record.end)
        self.assertEqual(observed, expected)

    def testGetName(self):
        expected = [
            "1:100[30]=T/<CNV:TR>",
            "1:2827694=CGTGGATGCGGGGAC/C",
            "2:321682[-205]=T/<DEL>",
            "2:321682[99999]=T/<INV>",
            "2:14477084[-297]=C/<DEL:ME:ALU>",
            "3:9425916[6027]=C/<INS:ME:L1>",
            "3:12665100[21100]=A/<DUP>",
            "4:18665128[76]=T/<DUP:TANDEM>",
            "17:41242722[-4157]=N/<DEL>"
        ]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                observed.append(record.getName())
        self.assertEqual(observed, expected)

    def testIsDeletion(self):
        expected = [True, True, False, True, False, False, False, True]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.alt[0].startswith("<CNV"):
                    with self.assertRaises(Exception) as context:
                        record.isDeletion()
                    self.assertTrue('CNV can be deletion or insertion' in str(context.exception))
                else:
                    observed.append(record.isDeletion())
        self.assertEqual(observed, expected)

    def testIsIndel(self):
        expected = [True, True, True, False, True, True, True, True, True]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                observed.append(record.isIndel())
        self.assertEqual(observed, expected)

    def testIsInsAndDel(self):
        expected = [False, False, False, True, False, False, False, False, False]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                observed.append(record.isInsAndDel())
        self.assertEqual(observed, expected)

    def testIsInsertion(self):
        expected = [False, False, False, False, True, True, True, False]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.alt[0].startswith("<CNV"):
                    with self.assertRaises(Exception) as context:
                        record.isInsertion()
                    self.assertTrue('CNV can be deletion or insertion' in str(context.exception))
                else:
                    observed.append(record.isInsertion())
        self.assertEqual(observed, expected)

    def testIsInversion(self):
        expected = [False, False, True, False, False, False, False, False]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.__class__ is VCFSymbAltRecord:
                    observed.append(record.isInversion())
        self.assertEqual(observed, expected)

    def testLitteralized(self):
        chrom = "chrA"
        fastIdxRef(chrom, "ATGCGAAAAAAATGT", self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # DEL
            var = VCFSymbAltRecord(chrom, 2, None, "T", ["<DEL>"], info={"SVLEN": [2], "END": 4})
            litteral_var = var.litteralized(seq_handler)
            self.assertEqual(
                (litteral_var.pos, litteral_var.ref, litteral_var.alt[0]),
                (2, "TGC", "T")
            )
            # DUP in homopolymer
            var = VCFSymbAltRecord(chrom, 5, None, "G", ["<DUP>"], info={"SVLEN": [3], "CIPOS": [0, 5]})
            litteral_var = var.litteralized(seq_handler)
            self.assertEqual(
                (litteral_var.pos, litteral_var.ref, litteral_var.alt[0]),
                (5, "G", "GAAA")
            )
            # INV
            var = VCFSymbAltRecord(chrom, 3, None, "G", ["<INV>"], info={"END": 11})
            litteral_var = var.litteralized(seq_handler)
            self.assertEqual(
                (litteral_var.pos, litteral_var.ref, litteral_var.alt[0]),
                (4, "CGAAAAAA", "AAAAAAGC")
            )
            # INS
            var = VCFSymbAltRecord(chrom, 14, None, "T", ["<INS>"], info={"SVLEN": [100]})
            with self.assertRaises(Exception) as context:
                var.litteralized(seq_handler)
            self.assertTrue('sequence is unknown' in str(context.exception))
            # IMPRECISE
            var = VCFSymbAltRecord(chrom, 2, None, "T", ["<DEL>"], info={"IMPRECISE": True})
            with self.assertRaises(Exception) as context:
                var.litteralized(seq_handler)
            self.assertTrue('real position is unknown' in str(context.exception))
            # CNV
            var = VCFSymbAltRecord(chrom, 14, None, "T", ["<CNV>"], info={"SVLEN": [100]})
            with self.assertRaises(Exception) as context:
                var.litteralized(seq_handler)
            self.assertTrue('CNV can be DUP or DEL' in str(context.exception))

    def testParsingClass(self):
        expected = [
            VCFSymbAltRecord, VCFRecord, VCFSymbAltRecord, VCFSymbAltRecord,
            VCFSymbAltRecord, VCFSymbAltRecord, VCFSymbAltRecord,
            VCFSymbAltRecord, VCFSymbAltRecord
        ]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                observed.append(record.__class__)
        self.assertEqual(observed, expected)

    def testRefEnd(self):
        expected = [2827708, 321887, 421681, 14477381, 9425916.5, 12665100.5, 18665128.5, 41246879]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.alt[0].startswith("<CNV"):
                    with self.assertRaises(Exception) as context:
                        record.refEnd()
                    self.assertTrue('it can be deletion or insertion' in str(context.exception))
                else:
                    observed.append(record.refEnd())
        self.assertEqual(observed, expected)

    def testRefStart(self):
        expected = [2827695, 321683, 321683, 14477085, 9425916.5, 12665100.5, 18665128.5, 41242723]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.alt[0].startswith("<CNV"):
                    with self.assertRaises(Exception) as context:
                        record.refStart()
                    self.assertTrue('it can be deletion or insertion' in str(context.exception))
                else:
                    observed.append(record.refStart())
        self.assertEqual(observed, expected)

    def testStandardize(self):
        chrom = "chrA"
        fastIdxRef(chrom, "ATGCGAAAAAAATGTCCGTGATGCAT", self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # IMPRECISE
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "IMPRECISE": True})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (10, "A", "<DEL>")
            )
            # DEL with CIPOS
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "CIPOS": [-5,0]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"]),
                (5, "G", "<DEL>", (0, 5))
            )
            # DEL without CIPOS
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (5, "G", "<DEL>")
            )
            # DEL with CIPOS and END
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "CIPOS": (-5, 0), "END": 12})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"], var.info["END"]),
                (5, "G", "<DEL>", (0, 5), 7)
            )
            # DEL without CIPOS and END
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "END": 12})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["END"]),
                (5, "G", "<DEL>", 7)
            )
            # DEL complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<DEL>"], info={"SVLEN": [3]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (18, "G", "<DEL>")
            )
            # DUP with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-4, 1)})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"]),
                (5, "G", "<DUP>", (0, 5))
            )
            # DUP without CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (5, "G", "<DUP>")
            )
            # DUP with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-4, 1), "END": 11})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"], var.info["END"]),
                (5, "G", "<DUP>", (0, 5), 7)
            )
            # DUP without CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "END": 11})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["END"]),
                (5, "G", "<DUP>", 7)
            )
            # DUP complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<DUP>"], info={"SVLEN": [3]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (18, "G", "<DUP>")
            )
            # INV with END
            var = VCFSymbAltRecord(chrom, 3, None, "G", ["<INV>"], info={"END": 11})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["END"]),
                (3, "G", "<INV>", 11)
            )
            # INS with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"CIPOS": [-4,1]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"]),
                (5, "G", "<INS>", (0, 5))
            )
            # INS without CIPOS => cannot be moved
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (9, "A", "<INS>")
            )
            # INS with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"CIPOS": (-4, 1), "END": 9})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"], var.info["END"]),
                (5, "G", "<INS>", (0, 5), 5)
            )
            # INS without CIPOS and with END => cannot be moved
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"END": 9})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["END"]),
                (9, "A", "<INS>", 9)
            )
            # CNV with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1)})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"]),
                (5, "G", "<CNV:TR>", (0, 5))
            )
            # CNV without CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (5, "G", "<CNV:TR>")
            )
            # CNV with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1), "END": 11})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["CIPOS"], var.info["END"]),
                (5, "G", "<CNV:TR>", (0, 5), 7)
            )
            # CNV without CIPOS and with END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "END": 11})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0], var.info["END"]),
                (5, "G", "<CNV:TR>", 7)
            )
            # CNV complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<CNV>"], info={"SVLEN": [3]})
            var.standardize(seq_handler)
            self.assertEqual(
                (var.pos, var.ref, var.alt[0]),
                (18, "G", "<CNV>")
            )

    def testSvLen(self):
        # With SVLEN
        expected = [30, -205, 99999, -297, 6027, 21100, 76, -4157]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.__class__ is VCFSymbAltRecord:
                    observed.append(record.sv_len)
        self.assertEqual(observed, expected)
        # Without SVLEN
        expected = [30, -205, 99999, -297, None, 21100, 76, -4157]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                if record.__class__ is VCFSymbAltRecord:
                    if "SVLEN" in record.info:
                        del record.info["SVLEN"]
                    observed.append(record.sv_len)
        self.assertEqual(observed, expected)

    def testType(self):
        expected = ["indel", "indel", "indel", "inv", "indel", "indel", "indel", "indel", "indel"]
        observed = list()
        with VCFIO(self.tmp_variants) as reader:
            for record in reader:
                observed.append(record.type())
        self.assertEqual(observed, expected)

    def testUpstreamed(self):
        chrom = "chrA"
        fastIdxRef(chrom, "ATGCGAAAAAAATGTCCGTGATGCAT", self.tmp_fa, self.tmp_fai)
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # IMPRECISE
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "IMPRECISE": True})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (var.pos, var.ref, var.alt[0])
            )
            # DEL with CIPOS
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "CIPOS": [-5,0]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"]),
                (5, "G", "<DEL>", (0, 5))
            )
            # DEL without CIPOS
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (5, "G", "<DEL>")
            )
            # DEL with CIPOS and END
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "CIPOS": (-5, 0), "END": 12})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"], up_var.info["END"]),
                (5, "G", "<DEL>", (0, 5), 7)
            )
            # DEL without CIPOS and END
            var = VCFSymbAltRecord(chrom, 10, None, "A", ["<DEL>"], info={"SVLEN": [2], "END": 12})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["END"]),
                (5, "G", "<DEL>", 7)
            )
            # DEL complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<DEL>"], info={"SVLEN": [3]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (18, "G", "<DEL>")
            )
            # DUP with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-4, 1)})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"]),
                (5, "G", "<DUP>", (0, 5))
            )
            # DUP without CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (5, "G", "<DUP>")
            )
            # DUP with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "CIPOS": (-4, 1), "END": 11})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"], up_var.info["END"]),
                (5, "G", "<DUP>", (0, 5), 7)
            )
            # DUP without CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<DUP>"], info={"SVLEN": [2], "END": 11})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["END"]),
                (5, "G", "<DUP>", 7)
            )
            # DUP complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<DUP>"], info={"SVLEN": [3]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (18, "G", "<DUP>")
            )
            # INV with END
            var = VCFSymbAltRecord(chrom, 3, None, "G", ["<INV>"], info={"END": 11})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["END"]),
                (3, "G", "<INV>", 11)
            )
            # INS with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"CIPOS": [-4,1]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"]),
                (5, "G", "<INS>", (0, 5))
            )
            # INS without CIPOS => cannot be moved
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (9, "A", "<INS>")
            )
            # INS with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"CIPOS": (-4, 1), "END": 9})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"], up_var.info["END"]),
                (5, "G", "<INS>", (0, 5), 5)
            )
            # INS without CIPOS and with END => cannot be moved
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<INS>"], info={"END": 9})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["END"]),
                (9, "A", "<INS>", 9)
            )
            # CNV with CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1)})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"]),
                (5, "G", "<CNV:TR>", (0, 5))
            )
            # CNV without CIPOS
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (5, "G", "<CNV:TR>")
            )
            # CNV with CIPOS and END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "CIPOS": (-4, 1), "END": 11})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["CIPOS"], up_var.info["END"]),
                (5, "G", "<CNV:TR>", (0, 5), 7)
            )
            # CNV without CIPOS and with END
            var = VCFSymbAltRecord(chrom, 9, None, "A", ["<CNV:TR>"], info={"SVLEN": [2], "END": 11})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0], up_var.info["END"]),
                (5, "G", "<CNV:TR>", 7)
            )
            # CNV complex move
            var = VCFSymbAltRecord(chrom, 20, None, "G", ["<CNV>"], info={"SVLEN": [3]})
            up_var = var.upstreamed(seq_handler)
            self.assertEqual(
                (up_var.pos, up_var.ref, up_var.alt[0]),
                (18, "G", "<CNV>")
            )


class TestRefWalker(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_fa = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + ".fa.fai")
        self.chrom = "test"
        self.ref = "nnNNATGCCAaTgATGTTtTaAGCCGAGCCGAT"  # length: 33
        #           | | | | | | | | | | | | | | | | |
        #           1 3 5 7 9 11| 15| 19| 23| 27| 31|
        #                       13  17  21  25  29  33
        fastIdxRef(self.chrom, self.ref, self.tmp_fa, self.tmp_fai)

    def tearDown(self):
        for curr_file in [self.tmp_fa, self.tmp_fai]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testDownIter(self):
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Short window, Large window and Too large window
            expected = [elt.upper() for elt in self.ref[14 - 7:14]][::-1]
            for buffer_size in [2, 10, 100]:
                observed = list()
                for idx, curr in enumerate(_refDownIter(seq_handler, self.chrom, 14, buffer_size)):
                    if idx == 7:
                        break
                    observed.append(curr)
                self.assertEqual(observed, expected)
            # Stop iterration
            expected = [elt.upper() for elt in self.ref[0:14]][::-1]
            observed = list()
            for curr in _refDownIter(seq_handler, self.chrom, 14, 3):
                observed.append(curr)
            self.assertEqual(observed, expected)

    def testUpIter(self):
        with IdxFastaIO(self.tmp_fa) as seq_handler:
            # Short window, Large window and Too large window
            expected = [elt.upper() for elt in self.ref[14:14 + 7]]
            for buffer_size in [2, 10, 100]:
                observed = list()
                for idx, curr in enumerate(_refUpIter(seq_handler, self.chrom, 15, buffer_size)):
                    if idx == 7:
                        break
                    observed.append(curr)
                self.assertEqual(observed, expected)
            # Stop iterration
            expected = [elt.upper() for elt in self.ref[14:]]
            observed = list()
            for curr in _refUpIter(seq_handler, self.chrom, 15, 3):
                observed.append(curr)
            self.assertEqual(observed, expected)


def fastIdxRef(chrom, seq, seq_path, idx_path):
    with open(seq_path, "w") as writer:
        writer.write(">{}\n".format(chrom))
        writer.write(seq)
    with open(idx_path, "w") as writer:
        writer.write("{}\t{}\t{}\t{}\t{}".format(chrom, len(seq), len(chrom) + 2, len(seq), len(seq) + 1))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
