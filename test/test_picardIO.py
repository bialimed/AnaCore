#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.picardIO import castCol, getColType, getValType, PicardReader


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestCastCol(unittest.TestCase):
    def testCastCol(self):
        # To str
        data = ["a", "b", "c", "10", "", "15.1", "16.3"]
        expected = ["a", "b", "c", "10", None, "15.1", "16.3"]
        castCol(data, "str")
        self.assertEqual(expected, data)
        # To int
        data = ["10", "", "15", "16"]
        expected = [10, None, 15, 16]
        castCol(data, "int")
        self.assertEqual(expected, data)
        # To float
        data = ["10", "15.1", "", "16.3"]
        expected = [10.0, 15.1, None, 16.3]
        castCol(data, "float")
        self.assertEqual(expected, data)


class TestGetColType(unittest.TestCase):
    def testGetColType(self):
        # Is None
        data = [None, ""]
        expected = None
        observed = getColType(data)
        self.assertEqual(expected, observed)
        # Is str
        data = ["a", "b", "c", "10", None, "", "15.1", "16.3"]
        expected = "str"
        observed = getColType(data)
        self.assertEqual(expected, observed)
        # Is int
        data = ["10", None, "", "15", "16"]
        expected = "int"
        observed = getColType(data)
        self.assertEqual(expected, observed)
        # Is float
        data = ["10", None, "", "15.1", "16"]
        expected = "float"
        observed = getColType(data)
        self.assertEqual(expected, observed)


class TestGetValType(unittest.TestCase):
    def testGetValType(self):
        data = ["a", "b", "c", "10", None, "", "15.1", "16.3"]
        expected = ["str", "str", "str", "int", None, "str", "float", "float"]
        observed = [getValType(elt) for elt in data]
        self.assertEqual(expected, observed)


class TestPicardReader(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        # CollectArraysVariantCallingMetrics
        self.tmp_variants_call = os.path.join(tmp_folder, unique_id + "_vcm.txt")
        with open(self.tmp_variants_call, "w") as writer:
            writer.write("""## htsjdk.samtools.metrics.StringHeader
# CollectArraysVariantCallingMetrics  --INPUT /cromwell_root/broad-gotc-dev-cromwell-execution/Arrays/9972b43c-df1e-414a-9017-24edc37b5250/call-IlluminaGenotypingArray/IlluminaGenotypingArray/805b2e59-add2-4794-ae89-350d69a0f7c3/call-MergePedIntoVcf/7991775143_R01C01.vcf.gz --OUTPUT 7991775143_R01C01 --CALL_RATE_PF_THRESHOLD 0.98 --DBSNP /cromwell_root/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz  --NUM_PROCESSORS 0 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --GA4GH_CLIENT_SECRETS client_secrets.json --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
## htsjdk.samtools.metrics.StringHeader
# Started on: Thu Mar 05 18:21:35 UTC 2020

## METRICS CLASS	picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingDetailMetrics
CHIP_WELL_BARCODE	SAMPLE_ALIAS	ANALYSIS_VERSION	CHIP_TYPE	AUTOCALL_PF	AUTOCALL_DATE	IMAGING_DATE	IS_ZCALLED	AUTOCALL_GENDER	FP_GENDER	REPORTED_GENDER	GENDER_CONCORDANCE_PF	HET_PCT	CLUSTER_FILE_NAME	P95_GREEN	P95_RED	AUTOCALL_VERSION	ZCALL_VERSION	EXTENDED_MANIFEST_VERSION	HET_HOMVAR_RATIO	SCANNER_NAME	NUM_ASSAYS	NUM_NON_FILTERED_ASSAYS	NUM_FILTERED_ASSAYS	NUM_ZEROED_OUT_ASSAYS	NUM_SNPS	NUM_INDELS	NUM_CALLS	NUM_AUTOCALL_CALLS	NUM_NO_CALLS	NUM_IN_DB_SNP	NOVEL_SNPS	PCT_DBSNP	CALL_RATE	AUTOCALL_CALL_RATE	NUM_SINGLETONS
7991775143_R01C01	NA12878	1	HumanExome-12v1-1_A	Y	2019-12-14T03:05:00+0000	2012-12-16T16:10:10+0000	Y	F	U	F	Y	0.042993	HumanExomev1_1_CEPH_A.egt	4251	5479	3.0.0	1.0.0.0	1.3	1.525553	N370	242513	242493	20	0	242450	43	242318	241250	175	238085	4365	0.981996	0.999278	0.994874	10418

""")
        # CollectWgsMetrics
        self.tmp_wgs = os.path.join(tmp_folder, unique_id + "_wgsm.txt")
        with open(self.tmp_wgs, "w") as writer:
            writer.write("""## htsjdk.samtools.metrics.StringHeader
# CollectWgsMetrics INPUT=/cromwell_root/broad-gotc-dev-cromwell-execution/WholeGenomeGermlineSingleSample/e36be5d3-ddbf-4ee6-a98d-2a1933826550/call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/da7029e5-cf73-4d9b-8369-23c3aea72f19/call-GatherBamFiles/NA12878_PLUMBING.bam OUTPUT=NA12878_PLUMBING.wgs_metrics INCLUDE_BQ_HISTOGRAM=true INTERVALS=/cromwell_root/gcp-public-data--broad-references/hg38/v0/wgs_coverage_regions.hg38.interval_list USE_FAST_ALGORITHM=true READ_LENGTH=250 VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=/cromwell_root/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta    MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 COVERAGE_CAP=250 LOCUS_ACCUMULATION_CAP=100000 STOP_AFTER=-1 COUNT_UNPAIRED=false SAMPLE_SIZE=10000 ALLELE_FRACTION=[0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5] VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
## htsjdk.samtools.metrics.StringHeader
# Started on: Fri Mar 20 20:49:11 UTC 2020

## METRICS CLASS	picard.analysis.WgsMetrics
GENOME_TERRITORY	MEAN_COVERAGE	SD_COVERAGE	MEDIAN_COVERAGE	MAD_COVERAGE	PCT_EXC_ADAPTER	PCT_EXC_MAPQ	PCT_EXC_DUPE	PCT_EXC_UNPAIRED	PCT_EXC_BASEQ	PCT_EXC_OVERLAP	PCT_EXC_CAPPED	PCT_EXC_TOTAL	PCT_1X	PCT_5X	PCT_10X	PCT_15X	PCT_20X	PCT_25X	PCT_30X	PCT_40X	PCT_50X	PCT_60X	PCT_70X	PCT_80X	PCT_90X	PCT_100X	HET_SNP_SENSITIVITY	HET_SNP_Q
2745186692	0.003328	0.059988	0	0	0.000045	0.044363	0.000482	0.00674	0.121016	0.090623	0	0.263269	0.003305	0.000001	0	0	0	0	0	0	0	0	0	0	0	0	0.000105	0

## HISTOGRAM	java.lang.Integer
coverage	high_quality_coverage_count	unfiltered_baseq_count
0	2736113146	0
1	9041943	0
2	24393	0
3	3102	730
4	1433	940
5	779	10184
6	861	1526
7	236	0
8	149	0
9	83	0
10	75	45052
11	54	0
12	0	0
13	0	0
14	0	0
15	0	0
16	0	0
17	1	0
18	1	0
19	10	0
20	0	67107
21	2	0
22	0	0
23	3	0
24	11	0
25	7	0
26	16	0
27	1	0
28	13	0
29	17	0
30	23	10192669
31	12	0
32	22	0
""")
        # CollectInsertSizeMetrics
        self.tmp_insert = os.path.join(tmp_folder, unique_id + "_insert.txt")
        with open(self.tmp_insert, "w") as writer:
            writer.write("""## htsjdk.samtools.metrics.StringHeader
# picard.analysis.CollectInsertSizeMetrics INPUT=/path/to/bam REFERENCE_SEQUENCE=/path/to/reference ASSUME_SORTED=true OUTPUT=/path/to/output H=/path/to/histogram STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
## htsjdk.samtools.metrics.StringHeader
# Started on: Fri Nov 13 20:25:43 CET 2020

## METRICS CLASS	picard.analysis.InsertSizeMetrics
MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
286	166	70	14	244976752	260.992864	87.727108	14114058	FR	23	43	67	99	141	193	227	243	269	299	437
98	50	1	96	105	1.1	1.1	1	RF	1	1	1	1	1	1	1	1	1	1	1
98	40	1	97	98	1.1	1.1	1	TANDEM	1	1	1	1	1	1	1	1	1	1	1

## HISTOGRAM	java.lang.Integer
insert_size	All_Reads.fr_count	All_Reads.rf_count	All_Reads.tandem_count
1	0	100	0
2	0	100	0
3	0	100	0
4	99	100	0
5	0	100	0
6	0	100	0
7	0	100	0
8	0	100	0
9	0	100	0
10	0	100	0
11	0	100	0
12	0	100	0
13	0	100	50
14	0	100	0
15	0	100	0
16	0	100	0
17	0	100	0
18	0	100	0
19	0	10	0
20	0	100	0""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_variants_call, self.tmp_wgs, self.tmp_insert]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testCollectInsertSizeMetrics(self):
        observed = PicardReader(self.tmp_insert)
        expected = {
            "command": "picard.analysis.CollectInsertSizeMetrics INPUT=/path/to/bam REFERENCE_SEQUENCE=/path/to/reference ASSUME_SORTED=true OUTPUT=/path/to/output H=/path/to/histogram STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false",
            "metrics": [
                {
                    'MEDIAN_INSERT_SIZE': 286,
                    'MODE_INSERT_SIZE': 166,
                    'MEDIAN_ABSOLUTE_DEVIATION': 70,
                    'MIN_INSERT_SIZE': 14,
                    'MAX_INSERT_SIZE': 244976752,
                    'MEAN_INSERT_SIZE': 260.992864,
                    'STANDARD_DEVIATION': 87.727108,
                    'READ_PAIRS': 14114058,
                    'PAIR_ORIENTATION': 'FR',
                    'WIDTH_OF_10_PERCENT': 23,
                    'WIDTH_OF_20_PERCENT': 43,
                    'WIDTH_OF_30_PERCENT': 67,
                    'WIDTH_OF_40_PERCENT': 99,
                    'WIDTH_OF_50_PERCENT': 141,
                    'WIDTH_OF_60_PERCENT': 193,
                    'WIDTH_OF_70_PERCENT': 227,
                    'WIDTH_OF_80_PERCENT': 243,
                    'WIDTH_OF_90_PERCENT': 269,
                    'WIDTH_OF_95_PERCENT': 299,
                    'WIDTH_OF_99_PERCENT': 437,
                    'SAMPLE': None,
                    'LIBRARY': None,
                    'READ_GROUP': None
                }, {
                    'MEDIAN_INSERT_SIZE': 98, 'MODE_INSERT_SIZE': 50, 'MEDIAN_ABSOLUTE_DEVIATION': 1, 'MIN_INSERT_SIZE': 96, 'MAX_INSERT_SIZE': 105, 'MEAN_INSERT_SIZE': 1.1, 'STANDARD_DEVIATION': 1.1, 'READ_PAIRS': 1, 'PAIR_ORIENTATION': 'RF', 'WIDTH_OF_10_PERCENT': 1, 'WIDTH_OF_20_PERCENT': 1, 'WIDTH_OF_30_PERCENT': 1, 'WIDTH_OF_40_PERCENT': 1, 'WIDTH_OF_50_PERCENT': 1, 'WIDTH_OF_60_PERCENT': 1, 'WIDTH_OF_70_PERCENT': 1, 'WIDTH_OF_80_PERCENT': 1, 'WIDTH_OF_90_PERCENT': 1, 'WIDTH_OF_95_PERCENT': 1, 'WIDTH_OF_99_PERCENT': 1, 'SAMPLE': None, 'LIBRARY': None, 'READ_GROUP': None
                }, {
                    'MEDIAN_INSERT_SIZE': 98, 'MODE_INSERT_SIZE': 40, 'MEDIAN_ABSOLUTE_DEVIATION': 1, 'MIN_INSERT_SIZE': 97, 'MAX_INSERT_SIZE': 98, 'MEAN_INSERT_SIZE': 1.1, 'STANDARD_DEVIATION': 1.1, 'READ_PAIRS': 1, 'PAIR_ORIENTATION': 'TANDEM', 'WIDTH_OF_10_PERCENT': 1, 'WIDTH_OF_20_PERCENT': 1, 'WIDTH_OF_30_PERCENT': 1, 'WIDTH_OF_40_PERCENT': 1, 'WIDTH_OF_50_PERCENT': 1, 'WIDTH_OF_60_PERCENT': 1, 'WIDTH_OF_70_PERCENT': 1, 'WIDTH_OF_80_PERCENT': 1, 'WIDTH_OF_90_PERCENT': 1, 'WIDTH_OF_95_PERCENT': 1, 'WIDTH_OF_99_PERCENT': 1, 'SAMPLE': None, 'LIBRARY': None, 'READ_GROUP': None
                }
            ],
            "histogram": {
                'insert_size': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                'All_Reads.fr_count': [0, 0, 0, 99, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'All_Reads.rf_count': [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 10, 100],
                'All_Reads.tandem_count': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0]
            }
        }
        self.assertEqual(expected["command"], observed.command)
        self.assertEqual(expected["metrics"], observed.metrics)
        self.assertEqual(expected["histogram"], observed.histogram)

    def testCollectArraysVariantCallingMetrics(self):
        observed = PicardReader(self.tmp_variants_call)
        expected = {
            "command": "CollectArraysVariantCallingMetrics  --INPUT /cromwell_root/broad-gotc-dev-cromwell-execution/Arrays/9972b43c-df1e-414a-9017-24edc37b5250/call-IlluminaGenotypingArray/IlluminaGenotypingArray/805b2e59-add2-4794-ae89-350d69a0f7c3/call-MergePedIntoVcf/7991775143_R01C01.vcf.gz --OUTPUT 7991775143_R01C01 --CALL_RATE_PF_THRESHOLD 0.98 --DBSNP /cromwell_root/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz  --NUM_PROCESSORS 0 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --GA4GH_CLIENT_SECRETS client_secrets.json --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false",
            "metrics": [
                {
                    'CHIP_WELL_BARCODE': '7991775143_R01C01',
                    'SAMPLE_ALIAS': 'NA12878',
                    'ANALYSIS_VERSION': 1,
                    'CHIP_TYPE': 'HumanExome-12v1-1_A',
                    'AUTOCALL_PF': 'Y',
                    'AUTOCALL_DATE': '2019-12-14T03:05:00+0000',
                    'IMAGING_DATE': '2012-12-16T16:10:10+0000',
                    'IS_ZCALLED': 'Y',
                    'AUTOCALL_GENDER': 'F',
                    'FP_GENDER': 'U',
                    'REPORTED_GENDER': 'F',
                    'GENDER_CONCORDANCE_PF': 'Y',
                    'HET_PCT': 0.042993,
                    'CLUSTER_FILE_NAME': 'HumanExomev1_1_CEPH_A.egt',
                    'P95_GREEN': 4251,
                    'P95_RED': 5479,
                    'AUTOCALL_VERSION': '3.0.0',
                    'ZCALL_VERSION': '1.0.0.0',
                    'EXTENDED_MANIFEST_VERSION': 1.3,
                    'HET_HOMVAR_RATIO': 1.525553,
                    'SCANNER_NAME': 'N370',
                    'NUM_ASSAYS': 242513,
                    'NUM_NON_FILTERED_ASSAYS': 242493,
                    'NUM_FILTERED_ASSAYS': 20,
                    'NUM_ZEROED_OUT_ASSAYS': 0,
                    'NUM_SNPS': 242450,
                    'NUM_INDELS': 43, 'NUM_CALLS': 242318,
                    'NUM_AUTOCALL_CALLS': 241250,
                    'NUM_NO_CALLS': 175,
                    'NUM_IN_DB_SNP': 238085,
                    'NOVEL_SNPS': 4365,
                    'PCT_DBSNP': 0.981996,
                    'CALL_RATE': 0.999278,
                    'AUTOCALL_CALL_RATE': 0.994874,
                    'NUM_SINGLETONS': 10418
                }
            ],
            "histogram": None
        }
        self.assertEqual(expected["command"], observed.command)
        self.assertEqual(expected["metrics"], observed.metrics)
        self.assertEqual(expected["histogram"], observed.histogram)

    def testCollectWgsMetrics(self):
        observed = PicardReader(self.tmp_wgs)
        expected = {
            "command": "CollectWgsMetrics INPUT=/cromwell_root/broad-gotc-dev-cromwell-execution/WholeGenomeGermlineSingleSample/e36be5d3-ddbf-4ee6-a98d-2a1933826550/call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/da7029e5-cf73-4d9b-8369-23c3aea72f19/call-GatherBamFiles/NA12878_PLUMBING.bam OUTPUT=NA12878_PLUMBING.wgs_metrics INCLUDE_BQ_HISTOGRAM=true INTERVALS=/cromwell_root/gcp-public-data--broad-references/hg38/v0/wgs_coverage_regions.hg38.interval_list USE_FAST_ALGORITHM=true READ_LENGTH=250 VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=/cromwell_root/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta    MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 COVERAGE_CAP=250 LOCUS_ACCUMULATION_CAP=100000 STOP_AFTER=-1 COUNT_UNPAIRED=false SAMPLE_SIZE=10000 ALLELE_FRACTION=[0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5] VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false",
            "metrics": [{'GENOME_TERRITORY': 2745186692, 'MEAN_COVERAGE': 0.003328, 'SD_COVERAGE': 0.059988, 'MEDIAN_COVERAGE': 0, 'MAD_COVERAGE': 0, 'PCT_EXC_ADAPTER': 4.5e-05, 'PCT_EXC_MAPQ': 0.044363, 'PCT_EXC_DUPE': 0.000482, 'PCT_EXC_UNPAIRED': 0.00674, 'PCT_EXC_BASEQ': 0.121016, 'PCT_EXC_OVERLAP': 0.090623, 'PCT_EXC_CAPPED': 0, 'PCT_EXC_TOTAL': 0.263269, 'PCT_1X': 0.003305, 'PCT_5X': 1e-06, 'PCT_10X': 0, 'PCT_15X': 0, 'PCT_20X': 0, 'PCT_25X': 0, 'PCT_30X': 0, 'PCT_40X': 0, 'PCT_50X': 0, 'PCT_60X': 0, 'PCT_70X': 0, 'PCT_80X': 0, 'PCT_90X': 0, 'PCT_100X': 0, 'HET_SNP_SENSITIVITY': 0.000105, 'HET_SNP_Q': 0}],
            "histogram": {
                'coverage': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
                'high_quality_coverage_count': [2736113146, 9041943, 24393, 3102, 1433, 779, 861, 236, 149, 83, 75, 54, 0, 0, 0, 0, 0, 1, 1, 10, 0, 2, 0, 3, 11, 7, 16, 1, 13, 17, 23, 12, 22],
                'unfiltered_baseq_count': [0, 0, 0, 730, 940, 10184, 1526, 0, 0, 0, 45052, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10192669, 0, 0]
            }
        }
        self.assertEqual(expected["command"], observed.command)
        self.assertEqual(expected["metrics"], observed.metrics)
        self.assertEqual(expected["histogram"], observed.histogram)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
