#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'


import os
import sys
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.elementbio.samplesheet import SampleSheet


class TestSampleSheet(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        self.tmp_file = os.path.join(tmp_folder, "RunManifest.csv")
        self.test_cases = [
            # {
            #     "name": "no_sample",
            #     "content": """[SAMPLES]"""
            # },
            {
                "name": "no_barcode",
                "content": """[SAMPLES]
SampleName
Sample-1"""
            }, {
                "name": "one_barcode",
                "content": """[SAMPLES]
SampleName,Index1
Sample-1,ACGTGTAGC
Sample-2,ATTGCCATA"""
            }, {
                "name": "only_sample_section_clean_comma",
                "content": """[SAMPLES]
SampleName,Index1,Index2,Lane
PhiX,ACGTGTAGC,GCTAGTGCA,1+2
PhiX,CACATGCTG,AGACACTGT,1+2

Sample-1,TGTAGGCCA,TCTAGCCTC,1
Sample-2,ATTGCCATA,ATATTGGTG,2
Sample-3,ATTGCCATA,ATATTGGTG,1+2"""
            }, {
                "name": "all_sections_comma",
                "content": """[SETTINGS],,,
SettingName,Value,Lane,
# Lines can be commented out with a leading # sign. ,,,
# Replace the example adapters below with the adapter used in the kit. ,,,
R1AdapterTrim,FALSE,1+2,
R2AdapterTrim,FALSE,2,
,,,
[RUNVALUES],,,
Keyname,Value,,
key,my_val,,
,,,
[SAMPLES],,,
SampleName,Index1,Index2,Lane
# Make sure that in the Index1 and Index2 columns the index sequences for each library is in the same orientation (samples and PhiX controls).,,,
PhiX,ACGTGTAGC,GCTAGTGCA,1+2
PhiX,CACATGCTG,AGACACTGT,1+2
,,,
# Fill in the correct sample schema associated with the Elevate Cloudbreak Preparation Workflow for all sequenced samples.,,,
Sample-1,TGTAGGCCA,TCTAGCCTC,1
Sample-2,ATTGCCATA,ATATTGGTG,2
Sample-3,ATTGCCATA,ATATTGGTG,1+2"""
            }, {
                "name": "sections_settings_samples_no_lane",
                "content": """[SETTINGS]
SettingName,Value
R1AdapterTrim,FALSE
R2AdapterTrim,FALSE

[SAMPLES]
SampleName,Index1,Index2,Design,Description
Sample-1,TGTAGGCCA,TCTAGCCTC,Enrich_Lymph_V1,
Sample-2,ATTGCCATA,ATATTGGTG,Enrich_Solid_V2,Melanoma"""
            }
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testIsValid(self):
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            self.assertTrue(
                SampleSheet.isValid(self.tmp_file)
            )
        test_cases_invalid = [
            {
                "name": "empty",
                "content": """"""
            },
            {
                "name": "missing_SampleName",
                "content": """[SAMPLES]
I1,I2
ATGC,ATTG
"""
            },
            {
                "name": "blank_SampleName",
                "content": """[SAMPLES]
SampleName,I1
,ATTG
"""
            },
            {
                "name": "unaivalable _section",
                "content": """[ODD_SECTION]
[SAMPLES]
SampleName,I1
splA,ATTG
"""
            },
            {
                "name": "inconsistent_number_fields",
                "content": """[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma
splB,GGTC,Enrich_V1,Lymphome,test
"""
            },
            {
                "name": "invalid_settings_titles",
                "content": """[SETTINGS]
key,val

[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma"""
            },
            {
                "name": "missing_run_values_titles",
                "content": """[RUNVALUES],,,
operator,John

[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma"""
            },
            {
                "name": "invalid_number_fields_run_values",
                "content": """[RUNVALUES]
KeyName,Value
operator,John,Smith

[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma"""
            },
            {
                "name": "invalid_number_fields_settings",
                "content": """[SETTINGS]
SettingName,Value
R1AdapterTrim,FALSE,1,8

[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma"""
            },
            {
                "name": "inconsistent_number_fields_settings",
                "content": """[SETTINGS]
SettingName,Value
R1AdapterTrim,FALSE,1

[SAMPLES]
SampleName,I1,Design,Disease
splA,ATTG,Enrich_V1,Melanoma"""
            },
        ]
        for curr_test in test_cases_invalid:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            self.assertFalse(
                SampleSheet.isValid(self.tmp_file)
            )

    def testParse(self):
        expected = {
            # "no_sample": {
            #     "header": {},
            #     "settings": {},
            #     "samples": [
            #         {
            #             'barcodes': dict(),
            #             'basename': 'Sample-1', ###################
            #             'description': None,
            #             'id': 'Sample-1',
            #             'metadata': {}
            #         }
            #     ]
            # },
            "no_barcode": {
                "header": {},
                "settings": {},
                "samples": [
                    {
                        'barcodes': dict(),
                        'basename': 'Sample-1',
                        'description': None,
                        'id': 'Sample-1',
                        'metadata': {}
                    }
                ]
            },
            "one_barcode": {
                "header": {},
                "settings": {},
                "samples": [
                    {
                        'barcodes': {"index": "ACGTGTAGC"},
                        'basename': 'Sample-1',
                        'description': None,
                        'id': 'Sample-1',
                        'metadata': {}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA"},
                        'basename': 'Sample-2',
                        'description': None,
                        'id': 'Sample-2',
                        'metadata': {}
                    },
                ]
            },
            "only_sample_section_clean_comma": {
                "header": {},
                "settings": {},
                "samples": [
                    {
                        'barcodes': {"index": "ACGTGTAGC", "index2": "GCTAGTGCA"},
                        'basename': 'PhiX',
                        'description': None,
                        'id': 'PhiX',
                        'metadata': {"Lane": ["1", "2"]}
                    },
                    {
                        'barcodes': {"index": "CACATGCTG", "index2": "AGACACTGT"},
                        'basename': 'PhiX',
                        'description': None,
                        'id': 'PhiX',
                        'metadata': {"Lane": ["1", "2"]}
                    },
                    {
                        'barcodes': {"index": "TGTAGGCCA", "index2": "TCTAGCCTC"},
                        'basename': 'Sample-1',
                        'description': None,
                        'id': 'Sample-1',
                        'metadata': {"Lane": ["1"]}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA", "index2": "ATATTGGTG"},
                        'basename': 'Sample-2',
                        'description': None,
                        'id': 'Sample-2',
                        'metadata': {"Lane": ["2"]}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA", "index2": "ATATTGGTG"},
                        'basename': 'Sample-3',
                        'description': None,
                        'id': 'Sample-3',
                        'metadata': {"Lane": ["1", "2"]}
                    }
                ]
            },
            "all_sections_comma": {
                "header": {"key": "my_val"},
                "settings": {
                    "R1AdapterTrim": {"1": "FALSE", "2": "FALSE"},
                    "R2AdapterTrim": {"2": "FALSE"}
                },
                "samples": [
                    {
                        'barcodes': {"index": "ACGTGTAGC", "index2": "GCTAGTGCA"},
                        'basename': 'PhiX',
                        'description': None,
                        'id': 'PhiX',
                        'metadata': {"Lane": ["1", "2"]}
                    },
                    {
                        'barcodes': {"index": "CACATGCTG", "index2": "AGACACTGT"},
                        'basename': 'PhiX',
                        'description': None,
                        'id': 'PhiX',
                        'metadata': {"Lane": ["1", "2"]}
                    },
                    {
                        'barcodes': {"index": "TGTAGGCCA", "index2": "TCTAGCCTC"},
                        'basename': 'Sample-1',
                        'description': None,
                        'id': 'Sample-1',
                        'metadata': {"Lane": ["1"]}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA", "index2": "ATATTGGTG"},
                        'basename': 'Sample-2',
                        'description': None,
                        'id': 'Sample-2',
                        'metadata': {"Lane": ["2"]}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA", "index2": "ATATTGGTG"},
                        'basename': 'Sample-3',
                        'description': None,
                        'id': 'Sample-3',
                        'metadata': {"Lane": ["1", "2"]}
                    }
                ]
            },
            "sections_settings_samples_no_lane": {
                "header": {},
                "settings": {
                    "R1AdapterTrim": {"1": "FALSE", "2": "FALSE"},
                    "R2AdapterTrim": {"1": "FALSE", "2": "FALSE"}
                },
                "samples": [
                    {
                        'barcodes': {"index": "TGTAGGCCA", "index2": "TCTAGCCTC"},
                        'basename': 'Sample-1',
                        'description': None,
                        'id': 'Sample-1',
                        'metadata': {"Design": "Enrich_Lymph_V1"}
                    },
                    {
                        'barcodes': {"index": "ATTGCCATA", "index2": "ATATTGGTG"},
                        'basename': 'Sample-2',
                        'description': "Melanoma",
                        'id': 'Sample-2',
                        'metadata': {"Design": "Enrich_Solid_V2"}
                    }
                ]
            }
        }
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            samplesheet = SampleSheet(self.tmp_file)
            observed = {
                "header": samplesheet.header,
                "samples": [spl.toDict() for spl in samplesheet.samples],
                "settings": samplesheet.settings
            }
            self.assertEqual(observed, expected[curr_test["name"]])


if __name__ == "__main__":
    unittest.main()
