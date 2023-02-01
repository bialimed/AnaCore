#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina.samplesheet import SampleSheetFactory, SampleSheetV1, SampleSheetV2


class TestSampleSheet(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_SampleSheet.csv")
        self.test_cases = [
            {
                "name": "TSO500_ADN-ARN_V2",
                "version": 2,
                "content": """[Header],,,,,,,,,,,,
FileFormatVersion,2,,,,,,,,,,,
RunName,Test_V2_SampleSheet,,,,,,,,,,,
InstrumentType,NovaSeq6000,,,,,,,,,,,
InstrumentPlatform,...,,,,,,,,,,,
[Reads],,,,,,,,,,,,
Read1Cycles,101,,,,,,,,,,,
Read2Cycles,101,,,,,,,,,,,
Index1Cycles,10,,,,,,,,,,,
Index2Cycles,10,,,,,,,,,,,
[TSO500S_Settings],,,,,,,,,,,,
SoftwareVersion,1.0.0,,,,,,,,,,,
AdapterRead1,CTGTCTCTTATACACATCTCCGAGCCCACGAGAC,,,,,,,,,,,
AdapterRead2,CTGTCTCTTATACACATCTGACGCTGCCGACGA,,,,,,,,,,,
AdapterBehavior,trim,,,,,,,,,,,
MinimumTrimmedReadLength,35,,,,,,,,,,,
MaskShortReads,35,,,,,,,,,,,
OverrideCycles,U7N1Y93;I10;I10;U7N1Y93,,,,,,,,,,,
[TSO500S_Data],,,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_ID,index,index2,I7_Index_ID,I5_Index_ID,Project,Description,Pair_ID,Sample_Type
Sample1-DNA,Sample1-DNA,,A12,UDP0089,GTCCGTAAGC,CGTGTATCTT,UDP0089,UDP0089,,,Sample1,DNA
Sample1-RNA,Sample1-RNA,,A10,UDP0073,CCTGCGGAAC,ATCATAGGCT,UDP0073,UDP0073,,,Sample1,RNA
[BCLConvert_Settings],,,,,,,,,,,,
SoftwareVersion,3.6.3,,,,,,,,,,,
AdapterRead1,CTGTCTCTTATACACATCTCCGAGCCCACGAGAC,,,,,,,,,,,
AdapterRead2,CTGTCTCTTATACACATCTGACGCTGCCGACGA,,,,,,,,,,,
MinimumTrimmedReadLength,35,,,,,,,,,,,
MaskShortReads,35,,,,,,,,,,,
OverrideCycles,U7N1Y93;I10;I10;U7N1Y93,,,,,,,,,,,
[BCLConvert_Data],,,,,,,,,,,,
Sample_ID,index,index2,,,,,,,,,,
Sample1-DNA,GTCCGTAAGC,CGTGTATCTT,,,,,,,,,,
Sample1-RNA,CCTGCGGAAC,ATCATAGGCT,,,,,,,,,,"""
            },
            {
                "name": "TSO500_ADN-ARN_V1",
                "version": 1,
                "content": """[Header],,,,,,,,
IEMFileVersion,4,,,,,,,
Investigator Name,User Name,,,,,,,
Experiment Name,Experiment,,,,,,,
Date,10/3/2018,,,,,,,
Workflow,GenerateFASTQ,,,,,,,
Application,NextSeq FASTQ Only,,,,,,,
Assay,,,,,,,,
Description,,,,,,,,
Chemistry,Default,,,,,,,
,,,,,,,,
[Reads],,,,,,,,
101,,,,,,,,
101,,,,,,,,
,,,,,,,,
[Settings],,,,,,,,
AdapterRead1,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,
AdapterBehavior,trim,,,,,,,
MinimumTrimmedReadLength,35,,,,,,,
MaskShortReads,35,,,,,,,
OverrideCycles,U7N1Y93;I8;I8;U7N1Y93,,,,,,,
,,,,,,,,
[Data],,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_ID,index,index2,Sample_Type,Pair_ID
Test_Sample_UP01,,,,UP01,TCCGGAGA,AGGATAGG,DNA,Test_Sample_UP01
Test_Sample_UP02,,,,UP02,CTGAAGCT,TCAGAGCC,DNA,Test_Sample_UP02"""
            },
            {
                "name": "simple_V2_SE",
                "version": 2,
                "content": """[Header],,,
FileFormatVersion,2,,
RunName,MyRun,,
InstrumentPlatform,NextSeq1k2k,,
InstrumentType,NextSeq2000,,
,,,
[Reads],,,
Read1Cycles,125,,
Index1Cycles,8,,
Index2Cycles,8,,
,,,
[BCLConvert_Settings],,,
SoftwareVersion,x.y.z,,
,,,
[BCLConvert_Data],,,
Lane,Sample_ID,index,index2
1,S01-TOO-12plex-P1-rep1,ATCCACTG,AGGTGCGT
1,S02-TOO-12plex-P1-rep2,GCTTGTCA,GAACATAC
,,,
,,,
,,,"""
            },
            {
                "name": "nextera_V1_PE_2-INDEX",
                "version": 1,
                "content": """[Header]
IEMFileVersion,4
InvestigatorName,user1
ProjectName,Sequences
Date,01/03/2021
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,Nextera
Description,Sequencing
Chemistry,Amplicon

[Reads]
151
100

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,I7_Index_ID,index,I5_Index_ID,index2
1,DG-1-A01,i7A01,GTACGTCA,i5A01,GAGGTAGT
2,DG-1-A02,i7A02,TGCAGTTA,i5A01,GAGGTAGT

[Clinical]
laboratory,Broad
trial_id,178

"""
            },
            {
                "name": "SureSelectXT_V1_SE_1-INDEX",
                "version": 1,
                "content": """[Header],,,,,,,,
IEM1FileVersion,4,,,,,,,
Investigator Name,jdoe,,,,,,,
Experiment Name,exp001,,,,,,,
Date,11/16/2017,,,,,,,
Workflow,SureSelectXT,,,,,,,
Application,NextSeq FASTQ Only,,,,,,,
Assay,SureSelectXT,,,,,,,
Description,A description of this flow cell,,,,,,,
Chemistry,Default,,,,,,,
,,,,,,,,
[Reads],,,,,,,,
151,,,,,,,,
,,,,,,,,
[Settings],,,,,,,,
CreateFastqForIndexReads,1,,,,,,,
BarcodeMismatches,2,,,,,,,
,,,,,,,,
[Data],,,,,,,,
Sample_ID,Sample_Name,index,Description,Library_ID,Read_Structure,Reference_Name,Sample_Project,Target_Set
1823A,1823A-tissue,GAATCTGA,0.5x treatment,2017-01-20,151T8B151T,mm10,exp001,Intervals-001
1823B,1823B-tissue,AGCAGGAA,0.5x treatment,2017-01-20,151T8B151T,mm10,exp001,Intervals-001
,,,,,,,,
,,,,,,,,
,,,,,,,,"""
            }
        ]

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testIsValid(self):
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            self.assertEqual(
                SampleSheetV1.isValid(self.tmp_file),
                curr_test["version"] == 1
            )
            self.assertEqual(
                SampleSheetV2.isValid(self.tmp_file),
                curr_test["version"] == 2
            )
            if curr_test["version"] == 2:
                with open(self.tmp_file, "w") as writer:
                    writer.write(curr_test["content"].replace("FileFormatVersion,2", "FileFormatVersion,1"))
                self.assertFalse(
                    SampleSheetV1.isValid(self.tmp_file)
                )
                self.assertFalse(
                    SampleSheetV2.isValid(self.tmp_file)
                )
            if curr_test["version"] == 1:
                with open(self.tmp_file, "w") as writer:
                    writer.write(curr_test["content"].replace("[Data]", "[BCLConvert_Data]"))
                self.assertFalse(
                    SampleSheetV1.isValid(self.tmp_file)
                )
                self.assertFalse(
                    SampleSheetV2.isValid(self.tmp_file)
                )

    def testParse(self):
        expected = {
            "TSO500_ADN-ARN_V2": {
                'extra': {
                    'TSO500S': {
                        'settings': {
                            'SoftwareVersion': '1.0.0',
                            'AdapterRead1': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
                            'AdapterRead2': 'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
                            'AdapterBehavior': 'trim',
                            'MinimumTrimmedReadLength': '35',
                            'MaskShortReads': '35',
                            'OverrideCycles': 'U7N1Y93;I10;I10;U7N1Y93'
                        },
                        'data': [
                            {
                                'Sample_ID': 'Sample1-DNA', 'Sample_Name': 'Sample1-DNA',
                                'Sample_Plate': '', 'Sample_Well': 'A12',
                                'Index_ID': 'UDP0089', 'index': 'GTCCGTAAGC',
                                'index2': 'CGTGTATCTT', 'I7_Index_ID': 'UDP0089',
                                'I5_Index_ID': 'UDP0089', 'Project': '',
                                'Description': '', 'Pair_ID': 'Sample1',
                                'Sample_Type': 'DNA'
                            },
                            {
                                'Sample_ID': 'Sample1-RNA', 'Sample_Name': 'Sample1-RNA',
                                'Sample_Plate': '', 'Sample_Well': 'A10',
                                'Index_ID': 'UDP0073', 'index': 'CCTGCGGAAC',
                                'index2': 'ATCATAGGCT', 'I7_Index_ID': 'UDP0073',
                                'I5_Index_ID': 'UDP0073', 'Project': '',
                                'Description': '', 'Pair_ID': 'Sample1',
                                'Sample_Type': 'RNA'
                            }
                        ]
                    }
                },
                'header': {
                    'FileFormatVersion': '2', 'RunName': 'Test_V2_SampleSheet',
                    'InstrumentType': 'NovaSeq6000', 'InstrumentPlatform': '...',
                    'Description': ''
                },
                'manifests': {},
                'reads': {
                    'phases': [
                        {'is_index': False, 'nb_cycles': 101},
                        {'is_index': True, 'nb_cycles': 10},
                        {'is_index': True, 'nb_cycles': 10},
                        {'is_index': False, 'nb_cycles': 101}
                    ]
                },
                'samples': [
                    {
                        'Sample_ID': 'Sample1-DNA', 'index': 'GTCCGTAAGC',
                        'index2': 'CGTGTATCTT', 'Sample_Basename': 'Sample1-DNA',
                        'Library_Basename': 'Sample1-DNA_S1'
                    },
                    {
                        'Sample_ID': 'Sample1-RNA', 'index': 'CCTGCGGAAC',
                        'index2': 'ATCATAGGCT', 'Sample_Basename': 'Sample1-RNA',
                        'Library_Basename': 'Sample1-RNA_S2'
                    }
                ]
            },
            "TSO500_ADN-ARN_V1": {
                'extra': {},
                'header': {
                    'IEMFileVersion': '4', 'Investigator Name': 'User Name',
                    'Experiment Name': 'Experiment', 'Date': '10/3/2018',
                    'Workflow': 'GenerateFASTQ', 'Application': 'NextSeq FASTQ Only',
                    'Assay': '', 'Description': '', 'Chemistry': 'Default'
                },
                'manifests': {},
                'reads': {'nb_cycles': {'R1': 101, 'R2': 101}},
                'samples': [
                    {
                        'Sample_ID': 'Test_Sample_UP01', 'Sample_Name': '',
                        'Sample_Plate': '', 'Sample_Well': '', 'Index_ID': 'UP01',
                        'index': 'TCCGGAGA', 'index2': 'AGGATAGG',
                        'Sample_Type': 'DNA', 'Pair_ID': 'Test_Sample_UP01',
                        'Sample_Basename': 'Test-Sample-UP01',
                        'Library_Basename': 'Test-Sample-UP01_S1'
                    },
                    {
                        'Sample_ID': 'Test_Sample_UP02', 'Sample_Name': '',
                        'Sample_Plate': '', 'Sample_Well': '', 'Index_ID': 'UP02',
                        'index': 'CTGAAGCT', 'index2': 'TCAGAGCC',
                        'Sample_Type': 'DNA', 'Pair_ID': 'Test_Sample_UP02',
                        'Sample_Basename': 'Test-Sample-UP02',
                        'Library_Basename': 'Test-Sample-UP02_S2'
                    }
                ]
            },
            "simple_V2_SE": {
                "extra": {},
                "header": {
                    "FileFormatVersion": "2",
                    "RunName": "MyRun",
                    "InstrumentPlatform": "NextSeq1k2k",
                    "InstrumentType": "NextSeq2000",
                    "Description": ""
                },
                "manifests": {},
                "reads": {
                    "phases": [
                        {'is_index': False, 'nb_cycles': 125},
                        {'is_index': True, 'nb_cycles': 8},
                        {'is_index': True, 'nb_cycles': 8},
                    ],
                },
                "samples":  [
                    {
                        "Lane": "1", "Sample_ID": "S01-TOO-12plex-P1-rep1",
                        "index": "ATCCACTG", "index2": "AGGTGCGT",
                        "Sample_Basename": "S01-TOO-12plex-P1-rep1",
                        "Library_Basename": "S01-TOO-12plex-P1-rep1_S1"
                    },
                    {
                        "Lane": "1", "Sample_ID": "S02-TOO-12plex-P1-rep2",
                        "index": "GCTTGTCA", "index2": "GAACATAC",
                        "Sample_Basename": "S02-TOO-12plex-P1-rep2",
                        "Library_Basename": "S02-TOO-12plex-P1-rep2_S2"
                    }
                ]
            },
            "nextera_V1_PE_2-INDEX": {
                'extra': {
                    "Clinical": {
                        "laboratory": "Broad",
                        "trial_id": "178"
                    }
                },
                'header': {
                    'IEMFileVersion': '4', 'InvestigatorName': 'user1',
                    'ProjectName': 'Sequences', 'Date': '01/03/2021',
                    'Workflow': 'GenerateFASTQ', 'Application': 'FASTQ Only',
                    'Assay': 'Nextera', 'Description': 'Sequencing',
                    'Chemistry': 'Amplicon'
                },
                'manifests': {},
                'reads': {'nb_cycles': {'R1': 151, 'R2': 100}},
                'samples': [
                    {
                        'Sample_ID': '1', 'Sample_Name': 'DG-1-A01',
                        'I7_Index_ID': 'i7A01', 'index': 'GTACGTCA',
                        'I5_Index_ID': 'i5A01', 'index2': 'GAGGTAGT',
                        'Sample_Basename': '1', 'Library_Basename': '1_S1'
                    },
                    {
                        'Sample_ID': '2', 'Sample_Name': 'DG-1-A02',
                        'I7_Index_ID': 'i7A02', 'index': 'TGCAGTTA',
                        'I5_Index_ID': 'i5A01', 'index2': 'GAGGTAGT',
                        'Sample_Basename': '2', 'Library_Basename': '2_S2'
                    }
                ]
            },
            "SureSelectXT_V1_SE_1-INDEX": {
                'extra': {},
                'header': {
                    'IEM1FileVersion': '4', 'Investigator Name': 'jdoe',
                    'Experiment Name': 'exp001', 'Date': '11/16/2017',
                    'Workflow': 'SureSelectXT', 'Application': 'NextSeq FASTQ Only',
                    'Assay': 'SureSelectXT', 'Description': 'A description of this flow cell',
                    'Chemistry': 'Default'
                },
                'manifests': {},
                'reads': {'nb_cycles': {'R1': 151}},
                'samples': [
                    {
                        'Sample_ID': '1823A', 'Sample_Name': '1823A-tissue',
                        'index': 'GAATCTGA', 'Description': '0.5x treatment',
                        'Library_ID': '2017-01-20', 'Read_Structure': '151T8B151T',
                        'Reference_Name': 'mm10', 'Sample_Project': 'exp001',
                        'Target_Set': 'Intervals-001', 'Sample_Basename': '1823A',
                        'Library_Basename': '1823A_S1'
                    },
                    {
                        'Sample_ID': '1823B', 'Sample_Name': '1823B-tissue',
                        'index': 'AGCAGGAA', 'Description': '0.5x treatment',
                        'Library_ID': '2017-01-20', 'Read_Structure': '151T8B151T',
                        'Reference_Name': 'mm10', 'Sample_Project': 'exp001',
                        'Target_Set': 'Intervals-001', 'Sample_Basename': '1823B',
                        'Library_Basename': '1823B_S2'
                    }
                ]
            }
        }
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as writer:
                writer.write(curr_test["content"])
            observed = SampleSheetFactory.get(self.tmp_file)
            observed = {k: v for k, v in observed.__dict__.items() if k != "filepath"}
            self.assertEqual(observed, expected[curr_test["name"]])


if __name__ == "__main__":
    unittest.main()
