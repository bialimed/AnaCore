#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest
from datetime import datetime

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina import Bcl2fastqLog, getInfFromSeqID, RTAComplete, RunInfo, RunParameters


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestBcl2fastqLog(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_bcl2fastq_log.txt")
        self.test_cases = [
            {  # Complete
                "content": """BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.

2021-06-03 17:18:10 [2afc55848400] Command-line invocation: bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2
2021-06-03 17:18:10 [2afc55848400] INFO: Minimum log level: INFO
2021-06-03 23:22:13 [2afc55848400] Processing completed with 0 errors and 1 warnings.""",
                "expected": {
                    "command": "bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2",
                    "complete": {
                        "nb_errors": 0,
                        "nb_warnings": 1
                    },
                    "version": "v2.20.0.422",
                    "start_time": datetime(2021, 6, 3, 17, 18, 10),
                    "end_time": datetime(2021, 6, 3, 23, 22, 13)
                }
            }, {  # Incomplete
                "content": """BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.

2021-06-03 17:18:10 [2afc55848400] Command-line invocation: bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2
2021-06-03 17:18:10 [2afc55848400] INFO: Minimum log level: INFO""",
                "expected": {
                    "command": "bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2",
                    "complete": None,
                    "version": "v2.20.0.422",
                    "start_time": datetime(2021, 6, 3, 17, 18, 10),
                    "end_time": None
                }
            }
        ]

    def testParse(self):
        expected = [curr_test["expected"] for curr_test in self.test_cases]
        observed = []
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as handle:
                handle.write(curr_test["content"] + "\n")
            res = Bcl2fastqLog(self.tmp_file)
            observed.append({
                "command": res.command,
                "complete": res.complete,
                "version": res.version,
                "start_time": res.start_time,
                "end_time": res.end_time
            })
        self.assertEqual(expected, observed)


class TestRTAComplete(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_RTAComplete.txt")
        self.test_cases = [
            {
                "content": "11/1/2017,15:11:43.174,Illumina RTA 1.18.54",
                "expected": {"RTA_version": "1.18.54", "end_date": datetime(2017, 11, 1, 15, 11, 43)}  # datetime.fromtimestamp(1484147503)
            }, {
                "content": "RTA 2.4.11 completed on 11/14/2019 4:56:45 AM",
                "expected": {"RTA_version": "2.4.11", "end_date": datetime(2019, 11, 14, 4, 56, 45)}  # datetime.fromtimestamp(1573707405)
            }, {
                "content": "RTA 2.4.11 completed on 11/14/2019 4:56:45 PM",
                "expected": {"RTA_version": "2.4.11", "end_date": datetime(2019, 11, 14, 16, 56, 45)}  # datetime.fromtimestamp(1573750605)
            }
        ]

    def testParse(self):
        expected = [curr_test["expected"] for curr_test in self.test_cases]
        observed = []
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as handle:
                handle.write(curr_test["content"] + "\n")
            res = RTAComplete(self.tmp_file)
            observed.append({"RTA_version": res.RTA_version, "end_date": res.end_date})
        self.assertEqual(expected, observed)


class TestRunInfo(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_RunInfo.xml")
        self.test_cases = [
            {  # NexSeq dual index
                "content": '''<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="4">
  <Run Id="191113_NB551452_0101_AHW7JJAFXY" Number="101">
    <Flowcell>HW7JJAFXY</Flowcell>
    <Instrument>NB551452</Instrument>
    <Date>191113</Date>
    <Reads>
      <Read Number="1" NumCycles="71" IsIndexedRead="N" />
      <Read Number="2" NumCycles="17" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="4" NumCycles="71" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="4" SurfaceCount="2" SwathCount="1" TileCount="12" SectionPerLane="3" LanePerSection="2">
      <TileSet TileNamingConvention="FiveDigit">
        <Tiles>
          <Tile>1_11101</Tile>
        </Tiles>
      </TileSet>
    </FlowcellLayout>
    <ImageDimensions Width="2592" Height="1944" />
    <ImageChannels>
      <Name>Red</Name>
      <Name>Green</Name>
    </ImageChannels>
  </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': 'HW7JJAFXY', 'layout': {'LaneCount': '4', 'SurfaceCount': '2', 'SwathCount': '1', 'TileCount': '12', 'SectionPerLane': '3', 'LanePerSection': '2'}},
                    "instrument": {'id': 'NB551452', 'platform': 'NextSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 71}, {'is_index': True, 'nb_cycles': 17}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 71}],
                    "run": {'number': '101', 'id': '191113_NB551452_0101_AHW7JJAFXY', 'start_date': datetime(2019, 11, 13, 0, 0)}
                },
            },
            {  # NexSeq dual index
                "content": '''<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="4">
  <Run Id="200316_NDX550421_RUO_0003_AHJK57BGXF" Number="3">
    <Flowcell>HJK57BGXF</Flowcell>
    <Instrument>NDX550421_RUO</Instrument>
    <Date>200316</Date>
    <Reads>
      <Read Number="1" NumCycles="151" IsIndexedRead="N" />
      <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="4" NumCycles="151" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="4" SurfaceCount="2" SwathCount="3" TileCount="12" SectionPerLane="3" LanePerSection="2">
      <TileSet TileNamingConvention="FiveDigit">
        <Tiles>
        </Tiles>
      </TileSet>
    </FlowcellLayout>
    <ImageDimensions Width="2592" Height="1944" />
    <ImageChannels>
      <Name>Red</Name>
      <Name>Green</Name>
    </ImageChannels>
  </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': 'HJK57BGXF', 'layout': {'LaneCount': '4', 'SurfaceCount': '2', 'SwathCount': '3', 'TileCount': '12', 'SectionPerLane': '3', 'LanePerSection': '2'}},
                    "instrument": {'id': 'NDX550421_RUO', 'platform': 'NextSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {'number': '3', 'id': '200316_NDX550421_RUO_0003_AHJK57BGXF', 'start_date': datetime(2020, 3, 16, 0, 0)}
                },
            },
            {  # MiSeq dual index
                "content": '''<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
  <Run Id="191120_M70265_0340_000000000-CRRBF" Number="340">
    <Flowcell>000000000-CRRBF</Flowcell>
    <Instrument>M70265</Instrument>
    <Date>191120</Date>
    <Reads>
      <Read NumCycles="151" Number="1" IsIndexedRead="N" />
      <Read NumCycles="8" Number="2" IsIndexedRead="Y" />
      <Read NumCycles="8" Number="3" IsIndexedRead="Y" />
      <Read NumCycles="151" Number="4" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="1" TileCount="19" />
  </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': '000000000-CRRBF', 'layout': {'LaneCount': '1', 'SurfaceCount': '2', 'SwathCount': '1', 'TileCount': '19'}},
                    "instrument": {'id': 'M70265', 'platform': 'MiSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {'number': '340', 'id': '191120_M70265_0340_000000000-CRRBF', 'start_date': datetime(2019, 11, 20, 0, 0)}
                }
            },
            {  # HiSeq 2500 single index
                "content": '''<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
  <Run Id="151014_D00154_0539_AHCTK2BCXX_HiSeq 2500_TruSeq Exome_12 replicates of NA12878" Number="539">
    <Flowcell>HCTK2BCXX</Flowcell>
    <Instrument>D00154</Instrument>
    <Date>151014</Date>
    <Reads>
      <Read Number="1" NumCycles="76" IsIndexedRead="N" />
      <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="76" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="2" SurfaceCount="2" SwathCount="2" TileCount="16" />
    <AlignToPhiX>
      <Lane>1</Lane>
      <Lane>2</Lane>
    </AlignToPhiX>
  </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': 'HCTK2BCXX', 'layout': {'LaneCount': '2', 'SurfaceCount': '2', 'SwathCount': '2', 'TileCount': '16'}},
                    "instrument": {'id': 'D00154', 'platform': 'HiSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 76}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 76}],
                    "run": {'number': '539', 'id': '151014_D00154_0539_AHCTK2BCXX_HiSeq 2500_TruSeq Exome_12 replicates of NA12878', 'start_date': datetime(2015, 10, 14, 0, 0)}
                }
            },
            {  # MiniSeq dual index
                "content": '''<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="4">
  <Run Id="160822_ML-P2-01_0042_A000H02MNK" Number="42">
    <Flowcell>000H02MNK</Flowcell>
    <Instrument>ML-P2-01</Instrument>
    <Date>160822</Date>
    <Reads>
      <Read Number="1" NumCycles="150" IsIndexedRead="N" />
      <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="4" NumCycles="150" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="3" TileCount="10" SectionPerLane="1" LanePerSection="1">
      <TileSet TileNamingConvention="FiveDigit">
        <Tiles>
          <Tile>1_11102</Tile>
          <Tile>1_21102</Tile>
        </Tiles>
      </TileSet>
    </FlowcellLayout>
    <ImageDimensions Width="2592" Height="1944" />
    <ImageChannels>
      <Name>Red</Name>
      <Name>Green</Name>
    </ImageChannels>
  </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': '000H02MNK', 'layout': {'LaneCount': '1', 'SurfaceCount': '2', 'SwathCount': '3', 'TileCount': '10', 'SectionPerLane': '1', 'LanePerSection': '1'}},
                    "instrument": {'id': 'ML-P2-01', 'platform': 'MiniSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 150}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 150}],
                    "run": {'number': '42', 'id': '160822_ML-P2-01_0042_A000H02MNK', 'start_date': datetime(2016, 8, 22, 0, 0)}
                }
            },
            {  # NovaSeq dual index from (https://github.com/broadinstitute/viral-ngs/tree/master/test/input/TestRunInfo)
                "content": '''<?xml version="1.0" encoding="utf-8"?>
<RunInfo>
    <Run Id="180627_SL-NVB_0098_AFCHCYTJDMXX" Number="98">
        <Flowcell>HCYTJDMXX</Flowcell>
        <Instrument>A00198</Instrument>
        <Date>6/27/2018 4:59:20 PM</Date>
        <Reads>
            <Read Number="1" NumCycles="101" IsIndexedRead="N"/>
            <Read Number="2" NumCycles="8" IsIndexedRead="Y"/>
            <Read Number="3" NumCycles="8" IsIndexedRead="Y"/>
            <Read Number="4" NumCycles="101" IsIndexedRead="N"/>
        </Reads>
        <FlowcellLayout LaneCount="2" SurfaceCount="2" SwathCount="4" TileCount="88" FlowcellSide="1">
            <TileSet TileNamingConvention="FourDigit">
                <Tiles>
                    <Tile>1_2101</Tile>
                </Tiles>
            </TileSet>
        </FlowcellLayout>
        <AlignToPhiX/>
        <ImageDimensions Width="3200" Height="3607"/>
        <ImageChannels>
            <Name>RED</Name>
            <Name>GREEN</Name>
        </ImageChannels>
    </Run>
</RunInfo>''',
                "expected": {
                    "flowcell": {'id': 'HCYTJDMXX', 'layout': {'LaneCount': '2', 'SurfaceCount': '2', 'SwathCount': '4', 'TileCount': '88', 'FlowcellSide': '1'}},
                    "instrument": {'id': 'A00198', 'platform': 'NovaSeq'},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 101}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 101}],
                    "run": {'number': '98', 'id': '180627_SL-NVB_0098_AFCHCYTJDMXX', 'start_date': datetime(2018, 6, 27, 16, 59, 20)}
                }
            }
        ]

    def testParse(self):
        expected = [curr_test["expected"] for curr_test in self.test_cases]
        observed = []
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as handle:
                handle.write(curr_test["content"] + "\n")
            res = RunInfo(self.tmp_file)
            observed.append({
                "flowcell": res.flowcell,
                "instrument": res.instrument,
                "reads_phases": res.reads_phases,
                "run": res.run
            })
        self.assertEqual(expected, observed)


class TestRunParameters(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_RunParameters.xml")
        self.test_cases = [
            {  # NexSeq with post-process
                "content": '''<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <Setup>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <ApplicationVersion>1.3.0.23</ApplicationVersion>
    <ApplicationName>NextSeq Control Software</ApplicationName>
    <NumTilesPerSwath>12</NumTilesPerSwath>
    <NumSwaths>3</NumSwaths>
    <NumLanes>4</NumLanes>
    <Read1>151</Read1>
    <Read2>151</Read2>
    <Index1Read>6</Index1Read>
    <Index2Read>0</Index2Read>
    <SectionPerLane>3</SectionPerLane>
    <LanePerSection>2</LanePerSection>
  </Setup>
  <RunID>141107_NS500413_0033_H14U5BGXX</RunID>
  <InstrumentID>NS500413</InstrumentID>
  <RunNumber>33</RunNumber>
  <RTAVersion>2.1.3</RTAVersion>
  <SystemSuiteVersion>1.3.0.9</SystemSuiteVersion>
  <FlowCellSerial>H14U5BGXX</FlowCellSerial>
  <PR2BottleSerial>NS2041161-BUFFR</PR2BottleSerial>
  <ReagentKitSerial>NS2044178-REAGT</ReagentKitSerial>
  <ReagentKitSerialWasEnteredInBaseSpace>false</ReagentKitSerialWasEnteredInBaseSpace>
  <ExperimentName>NEXTSEQ_HU01_CPCT_CPCT34WES</ExperimentName>
  <LibraryID>LibrariesPool_CPCT34WES</LibraryID>
  <Chemistry>NextSeq High</Chemistry>
  <SelectedTiles>
    <Tile>1_11101</Tile>
    <Tile>2_23612</Tile>
  </SelectedTiles>
  <RunFolder>D:\\Illumina\\NextSeq Control Software Temp\\141107_NS500413_0033_H14U5BGXX\\</RunFolder>
  <OutputFolder>\\NEXTSEQ500_HU01\\141107_NS500413_0033_H14U5BGXX\\</OutputFolder>
  <RecipeFolder>C:\\Illumina\\NextSeq Control Software\\Recipe\\High\\v1.1</RecipeFolder>
  <SimulationFolder />
  <RunStartDate>141107</RunStartDate>
  <FocusMethod>IXFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>false</SaveFocusImages>
  <SaveScanImages>false</SaveScanImages>
  <SelectiveSave>true</SelectiveSave>
  <IsPairedEnd>true</IsPairedEnd>
  <AnalysisWorkflowType>FastQ</AnalysisWorkflowType>
  <CustomReadOnePrimer>BP10</CustomReadOnePrimer>
  <CustomReadTwoPrimer>BP11</CustomReadTwoPrimer>
  <CustomIndexPrimer>BP12</CustomIndexPrimer>
  <CustomIndexTwoPrimer>EXT</CustomIndexTwoPrimer>
  <UsesCustomReadOnePrimer>false</UsesCustomReadOnePrimer>
  <UsesCustomReadTwoPrimer>false</UsesCustomReadTwoPrimer>
  <UsesCustomIndexPrimer>false</UsesCustomIndexPrimer>
  <UsesCustomIndexTwoPrimer>false</UsesCustomIndexTwoPrimer>
  <RunManagementType>BaseSpaceCloud</RunManagementType>
  <BaseSpaceRunId>9092087</BaseSpaceRunId>
  <BaseSpaceRunMode>UserAndIlluminaHealth</BaseSpaceRunMode>
  <ComputerName>NEXTSEQ500_HU01</ComputerName>
  <SequencingStarted>true</SequencingStarted>
  <PlannedRead1Cycles>151</PlannedRead1Cycles>
  <PlannedRead2Cycles>151</PlannedRead2Cycles>
  <PlannedIndex1ReadCycles>6</PlannedIndex1ReadCycles>
  <PlannedIndex2ReadCycles>0</PlannedIndex2ReadCycles>
  <IsRehyb>false</IsRehyb>
  <MaxCyclesSupportedByReagentKit>318</MaxCyclesSupportedByReagentKit>
</RunParameters>''',
                "expected": {
                    "instrument": {"id": "NS500413", "platform": "NextSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 6}, {'is_index': True, 'nb_cycles': 0}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {"number": "33", "id": "141107_NS500413_0033_H14U5BGXX", "start_date": datetime(2014, 11, 7)},
                    "kit": {"flowcell_id": "H14U5BGXX", "reagent_kit_id": "NS2044178-REAGT"},
                    "post_process": "FastQ",
                    "software": {"RTA": "2.1.3", "CS": "1.3.0.23"}
                }
            },
            {  # NexSeq without post-process
                "content": '''<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <RunParametersVersion>NextSeq_4_0_0</RunParametersVersion>
  <Setup>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <ApplicationVersion>4.0.1.41</ApplicationVersion>
    <ApplicationName>NextSeq Control Software</ApplicationName>
    <NumTilesPerSwath>12</NumTilesPerSwath>
    <NumSwaths>1</NumSwaths>
    <NumLanes>4</NumLanes>
    <Read1>151</Read1>
    <Read2>151</Read2>
    <Index1Read>8</Index1Read>
    <Index2Read>8</Index2Read>
    <SectionPerLane>3</SectionPerLane>
    <LanePerSection>2</LanePerSection>
  </Setup>
  <RunID>200318_NDX550421_RUO_0004_AH5THYAFX2</RunID>
  <CopyServiceRunId>16A9530EA81D22F4</CopyServiceRunId>
  <InstrumentID>NDX550421_RUO</InstrumentID>
  <RunNumber>4</RunNumber>
  <RTAVersion>2.11.3</RTAVersion>
  <SystemSuiteVersion>4.0.1.41</SystemSuiteVersion>
  <LocalRunManagerVersion>2.4.0.2385</LocalRunManagerVersion>
  <RecipeVersion>4.0.0.1</RecipeVersion>
  <FirmwareVersion>r.2019.2.0</FirmwareVersion>
  <FlowCellRfidTag>
    <SerialNumber>H5THYAFX2</SerialNumber>
    <PartNumber>20022409</PartNumber>
    <LotNumber>20428008</LotNumber>
    <ExpirationDate>2022-02-04T00:00:00</ExpirationDate>
  </FlowCellRfidTag>
  <PR2BottleRfidTag>
    <SerialNumber>NS4071554-BUFFR</SerialNumber>
    <PartNumber>15057941</PartNumber>
    <LotNumber>20429662</LotNumber>
    <ExpirationDate>2021-02-09T00:00:00</ExpirationDate>
  </PR2BottleRfidTag>
  <ReagentKitRfidTag>
    <SerialNumber>NS4065291-REAGT</SerialNumber>
    <PartNumber>15057939</PartNumber>
    <LotNumber>20426104</LotNumber>
    <ExpirationDate>2021-02-04T00:00:00</ExpirationDate>
  </ReagentKitRfidTag>
  <FlowCellSerial>H5THYAFX2</FlowCellSerial>
  <PR2BottleSerial>NS4071554-BUFFR</PR2BottleSerial>
  <ReagentKitSerial>NS4065291-REAGT</ReagentKitSerial>
  <ReagentKitSerialWasEnteredInBaseSpace>false</ReagentKitSerialWasEnteredInBaseSpace>
  <ExperimentName>200318</ExperimentName>
  <LibraryID />
  <StateDescription />
  <Chemistry>NextSeq Mid</Chemistry>
  <ChemistryVersion>2.1</ChemistryVersion>
  <SelectedTiles>
    <Tile>1_11101</Tile>
    <Tile>4_21612</Tile>
  </SelectedTiles>
  <RunFolder>D:\\Illumina\\NextSeq Control Software Temp\\200318_NDX550421_RUO_0004_AH5THYAFX2\\</RunFolder>
  <RTALogsFolder />
  <PreRunFolderRoot>D:\\Illumina\\NextSeq Control Software PreRun</PreRunFolderRoot>
  <PreRunFolder>D:\\Illumina\\NextSeq Control Software PreRun\\NDX550421_RUO_2020-03-18__16_36_39</PreRunFolder>
  <OutputFolder>\\NextSeq3\\200318_NDX550421_RUO_0004_AH5THYAFX2\\</OutputFolder>
  <RecipeFolder>C:\\Program Files\\Illumina\\NextSeq Control Software\\Recipe\\Mid\\v2.5</RecipeFolder>
  <SimulationFolder />
  <RunStartDate>200318</RunStartDate>
  <BaseSpaceUserName />
  <LocalRunManagerUserName />
  <FocusMethod>IXFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>false</SaveFocusImages>
  <SaveScanImages>false</SaveScanImages>
  <SelectiveSave>true</SelectiveSave>
  <IsPairedEnd>true</IsPairedEnd>
  <AnalysisWorkflowType />
  <CustomReadOnePrimer>BP10</CustomReadOnePrimer>
  <CustomReadTwoPrimer>BP11</CustomReadTwoPrimer>
  <CustomIndexOnePrimer>BP14</CustomIndexOnePrimer>
  <CustomIndexTwoPrimer>BP14</CustomIndexTwoPrimer>
  <UsesCustomReadOnePrimer>false</UsesCustomReadOnePrimer>
  <UsesCustomReadTwoPrimer>false</UsesCustomReadTwoPrimer>
  <UsesCustomIndexPrimer>false</UsesCustomIndexPrimer>
  <UsesCustomIndexTwoPrimer>false</UsesCustomIndexTwoPrimer>
  <BaseSpaceRunId>194270150</BaseSpaceRunId>
  <LocalRunManagerRunId xsi:nil="true" />
  <RunSetupType>ManualRunSetup</RunSetupType>
  <RunMode>PerformanceDataOnly</RunMode>
  <ComputerName>NDX550421-RUO</ComputerName>
  <SequencingStarted>true</SequencingStarted>
  <PlannedRead1Cycles>151</PlannedRead1Cycles>
  <PlannedRead2Cycles>151</PlannedRead2Cycles>
  <PlannedIndex1ReadCycles>8</PlannedIndex1ReadCycles>
  <PlannedIndex2ReadCycles>8</PlannedIndex2ReadCycles>
  <IsRehyb>false</IsRehyb>
  <PurgeConsumables>true</PurgeConsumables>
  <MaxCyclesSupportedByReagentKit>318</MaxCyclesSupportedByReagentKit>
  <ModuleName />
  <ModuleVersion />
</RunParameters>''',
                "expected": {
                    "instrument": {"id": "NDX550421_RUO", "platform": "NextSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {"number": "4", "id": "200318_NDX550421_RUO_0004_AH5THYAFX2", "start_date": datetime(2020, 3, 18)},
                    "kit": {"flowcell_id": "H5THYAFX2", "reagent_kit_id": "NS4065291-REAGT"},
                    "post_process": None,
                    "software": {"RTA": "2.11.3", "CS": "4.0.1.41"}
                }
            },
            {  # MiSeq
                "content": '''<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EnableCloud>true</EnableCloud>
  <RunParametersVersion>MiSeq_1_1</RunParametersVersion>
  <CopyManifests>true</CopyManifests>
  <FlowcellRFIDTag>
    <SerialNumber>000000000-AAC1Y</SerialNumber>
    <PartNumber>15028382</PartNumber>
    <ExpirationDate>2015-05-28T00:00:00</ExpirationDate>
  </FlowcellRFIDTag>
  <PR2BottleRFIDTag>
    <SerialNumber>MS2313945-00PR2</SerialNumber>
    <PartNumber>15041807</PartNumber>
    <ExpirationDate>2015-05-18T00:00:00</ExpirationDate>
  </PR2BottleRFIDTag>
  <ReagentKitRFIDTag>
    <SerialNumber>MS2740035-300V2</SerialNumber>
    <PartNumber>15033572</PartNumber>
    <ExpirationDate>2015-06-12T00:00:00</ExpirationDate>
  </ReagentKitRFIDTag>
  <Resumable>true</Resumable>
  <ManifestFiles />
  <AfterRunWashMethod>Post-Run Wash</AfterRunWashMethod>
  <Setup>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <ApplicationVersion>2.5.0.5</ApplicationVersion>
    <NumTilesPerSwath>14</NumTilesPerSwath>
    <NumSwaths>1</NumSwaths>
    <NumLanes>1</NumLanes>
    <ApplicationName>MiSeq Control Software</ApplicationName>
  </Setup>
  <RunID>141105_M01102_0149_000000000-AAC1Y</RunID>
  <ScannerID>M01102</ScannerID>
  <RunNumber>149</RunNumber>
  <FPGAVersion>9.5.12</FPGAVersion>
  <MCSVersion>2.5.0.5</MCSVersion>
  <RTAVersion>1.18.54</RTAVersion>
  <Barcode>000000000-AAC1Y</Barcode>
  <PR2BottleBarcode>MS2313945-00PR2</PR2BottleBarcode>
  <ReagentKitPartNumberEntered>15033572</ReagentKitPartNumberEntered>
  <ReagentKitVersion>Version2</ReagentKitVersion>
  <ReagentKitBarcode>MS2740035-300V2</ReagentKitBarcode>
  <PreviousReagentKitBarcode>MS2740035-300V2</PreviousReagentKitBarcode>
  <ExperimentName>MISEQ_UMC01_CVERMEULEN_WDL004</ExperimentName>
  <Chemistry>Amplicon</Chemistry>
  <Username>sbsuser</Username>
  <Workflow>
    <Analysis>GenerateFASTQ</Analysis>
  </Workflow>
  <EnableAnalysis>false</EnableAnalysis>
  <Reads>
    <RunInfoRead Number="1" NumCycles="151" IsIndexedRead="N" />
    <RunInfoRead Number="2" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="3" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="4" NumCycles="151" IsIndexedRead="N" />
  </Reads>
  <TempFolder>D:\\Illumina\\MiSeqTemp\\141105_M01102_0149_000000000-AAC1Y</TempFolder>
  <AnalysisFolder>D:\\Illumina\\MiSeqAnalysis\\141105_M01102_0149_000000000-AAC1Y</AnalysisFolder>
  <RunStartDate>141105</RunStartDate>
  <MostRecentWashType>PostRun</MostRecentWashType>
  <RecipeFolder>D:\\Illumina\\MiSeq Control Software\\CustomRecipe</RecipeFolder>
  <ILMNOnlyRecipeFolder>C:\\Illumina\MiSeq Control Software\\Recipe</ILMNOnlyRecipeFolder>
  <SampleSheetName>MISEQ_UMC01_CVERMEULEN_WDL004</SampleSheetName>
  <SampleSheetFolder>C:\\Users\\sbsuser\\AppData\\Roaming\\Illumina\\Illumina Experiment Manager\\Sample Sheets</SampleSheetFolder>
  <ManifestFolder>D:\\Illumina\\MiSeq Control Software\\Manifests</ManifestFolder>
  <OutputFolder>\\results_miseq_umc01\\141105_M01102_0149_000000000-AAC1Y</OutputFolder>
  <FocusMethod>AutoFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>true</SaveFocusImages>
  <SaveScanImages>true</SaveScanImages>
  <CloudUsername>user@institute.eu</CloudUsername>
  <RunManagementType>Standalone</RunManagementType>
  <CloudRunId>9058072</CloudRunId>
  <SendInstrumentHealthToILMN>false</SendInstrumentHealthToILMN>
</RunParameters>''',
                "expected": {
                    "instrument": {"id": "M01102", "platform": "MiSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {"number": "149", "id": "141105_M01102_0149_000000000-AAC1Y", "start_date": datetime(2014, 11, 5)},
                    "kit": {"flowcell_id": "000000000-AAC1Y", "reagent_kit_id": "MS2740035-300V2"},
                    "post_process": "GenerateFASTQ",
                    "software": {"RTA": "1.18.54", "CS": "2.5.0.5"}
                }
            },
            {  # HiSeq
                "content": '''<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <Setup>
    <ExperimentName>HiSeq_UMC01_B_Diag_HAMVRADXX</ExperimentName>
    <ScanID>D00267</ScanID>
    <FCPosition>B</FCPosition>
    <WorkFlowType>SINGLEINDEX</WorkFlowType>
    <PairEndFC>true</PairEndFC>
    <Read1>101</Read1>
    <IndexRead1>7</IndexRead1>
    <IndexRead2>0</IndexRead2>
    <Read2>101</Read2>
    <OutputFolder>Z:\\</OutputFolder>
    <CompressBcls>true</CompressBcls>
    <RemapQScores />
    <NearNeighborCorrection>false</NearNeighborCorrection>
    <PeriodicSave>Save All Thumbnails</PeriodicSave>
    <Flowcell>HiSeq Rapid Flow Cell v1</Flowcell>
    <FirstBaseConfirmation>false</FirstBaseConfirmation>
    <SampleSheet>C:\\Users\\sbsuser\\Desktop\\sample sheets\\Diagnostic Sample Sheets\\HiSeq_UMC01_B_Diag_HAMVRADXX.csv</SampleSheet>
    <KeepIntensityFiles>false</KeepIntensityFiles>
    <Sbs>TruSeq Rapid SBS Kit</Sbs>
    <Pe>TruSeq Rapid PE Cluster Kit</Pe>
    <Index>TruSeq Rapid PE Cluster Kit</Index>
    <AlignToPhiX>
      <Lane>1</Lane>
      <Lane>2</Lane>
    </AlignToPhiX>
    <ClusteringChoice>OnBoardClustering</ClusteringChoice>
    <RunMode>RapidRun</RunMode>
    <Rehyb>None</Rehyb>
    <PerformPreRunFluidicsCheck>false</PerformPreRunFluidicsCheck>
    <ServiceRun>false</ServiceRun>
    <ApplicationName>HiSeq Control Software</ApplicationName>
    <ApplicationVersion>2.2.38</ApplicationVersion>
    <RunID>141110_D00267_0193_BHAMVRADXX</RunID>
    <RunStartDate>141110</RunStartDate>
    <IntegrationMode>BaseSpace</IntegrationMode>
    <BaseSpaceSettings>
      <RunMonitoringOnly>true</RunMonitoringOnly>
      <Username>user@institute.nl</Username>
      <PlannedRun>false</PlannedRun>
      <RunId>9136142</RunId>
      <TempFolder>E:\\Illumina\\BaseSpaceTemp</TempFolder>
      <SendInstrumentHealthToILMN>true</SendInstrumentHealthToILMN>
    </BaseSpaceSettings>
    <ScannerID>D00267</ScannerID>
    <ScanNumber>193</ScanNumber>
    <ComputerName>HWI-D00267</ComputerName>
    <FPGAVersion>7.9.7</FPGAVersion>
    <CPLDVersion>3.0.0</CPLDVersion>
    <RTAVersion>1.18.61</RTAVersion>
    <ChemistryVersion>Illumina,Bruno Fluidics Controller,0,v2.0340</ChemistryVersion>
    <CameraFirmware>2.01-F20-R02</CameraFirmware>
    <CameraDriver>6.45.20.3690</CameraDriver>
    <FocusCameraFirmware />
    <Barcode>HAMVRADXX</Barcode>
    <Username>sbsuser</Username>
    <SelectedSections>
      <Section Name="A_1" />
      <Section Name="B_1" />
    </SelectedSections>
    <FocusMethod>DynamicITF</FocusMethod>
    <SelectedSurface>BothLaneSurfaces</SelectedSurface>
    <SwathScanMode>DualSwathFC</SwathScanMode>
    <EnableLft>true</EnableLft>
    <AutoTiltOnce>true</AutoTiltOnce>
    <EnableAutoCenter>true</EnableAutoCenter>
    <EnableAnalysis>true</EnableAnalysis>
    <EnableBasecalling>true</EnableBasecalling>
    <EnableCameraLogging>false</EnableCameraLogging>
    <AdapterPlate>HiSeq Adapter Plate</AdapterPlate>
    <SlideHolder>HiSeq Flow Cell Holder</SlideHolder>
    <TemplateCycleCount>4</TemplateCycleCount>
    <NumAnalysisThreads>8</NumAnalysisThreads>
    <FPGADynamicFocusSettings>
      <MaxInitialZJumpHalfUm>3</MaxInitialZJumpHalfUm>
      <MaxSubsequentZJumpHalfUm>7</MaxSubsequentZJumpHalfUm>
      <NumberOfInitialZJumps>0</NumberOfInitialZJumps>
      <CVGainStart>500</CVGainStart>
      <CVGainPosLocked>500</CVGainPosLocked>
      <Offset>250</Offset>
      <HotPixel>350</HotPixel>
      <MotorDelayFrames>10</MotorDelayFrames>
      <SoftwareLaserLag>200</SoftwareLaserLag>
      <DitherSize>100</DitherSize>
      <GroupSize>50</GroupSize>
      <DitherShift>0</DitherShift>
      <IntensityCeiling>65535</IntensityCeiling>
      <IGain>100</IGain>
      <IHistory>4</IHistory>
    </FPGADynamicFocusSettings>
    <TileWidth>2048</TileWidth>
    <TileHeight>10000</TileHeight>
    <ImageWidth>2048</ImageWidth>
    <ImageHeight>160000</ImageHeight>
    <AreaPerPixelmm2>1.40625E-07</AreaPerPixelmm2>
    <LaneLength>60</LaneLength>
    <NumTilesPerSwath>16</NumTilesPerSwath>
    <NumSwaths>2</NumSwaths>
    <UseExistingRecipe>false</UseExistingRecipe>
    <Reads>
      <Read Number="1" NumCycles="101" IsIndexedRead="N" />
      <Read Number="2" NumCycles="7" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="101" IsIndexedRead="N" />
    </Reads>
    <EnableNotifications>false</EnableNotifications>
    <ReagentKits>
      <Sbs>
        <SbsReagentKit>
          <ID>9806021</ID>
          <Prime>false</Prime>
          <NumberCyclesRemaining>225</NumberCyclesRemaining>
          <IsNew50Cycle>false</IsNew50Cycle>
          <IsNew200Cycle>true</IsNew200Cycle>
        </SbsReagentKit>
      </Sbs>
      <Index>
        <ReagentKit>
          <ID>9991543</ID>
        </ReagentKit>
      </Index>
      <Pe />
      <Rehyb />
    </ReagentKits>
    <ReagentBottles>
      <Sbs />
    </ReagentBottles>
    <Resume>false</Resume>
    <ResumeCycle>0</ResumeCycle>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <TempFolder>E:\\Illumina\\HiSeqTemp\\141110_D00267_0193_BHAMVRADXX</TempFolder>
    <RecipeFragmentVersion>1.5.14.0</RecipeFragmentVersion>
    <PromptForPeReagents>false</PromptForPeReagents>
    <MockRun>false</MockRun>
    <ScannedBarcode />
  </Setup>
  <Version>1</Version>
</RunParameters>''',
                "expected": {
                    "instrument": {"id": "D00267", "platform": "HiSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 101}, {'is_index': True, 'nb_cycles': 7}, {'is_index': False, 'nb_cycles': 101}],
                    "run": {"number": "193", "id": "141110_D00267_0193_BHAMVRADXX", "start_date": datetime(2014, 11, 10)},
                    "kit": {"flowcell_id": "HAMVRADXX", "reagent_kit_id": None},
                    "post_process": None,
                    "software": {"RTA": "1.18.61", "CS": "2.2.38"}
                }
            }
        ]

    def testParse(self):
        expected = [curr_test["expected"] for curr_test in self.test_cases]
        observed = []
        for curr_test in self.test_cases:
            with open(self.tmp_file, "w") as handle:
                handle.write(curr_test["content"] + "\n")
            res = RunParameters(self.tmp_file)
            observed.append({
                "instrument": res.instrument,
                "reads_phases": res.reads_phases,
                "run": res.run,
                "kit": res.kit,
                "post_process": res.post_process,
                "software": res.software
            })
        self.assertEqual(expected, observed)


class TestFunctions(unittest.TestCase):
    def testGetInfFromSeqID(self):
        # Whithout UMI
        expected = {
            "sequencer_id": "NDX550421_RUO",
            "run_id": "15",
            "flowcell_id": "HY53NBGXC",
            "lane_id": 1,
            "tile_id": 11101,
            "x_pos": 9518,
            "y_pos": 1037,
            "umi": None
        }
        observed = getInfFromSeqID("NDX550421_RUO:15:HY53NBGXC:1:11101:9518:1037")
        self.assertEqual(expected, observed)
        # With UMI
        expected = {
            "sequencer_id": "NDX550421_RUO",
            "run_id": "15",
            "flowcell_id": "HY53NBGXC",
            "lane_id": 1,
            "tile_id": 11101,
            "x_pos": 9518,
            "y_pos": 1037,
            "umi": "TTGCCNG+CCTCATA"
        }
        observed = getInfFromSeqID("NDX550421_RUO:15:HY53NBGXC:1:11101:9518:1037:TTGCCNG+CCTCATA")
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
