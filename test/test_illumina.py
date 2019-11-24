#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
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

from anacore.illumina import RTAComplete, RunInfo


########################################################################
#
# FUNCTIONS
#
########################################################################
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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
