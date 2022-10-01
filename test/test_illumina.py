#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.5.0'
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

from anacore.illumina import Bcl2fastqLog, DemultStat, getInfFromSeqID, RTAComplete, Run, RunInfo, RunParameters


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


class TestDemultStat(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_demuliplexStats.json")
        with open(self.tmp_file, "w") as writer:
            writer.write('{"Flowcell":"HFF5HBGXL","RunNumber":256,"RunId":"220620_NDX550421_RUO_0256_AHFF5HBGXL","ReadInfosForLanes":[{"LaneNumber":1,"ReadInfos":[{"Number":1,"NumCycles":101,"IsIndexedRead":false},{"Number":1,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":101,"IsIndexedRead":false}]},{"LaneNumber":2,"ReadInfos":[{"Number":1,"NumCycles":101,"IsIndexedRead":false},{"Number":1,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":101,"IsIndexedRead":false}]}],"ConversionResults":[{"LaneNumber":1,"TotalClustersRaw":118154091,"TotalClustersPF":110173898,"Yield":22255127396,"DemuxResults":[{"SampleId":"splA-DNA","SampleName":"splA-DNA","IndexMetrics":[{"IndexSequence":"TTAATCAG+CTTCGCCT","MismatchCounts":{"0":9478792,"1":141994}}],"NumberReads":9620786,"Yield":1943398772,"ReadMetrics":[{"ReadNumber":1,"Yield":971699386,"YieldQ30":924128846,"QualityScoreSum":33738537475,"TrimmedBases":48420782},{"ReadNumber":2,"Yield":971699386,"YieldQ30":918562030,"QualityScoreSum":33615137854,"TrimmedBases":47944431}]},{"SampleId":"splA-RNA","SampleName":"splA-RNA","IndexMetrics":[{"IndexSequence":"TCCGGAGA+AGGATAGG","MismatchCounts":{"0":2699794,"1":40900}}],"NumberReads":2740694,"Yield":553620188,"ReadMetrics":[{"ReadNumber":1,"Yield":276810094,"YieldQ30":261295087,"QualityScoreSum":9569976358,"TrimmedBases":34397480},{"ReadNumber":2,"Yield":276810094,"YieldQ30":237377116,"QualityScoreSum":9065991028,"TrimmedBases":31656823}]},{"SampleId":"splB-DNA","SampleName":"splB-DNA","IndexMetrics":[{"IndexSequence":"CGCTCATT+TAAGATTA","MismatchCounts":{"0":9564216,"1":131182}}],"NumberReads":9695398,"Yield":1958470396,"ReadMetrics":[{"ReadNumber":1,"Yield":979235198,"YieldQ30":935923618,"QualityScoreSum":34097986656,"TrimmedBases":37858895},{"ReadNumber":2,"Yield":979235198,"YieldQ30":929589854,"QualityScoreSum":33958376843,"TrimmedBases":37570166}]},{"SampleId":"splB-RNA","SampleName":"splB-RNA","IndexMetrics":[{"IndexSequence":"CTGAAGCT+TCAGAGCC","MismatchCounts":{"0":2597339,"1":49487}}],"NumberReads":2646826,"Yield":534658852,"ReadMetrics":[{"ReadNumber":1,"Yield":267329426,"YieldQ30":254469868,"QualityScoreSum":9286867469,"TrimmedBases":8407441},{"ReadNumber":2,"Yield":267329426,"YieldQ30":250039144,"QualityScoreSum":9191210821,"TrimmedBases":6605292}]},{"SampleId":"splC-DNA","SampleName":"splC-DNA","IndexMetrics":[{"IndexSequence":"TCCGCGAA+AGTAAGTA","MismatchCounts":{"0":9152620,"1":121190}}],"NumberReads":9273810,"Yield":1873309620,"ReadMetrics":[{"ReadNumber":1,"Yield":936654810,"YieldQ30":892413065,"QualityScoreSum":32555394461,"TrimmedBases":31828023},{"ReadNumber":2,"Yield":936654810,"YieldQ30":885028001,"QualityScoreSum":32393432900,"TrimmedBases":31576068}]},{"SampleId":"splD-DNA","SampleName":"splD-DNA","IndexMetrics":[{"IndexSequence":"ATTACTCG+GACTTCCT","MismatchCounts":{"0":9922693,"1":138354}}],"NumberReads":10061047,"Yield":2032331494,"ReadMetrics":[{"ReadNumber":1,"Yield":1016165747,"YieldQ30":966795535,"QualityScoreSum":35290830751,"TrimmedBases":53576473},{"ReadNumber":2,"Yield":1016165747,"YieldQ30":959254789,"QualityScoreSum":35125522244,"TrimmedBases":53052819}]}],"Undetermined":{"NumberReads":10126886,"Yield":2045630972,"ReadMetrics":[{"ReadNumber":1,"Yield":1022815486,"YieldQ30":942308528,"QualityScoreSum":34857988295,"TrimmedBases":50185140},{"ReadNumber":2,"Yield":1022815486,"YieldQ30":919589580,"QualityScoreSum":34250317568,"TrimmedBases":47568306}]}},{"LaneNumber":2,"TotalClustersRaw":116151639,"TotalClustersPF":108452034,"Yield":21907310868,"DemuxResults":[{"SampleId":"splA-DNA","SampleName":"splA-DNA","IndexMetrics":[{"IndexSequence":"TTAATCAG+CTTCGCCT","MismatchCounts":{"0":9312034,"1":144771}}],"NumberReads":9456805,"Yield":1910274610,"ReadMetrics":[{"ReadNumber":1,"Yield":955137305,"YieldQ30":908955716,"QualityScoreSum":33173997089,"TrimmedBases":47382071},{"ReadNumber":2,"Yield":955137305,"YieldQ30":903901935,"QualityScoreSum":33061153007,"TrimmedBases":46935857}]},{"SampleId":"splA-RNA","SampleName":"splA-RNA","IndexMetrics":[{"IndexSequence":"TCCGGAGA+AGGATAGG","MismatchCounts":{"0":2653660,"1":41311}}],"NumberReads":2694971,"Yield":544384142,"ReadMetrics":[{"ReadNumber":1,"Yield":272192071,"YieldQ30":257193370,"QualityScoreSum":9415122745,"TrimmedBases":33684765},{"ReadNumber":2,"Yield":272192071,"YieldQ30":233971488,"QualityScoreSum":8925076315,"TrimmedBases":31059435}]},{"SampleId":"splB-DNA","SampleName":"splB-DNA","IndexMetrics":[{"IndexSequence":"CGCTCATT+TAAGATTA","MismatchCounts":{"0":9389668,"1":133733}}],"NumberReads":9523401,"Yield":1923727002,"ReadMetrics":[{"ReadNumber":1,"Yield":961863501,"YieldQ30":919992361,"QualityScoreSum":33505686273,"TrimmedBases":37143142},{"ReadNumber":2,"Yield":961863501,"YieldQ30":914031955,"QualityScoreSum":33373758590,"TrimmedBases":36872757}]},{"SampleId":"splB-RNA","SampleName":"splB-RNA","IndexMetrics":[{"IndexSequence":"CTGAAGCT+TCAGAGCC","MismatchCounts":{"0":2552513,"1":49897}}],"NumberReads":2602410,"Yield":525686820,"ReadMetrics":[{"ReadNumber":1,"Yield":262843410,"YieldQ30":250401718,"QualityScoreSum":9134727206,"TrimmedBases":8237441},{"ReadNumber":2,"Yield":262843410,"YieldQ30":246175047,"QualityScoreSum":9043525469,"TrimmedBases":6480551}]},{"SampleId":"splC-DNA","SampleName":"splC-DNA","IndexMetrics":[{"IndexSequence":"TCCGCGAA+AGTAAGTA","MismatchCounts":{"0":8997000,"1":121956}}],"NumberReads":9118956,"Yield":1842029112,"ReadMetrics":[{"ReadNumber":1,"Yield":921014556,"YieldQ30":878282672,"QualityScoreSum":32026485795,"TrimmedBases":31172049},{"ReadNumber":2,"Yield":921014556,"YieldQ30":871422340,"QualityScoreSum":31875534840,"TrimmedBases":30932132}]},{"SampleId":"splD-DNA","SampleName":"splD-DNA","IndexMetrics":[{"IndexSequence":"ATTACTCG+GACTTCCT","MismatchCounts":{"0":9755783,"1":141064}}],"NumberReads":9896847,"Yield":1999163094,"ReadMetrics":[{"ReadNumber":1,"Yield":999581547,"YieldQ30":951751473,"QualityScoreSum":34728541382,"TrimmedBases":52608236},{"ReadNumber":2,"Yield":999581547,"YieldQ30":944728842,"QualityScoreSum":34573985187,"TrimmedBases":52120720}]}],"Undetermined":{"NumberReads":10144184,"Yield":2049125168,"ReadMetrics":[{"ReadNumber":1,"Yield":1024562584,"YieldQ30":944586733,"QualityScoreSum":34933559142,"TrimmedBases":50457930},{"ReadNumber":2,"Yield":1024562584,"YieldQ30":922237756,"QualityScoreSum":34343749503,"TrimmedBases":47618888}]}}],"UnknownBarcodes":[{"Lane":1,"Barcodes":{"GGGGGGGG+AGATCTCG":3267780,"ACTGCTTA+AGATCTCG":304460,"ATGCGGCT+AGATCTCG":280440,"AAAAAAAT+TAGCCGCG":234920,"AAAAAAAT+TTCGTAGG":223900,"AAAAAAAT+AGTAAGTA":221040}},{"Lane":2,"Barcodes":{"GGGGGGGG+AGATCTCG":3216540,"ACTGCTTA+AGATCTCG":292100,"ATGCGGCT+AGATCTCG":276940,"AAAAAAAT+TAGCCGCG":228260,"AAAAAAAT+TTCGTAGG":220120,"AAAAAAAT+AGTAAGTA":216460}}]}')

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testSamplesCounts(self):
        expected = {
            "splA-DNA": 19077591,
            "splA-RNA": 5435665,
            "splB-DNA": 19218799,
            "splB-RNA": 5249236,
            "splC-DNA": 18392766,
            "splD-DNA": 19957894
        }
        stats = DemultStat(self.tmp_file)
        self.assertEqual(stats.samplesCounts(), expected)

    def testUndeterminedCounts(self):
        expected = {
            "GGGGGGGG+AGATCTCG": 6484320,
            "ACTGCTTA+AGATCTCG": 596560,
            "ATGCGGCT+AGATCTCG": 557380,
            "AAAAAAAT+TAGCCGCG": 463180,
            "AAAAAAAT+TTCGTAGG": 444020,
            "AAAAAAAT+AGTAAGTA": 437500
        }
        stats = DemultStat(self.tmp_file)
        self.assertEqual(stats.undeterminedCounts(), expected)

    def testUnexpectedBarcodes(self):
        expected = []
        stats = DemultStat(self.tmp_file)
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # RNA under
        for lane in stats.data["UnknownBarcodes"]:
            for barcode, ct in lane["Barcodes"].items():
                lane["Barcodes"][barcode] = ct * 10
        expected = [  # "GGGGGGGG+AGATCTCG" is skipped because it corresponds to without UDI
            {"spl": "ACTGCTTA+AGATCTCG", "ct": 5965600},
            {"spl": "ATGCGGCT+AGATCTCG", "ct": 5573800}
        ]
        self.assertEqual(stats.unexpectedBarcodes(), expected)
        # RNA under with skipped
        expected = [  # "GGGGGGGG+AGATCTCG" is skipped because it corresponds to without UDI
            {"spl": "ACTGCTTA+AGATCTCG", "ct": 5965600},
            {"spl": "ATGCGGCT+AGATCTCG", "ct": 5573800}
        ]
        self.assertEqual(stats.unexpectedBarcodes({"splA-RNA"}), expected)
        expected = []
        self.assertEqual(stats.unexpectedBarcodes({"splA-RNA", "splB-RNA"}), expected)


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

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

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

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

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


class TestRun(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_folder = os.path.join(tmp_folder, unique_id + "_run")
        os.makedirs(self.tmp_folder)

    def tearDown(self):
        if os.path.exists(self.tmp_folder):
            for filename in os.listdir(self.tmp_folder):
                os.remove(os.path.join(self.tmp_folder, filename))
            os.rmdir(self.tmp_folder)

    def testIsCopied(self):
        run = Run(self.tmp_folder)
        self.assertFalse(run.isCopied())
        with open(os.path.join(self.tmp_folder, "RTAComplete.txt"), "w") as writer:
            writer.write("")
        # HiSeq
        with open(os.path.join(self.tmp_folder, "runParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["HiSeq"])
        self.assertTrue(run.isCopied())
        os.remove(os.path.join(self.tmp_folder, "runParameters.xml"))
        # MiSeq
        with open(os.path.join(self.tmp_folder, "runParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["MiSeq"])
        self.assertTrue(run.isCopied())
        os.remove(os.path.join(self.tmp_folder, "runParameters.xml"))
        # NextSeq old NCS
        with open(os.path.join(self.tmp_folder, "RunParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["NextSeq_withPostProcess"])
        self.assertTrue(run.isCopied())
        # NextSeq NCS >= 4
        with open(os.path.join(self.tmp_folder, "RunParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["NextSeq_woutPostProcess"])
        self.assertFalse(run.isCopied())
        with open(os.path.join(self.tmp_folder, "CopyComplete.txt"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isCopied())
        os.remove(os.path.join(self.tmp_folder, "CopyComplete.txt"))

    def testIsEnded(self):
        run = Run(self.tmp_folder)
        with open(os.path.join(self.tmp_folder, "RTAComplete.txt"), "w") as writer:
            writer.write("")
        with open(os.path.join(self.tmp_folder, "CopyComplete.txt"), "w") as writer:
            writer.write("")
        # NextSeq without post-process
        with open(os.path.join(self.tmp_folder, "RunParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["NextSeq_woutPostProcess"])
        self.assertTrue(run.isEnded())
        os.remove(os.path.join(self.tmp_folder, "RunParameters.xml"))
        # NextSeq with post-process
        with open(os.path.join(self.tmp_folder, "RunParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["NextSeq_withPostProcess"])
        self.assertFalse(run.isEnded())
        with open(os.path.join(self.tmp_folder, "CompletedJobInfo.xml"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isEnded())
        os.remove(os.path.join(self.tmp_folder, "RunParameters.xml"))
        os.remove(os.path.join(self.tmp_folder, "CompletedJobInfo.xml"))
        # MiSeq
        with open(os.path.join(self.tmp_folder, "runParameters.xml"), "w") as writer:
            writer.write(TestRunParameters.CONTENT["MiSeq"])
        self.assertFalse(run.isEnded())
        with open(os.path.join(self.tmp_folder, "CompletedJobInfo.xml"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isEnded())
        os.remove(os.path.join(self.tmp_folder, "runParameters.xml"))
        os.remove(os.path.join(self.tmp_folder, "CompletedJobInfo.xml"))

    def testIsSequenced(self):
        run = Run(self.tmp_folder)
        self.assertFalse(run.isSequenced())
        with open(os.path.join(self.tmp_folder, "RTAComplete.txt"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isSequenced())


class TestRunParameters(unittest.TestCase):
    CONTENT = {
        "NextSeq_withPostProcess": '''<?xml version="1.0"?>
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
        "NextSeq_woutPostProcess": '''<?xml version="1.0"?>
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
        "MiSeq": '''<?xml version="1.0"?>
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
        "HiSeq": '''<?xml version="1.0"?>
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
            "MiSeq_invalid_RFID": """<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EnableCloud>false</EnableCloud>
  <RunParametersVersion>MiSeq_1_1</RunParametersVersion>
  <CopyManifests>true</CopyManifests>
  <FlowcellRFIDTag>
    <PartNumber>20660742</PartNumber>
    <ExpirationDate>0001-01-01T00:00:00</ExpirationDate>
  </FlowcellRFIDTag>
  <PR2BottleRFIDTag>
    <PartNumber>20663557</PartNumber>
    <ExpirationDate>0001-01-01T00:00:00</ExpirationDate>
  </PR2BottleRFIDTag>
  <ReagentKitRFIDTag>
    <SerialNumber>MS3133997-500V2</SerialNumber>
    <PartNumber>15033573</PartNumber>
    <ExpirationDate>2023-05-02T00:00:00</ExpirationDate>
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
  <RunID>220913_M02726_0919_KKY8G</RunID>
  <ScannerID>M02726</ScannerID>
  <RunNumber>919</RunNumber>
  <FPGAVersion>9.5.12</FPGAVersion>
  <MCSVersion>2.5.0.5</MCSVersion>
  <RTAVersion>1.18.54</RTAVersion>
  <Barcode>KKY8G</Barcode>
  <PR2BottleBarcode>15041807</PR2BottleBarcode>
  <ReagentKitPartNumberEntered>15033573</ReagentKitPartNumberEntered>
  <ReagentKitVersion>Version2</ReagentKitVersion>
  <ReagentKitBarcode>MS3133997-500V2</ReagentKitBarcode>
  <PreviousPR2BottleBarcode>15041807</PreviousPR2BottleBarcode>
  <PreviousReagentKitBarcode />
  <ExperimentName>Run_Chim</ExperimentName>
  <Chemistry>Amplicon</Chemistry>
  <Username>sbsuser</Username>
  <Workflow>
    <Analysis>GenerateFASTQ</Analysis>
  </Workflow>
  <EnableAnalysis>false</EnableAnalysis>
  <Reads>
    <RunInfoRead Number="1" NumCycles="251" IsIndexedRead="N" />
    <RunInfoRead Number="2" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="3" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="4" NumCycles="251" IsIndexedRead="N" />
  </Reads>
  <TempFolder>D:\Illumina\MiSeqTemp\220913_M02726_0919_KKY8G</TempFolder>
  <AnalysisFolder>D:\Illumina\MiSeqAnalysis\220913_M02726_0919_KKY8G</AnalysisFolder>
  <RunStartDate>220913</RunStartDate>
  <MostRecentWashType>PostRun</MostRecentWashType>
  <RecipeFolder>D:\Illumina\MiSeq Control Software\CustomRecipe</RecipeFolder>
  <ILMNOnlyRecipeFolder>C:\Illumina\MiSeq Control Software\Recipe</ILMNOnlyRecipeFolder>
  <SampleSheetName>Run_Chim_500</SampleSheetName>
  <SampleSheetFolder>D:\lab</SampleSheetFolder>
  <ManifestFolder>D:\Illumina\MiSeq Control Software\Manifests</ManifestFolder>
  <OutputFolder>\\sequencers\MiSeq\220913_M02726_0919_KKY8G</OutputFolder>
  <FocusMethod>AutoFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>false</SaveFocusImages>
  <SaveScanImages>true</SaveScanImages>
  <CloudUsername />
  <RunManagementType>Standalone</RunManagementType>
  <CloudRunId />
  <SendInstrumentHealthToILMN>false</SendInstrumentHealthToILMN>
</RunParameters>"""
        }

    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(tmp_folder, unique_id + "_RunParameters.xml")
        self.test_cases = [
            {  # NexSeq with post-process
                "content": TestRunParameters.CONTENT["NextSeq_withPostProcess"],
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
                "content": TestRunParameters.CONTENT["NextSeq_woutPostProcess"],
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
                "content": TestRunParameters.CONTENT["MiSeq"],
                "expected": {
                    "instrument": {"id": "M01102", "platform": "MiSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 151}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 151}],
                    "run": {"number": "149", "id": "141105_M01102_0149_000000000-AAC1Y", "start_date": datetime(2014, 11, 5)},
                    "kit": {"flowcell_id": "000000000-AAC1Y", "reagent_kit_id": "MS2740035-300V2"},
                    "post_process": "GenerateFASTQ",
                    "software": {"RTA": "1.18.54", "CS": "2.5.0.5"}
                }
            },
            {  # MiSeq invalid RFID
                "content": TestRunParameters.CONTENT["MiSeq_invalid_RFID"],
                "expected": {
                    "instrument": {"id": "M02726", "platform": "MiSeq"},
                    "reads_phases": [{'is_index': False, 'nb_cycles': 251}, {'is_index': True, 'nb_cycles': 8}, {'is_index': True, 'nb_cycles': 8}, {'is_index': False, 'nb_cycles': 251}],
                    "run": {"number": "919", "id": "220913_M02726_0919_KKY8G", "start_date": datetime(2022, 9, 13)},
                    "kit": {"flowcell_id": "KKY8G", "reagent_kit_id": "MS3133997-500V2"},
                    "post_process": "GenerateFASTQ",
                    "software": {"RTA": "1.18.54", "CS": "2.5.0.5"}
                }
            },
            {  # HiSeq
                "content": TestRunParameters.CONTENT["HiSeq"],
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

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

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
