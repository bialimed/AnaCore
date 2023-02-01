#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import datetime
import json
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina.demultiplex.bcl2fastq import DemultLog, DemultStat


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
                    "version": "v2.20.0.422",
                    "start_time": datetime.datetime(2021, 6, 3, 17, 18, 10),
                    "status": "COMPLETED",
                    "end_time": datetime.datetime(2021, 6, 3, 23, 22, 13)
                }
            }, {  # Incomplete
                "content": """BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.

2021-06-03 17:18:10 [2afc55848400] Command-line invocation: bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2
2021-06-03 17:18:10 [2afc55848400] INFO: Minimum log level: INFO""",
                "expected": {
                    "command": "bcl2fastq --no-lane-splitting --loading-threads 1 --writing-threads 1 --processing-threads 16 --runfolder-dir /NextSeq3/210602_NDX550500_RUO_0063_AHVW3TAFX2",
                    "version": "v2.20.0.422",
                    "start_time": datetime.datetime(2021, 6, 3, 17, 18, 10),
                    "status": "RUNNING",
                    "end_time": None
                }
            }, {  # Error
                "content": """BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.

2023-01-31 11:30:53 [7feae79637c0] Command-line invocation: bcl2fastq --no-lane-splitting --loading-threads 2 --writing-threads 2 --processing-threads 28 --runfolder-dir /seq/230110_NDX550421_RUO_0377_AH5JKKAFX5
2023-01-31 11:30:53 [7feae79637c0] INFO: Minimum log level: INFO
2023-01-31 11:30:53 [7feae79637c0] INFO: Sample sheet: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/SampleSheet.csv'
2023-01-31 11:30:53 [7feae79637c0] INFO: Runfolder path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5'
2023-01-31 11:30:53 [7feae79637c0] INFO: Input path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/Data/Intensities/BaseCalls/'
2023-01-31 11:30:53 [7feae79637c0] INFO: Intensities path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/Data/Intensities/'
2023-01-31 11:30:53 [7feae79637c0] INFO: Output path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/Data/Intensities/BaseCalls/'
2023-01-31 11:30:53 [7feae79637c0] INFO: InterOp path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/InterOp/'
2023-01-31 11:30:53 [7feae79637c0] INFO: Stats path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/Data/Intensities/BaseCalls/Stats/'
2023-01-31 11:30:53 [7feae79637c0] INFO: Reports path: '/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/Data/Intensities/BaseCalls/Reports/'
2023-01-31 11:30:53 [7feae79637c0] INFO: Detected CPUs: 128
2023-01-31 11:30:53 [7feae79637c0] INFO: Loading threads: 2
2023-01-31 11:30:53 [7feae79637c0] INFO: Processing threads: 28
2023-01-31 11:30:53 [7feae79637c0] INFO: Writing threads: 2
2023-01-31 11:30:53 [7feae79637c0] INFO: Allowed barcode mismatches: 1
2023-01-31 11:30:53 [7feae79637c0] INFO: Tiles: <ALL>
2023-01-31 11:30:53 [7feae79637c0] INFO: Minimum trimmed read length: 35
2023-01-31 11:30:53 [7feae79637c0] INFO: Use bases masks: <NONE>
2023-01-31 11:30:53 [7feae79637c0] INFO: Mask short adapter reads: 22
2023-01-31 11:30:53 [7feae79637c0] INFO: Adapter stringency: 0.9
2023-01-31 11:30:53 [7feae79637c0] INFO: Adapter trimming method: Allow matches with indels
2023-01-31 11:30:53 [7feae79637c0] INFO: Ignore missing BCLs: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Ignore missing filters: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Ignore missing positions: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Ignore missing controls: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Include non-PF clusters: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Create FASTQs for index reads: NO
2023-01-31 11:30:53 [7feae79637c0] INFO: Use bgzf compression for FASTQ files: YES
2023-01-31 11:30:53 [7feae79637c0] INFO: FASTQ compression level: 4
2023-01-31 11:30:53 [7feae79637c0] INFO: RunInfo.xml: '"/seq/230110_NDX550421_RUO_0377_AH5JKKAFX5/RunInfo.xml"'
2023-01-31 11:30:53 [7feae79637c0] INFO: Lane: 1
2023-01-31 11:30:53 [7feae79637c0] INFO: Mask:
2023-01-31 11:30:53 [7feae79637c0] INFO:   Sample: #0 'unknown' 'Undetermined' [default]
2023-01-31 11:30:53 [7feae79637c0] INFO:   Sample: #1 'PBMC32-40ng' 'PBMC32-40ng' [default]
2023-01-31 11:30:53 [7feae79637c0] INFO:     Barcode: 'CTGATCGT'
2023-01-31 11:30:53 [7feae79637c0] INFO:   Sample: #2 'KARPASrep1-40ng' 'KARPASrep1-40ng' [default]
2023-01-31 11:30:53 [7feae79637c0] INFO:     Barcode: 'ACTCTCGA+CTGTACCA'
2023-01-31 11:30:53 [7feae79637c0] ERROR: bcl2fastq::common::Exception: 2023-Jan-31 11:30:53: Success (0): /soft/bcl2fastq/v2.20.0.422/src/cxx/lib/layout/BarcodeCollisionDetector.cpp(108): Throw in function void bcl2fastq::layout::BarcodeCollisionDetector::validateNewBarcodeSizesAgainstExisting(const std::vector<long unsigned int>&) const
Dynamic exception type: boost::exception_detail::clone_impl<bcl2fastq::layout::BarcodeCollisionError>
std::exception::what: Barcodes have an unequal number of components.""",
                "expected": {
                    "command": "bcl2fastq --no-lane-splitting --loading-threads 2 --writing-threads 2 --processing-threads 28 --runfolder-dir /seq/230110_NDX550421_RUO_0377_AH5JKKAFX5",
                    "version": "v2.20.0.422",
                    "start_time": datetime.datetime(2023, 1, 31, 11, 30, 53),
                    "status": "FAILED",
                    "end_time": datetime.datetime(2023, 1, 31, 11, 30, 53)
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
            res = DemultLog(self.tmp_file)
            observed.append({
                "command": res.command,
                "version": res.version,
                "start_time": res.start_time,
                "status": res.status,
                "end_time": res.end_time
            })
        self.assertEqual(expected, observed)


class TestDemultStatBcl2fastq(unittest.TestCase):
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
        # Increase count of undetermined
        with open(self.tmp_file) as reader:
            data = json.load(reader)
            for lane in data["UnknownBarcodes"]:
                for barcode, ct in lane["Barcodes"].items():
                    lane["Barcodes"][barcode] = ct * 10
        with open(self.tmp_file, "w") as writer:
            json.dump(data, writer)
        stats = DemultStat(self.tmp_file)
        # RNA under
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


if __name__ == "__main__":
    unittest.main()
