#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import shutil
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.illumina.demultiplex import DemultStatFactory
from anacore.illumina.demultiplex.bcl2fastq import DemultStat as DemultStatBcl2fastq
from anacore.illumina.demultiplex.bclconvert import DemultStat as DemultStatBclConvert


class TestDemultStatFactory(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in_folder = os.path.join(tmp_folder, unique_id)

    def testGetBcl2fastq(self):
        os.makedirs(self.tmp_in_folder)
        os.makedirs(os.path.join(self.tmp_in_folder, "Stats"))
        with open(os.path.join(self.tmp_in_folder, "Stats", "Stats.json"), "w") as writer:
            writer.write('{"Flowcell":"HFF5HBGXL","RunNumber":256,"RunId":"220620_NDX550421_RUO_0256_AHFF5HBGXL","ReadInfosForLanes":[{"LaneNumber":1,"ReadInfos":[{"Number":1,"NumCycles":101,"IsIndexedRead":false},{"Number":1,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":101,"IsIndexedRead":false}]},{"LaneNumber":2,"ReadInfos":[{"Number":1,"NumCycles":101,"IsIndexedRead":false},{"Number":1,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":8,"IsIndexedRead":true},{"Number":2,"NumCycles":101,"IsIndexedRead":false}]}],"ConversionResults":[{"LaneNumber":1,"TotalClustersRaw":118154091,"TotalClustersPF":110173898,"Yield":22255127396,"DemuxResults":[{"SampleId":"splA-DNA","SampleName":"splA-DNA","IndexMetrics":[{"IndexSequence":"TTAATCAG+CTTCGCCT","MismatchCounts":{"0":9478792,"1":141994}}],"NumberReads":9620786,"Yield":1943398772,"ReadMetrics":[{"ReadNumber":1,"Yield":971699386,"YieldQ30":924128846,"QualityScoreSum":33738537475,"TrimmedBases":48420782},{"ReadNumber":2,"Yield":971699386,"YieldQ30":918562030,"QualityScoreSum":33615137854,"TrimmedBases":47944431}]},{"SampleId":"splA-RNA","SampleName":"splA-RNA","IndexMetrics":[{"IndexSequence":"TCCGGAGA+AGGATAGG","MismatchCounts":{"0":2699794,"1":40900}}],"NumberReads":2740694,"Yield":553620188,"ReadMetrics":[{"ReadNumber":1,"Yield":276810094,"YieldQ30":261295087,"QualityScoreSum":9569976358,"TrimmedBases":34397480},{"ReadNumber":2,"Yield":276810094,"YieldQ30":237377116,"QualityScoreSum":9065991028,"TrimmedBases":31656823}]},{"SampleId":"splB-DNA","SampleName":"splB-DNA","IndexMetrics":[{"IndexSequence":"CGCTCATT+TAAGATTA","MismatchCounts":{"0":9564216,"1":131182}}],"NumberReads":9695398,"Yield":1958470396,"ReadMetrics":[{"ReadNumber":1,"Yield":979235198,"YieldQ30":935923618,"QualityScoreSum":34097986656,"TrimmedBases":37858895},{"ReadNumber":2,"Yield":979235198,"YieldQ30":929589854,"QualityScoreSum":33958376843,"TrimmedBases":37570166}]},{"SampleId":"splB-RNA","SampleName":"splB-RNA","IndexMetrics":[{"IndexSequence":"CTGAAGCT+TCAGAGCC","MismatchCounts":{"0":2597339,"1":49487}}],"NumberReads":2646826,"Yield":534658852,"ReadMetrics":[{"ReadNumber":1,"Yield":267329426,"YieldQ30":254469868,"QualityScoreSum":9286867469,"TrimmedBases":8407441},{"ReadNumber":2,"Yield":267329426,"YieldQ30":250039144,"QualityScoreSum":9191210821,"TrimmedBases":6605292}]},{"SampleId":"splC-DNA","SampleName":"splC-DNA","IndexMetrics":[{"IndexSequence":"TCCGCGAA+AGTAAGTA","MismatchCounts":{"0":9152620,"1":121190}}],"NumberReads":9273810,"Yield":1873309620,"ReadMetrics":[{"ReadNumber":1,"Yield":936654810,"YieldQ30":892413065,"QualityScoreSum":32555394461,"TrimmedBases":31828023},{"ReadNumber":2,"Yield":936654810,"YieldQ30":885028001,"QualityScoreSum":32393432900,"TrimmedBases":31576068}]},{"SampleId":"splD-DNA","SampleName":"splD-DNA","IndexMetrics":[{"IndexSequence":"ATTACTCG+GACTTCCT","MismatchCounts":{"0":9922693,"1":138354}}],"NumberReads":10061047,"Yield":2032331494,"ReadMetrics":[{"ReadNumber":1,"Yield":1016165747,"YieldQ30":966795535,"QualityScoreSum":35290830751,"TrimmedBases":53576473},{"ReadNumber":2,"Yield":1016165747,"YieldQ30":959254789,"QualityScoreSum":35125522244,"TrimmedBases":53052819}]}],"Undetermined":{"NumberReads":10126886,"Yield":2045630972,"ReadMetrics":[{"ReadNumber":1,"Yield":1022815486,"YieldQ30":942308528,"QualityScoreSum":34857988295,"TrimmedBases":50185140},{"ReadNumber":2,"Yield":1022815486,"YieldQ30":919589580,"QualityScoreSum":34250317568,"TrimmedBases":47568306}]}},{"LaneNumber":2,"TotalClustersRaw":116151639,"TotalClustersPF":108452034,"Yield":21907310868,"DemuxResults":[{"SampleId":"splA-DNA","SampleName":"splA-DNA","IndexMetrics":[{"IndexSequence":"TTAATCAG+CTTCGCCT","MismatchCounts":{"0":9312034,"1":144771}}],"NumberReads":9456805,"Yield":1910274610,"ReadMetrics":[{"ReadNumber":1,"Yield":955137305,"YieldQ30":908955716,"QualityScoreSum":33173997089,"TrimmedBases":47382071},{"ReadNumber":2,"Yield":955137305,"YieldQ30":903901935,"QualityScoreSum":33061153007,"TrimmedBases":46935857}]},{"SampleId":"splA-RNA","SampleName":"splA-RNA","IndexMetrics":[{"IndexSequence":"TCCGGAGA+AGGATAGG","MismatchCounts":{"0":2653660,"1":41311}}],"NumberReads":2694971,"Yield":544384142,"ReadMetrics":[{"ReadNumber":1,"Yield":272192071,"YieldQ30":257193370,"QualityScoreSum":9415122745,"TrimmedBases":33684765},{"ReadNumber":2,"Yield":272192071,"YieldQ30":233971488,"QualityScoreSum":8925076315,"TrimmedBases":31059435}]},{"SampleId":"splB-DNA","SampleName":"splB-DNA","IndexMetrics":[{"IndexSequence":"CGCTCATT+TAAGATTA","MismatchCounts":{"0":9389668,"1":133733}}],"NumberReads":9523401,"Yield":1923727002,"ReadMetrics":[{"ReadNumber":1,"Yield":961863501,"YieldQ30":919992361,"QualityScoreSum":33505686273,"TrimmedBases":37143142},{"ReadNumber":2,"Yield":961863501,"YieldQ30":914031955,"QualityScoreSum":33373758590,"TrimmedBases":36872757}]},{"SampleId":"splB-RNA","SampleName":"splB-RNA","IndexMetrics":[{"IndexSequence":"CTGAAGCT+TCAGAGCC","MismatchCounts":{"0":2552513,"1":49897}}],"NumberReads":2602410,"Yield":525686820,"ReadMetrics":[{"ReadNumber":1,"Yield":262843410,"YieldQ30":250401718,"QualityScoreSum":9134727206,"TrimmedBases":8237441},{"ReadNumber":2,"Yield":262843410,"YieldQ30":246175047,"QualityScoreSum":9043525469,"TrimmedBases":6480551}]},{"SampleId":"splC-DNA","SampleName":"splC-DNA","IndexMetrics":[{"IndexSequence":"TCCGCGAA+AGTAAGTA","MismatchCounts":{"0":8997000,"1":121956}}],"NumberReads":9118956,"Yield":1842029112,"ReadMetrics":[{"ReadNumber":1,"Yield":921014556,"YieldQ30":878282672,"QualityScoreSum":32026485795,"TrimmedBases":31172049},{"ReadNumber":2,"Yield":921014556,"YieldQ30":871422340,"QualityScoreSum":31875534840,"TrimmedBases":30932132}]},{"SampleId":"splD-DNA","SampleName":"splD-DNA","IndexMetrics":[{"IndexSequence":"ATTACTCG+GACTTCCT","MismatchCounts":{"0":9755783,"1":141064}}],"NumberReads":9896847,"Yield":1999163094,"ReadMetrics":[{"ReadNumber":1,"Yield":999581547,"YieldQ30":951751473,"QualityScoreSum":34728541382,"TrimmedBases":52608236},{"ReadNumber":2,"Yield":999581547,"YieldQ30":944728842,"QualityScoreSum":34573985187,"TrimmedBases":52120720}]}],"Undetermined":{"NumberReads":10144184,"Yield":2049125168,"ReadMetrics":[{"ReadNumber":1,"Yield":1024562584,"YieldQ30":944586733,"QualityScoreSum":34933559142,"TrimmedBases":50457930},{"ReadNumber":2,"Yield":1024562584,"YieldQ30":922237756,"QualityScoreSum":34343749503,"TrimmedBases":47618888}]}}],"UnknownBarcodes":[{"Lane":1,"Barcodes":{"GGGGGGGG+AGATCTCG":3267780,"ACTGCTTA+AGATCTCG":304460,"ATGCGGCT+AGATCTCG":280440,"AAAAAAAT+TAGCCGCG":234920,"AAAAAAAT+TTCGTAGG":223900,"AAAAAAAT+AGTAAGTA":221040}},{"Lane":2,"Barcodes":{"GGGGGGGG+AGATCTCG":3216540,"ACTGCTTA+AGATCTCG":292100,"ATGCGGCT+AGATCTCG":276940,"AAAAAAAT+TAGCCGCG":228260,"AAAAAAAT+TTCGTAGG":220120,"AAAAAAAT+AGTAAGTA":216460}}]}')
        stat = DemultStatFactory.get(self.tmp_in_folder)
        self.assertEqual(stat.__class__, DemultStatBcl2fastq)

    def testGetBclConvert(self):
        os.makedirs(self.tmp_in_folder)
        os.makedirs(os.path.join(self.tmp_in_folder, "Reports"))
        with open(os.path.join(self.tmp_in_folder, "Reports", "Demultiplex_Stats.csv"), "w") as writer:
            writer.write("""Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads
1,splA,CTGATCGT-GCGCATAT,1546203,1511172,35031,0,0.135021,0.176742,0.181226,0.000000
1,splB,ACTCTCGA-CTGTACCA,2585489,2522498,62991,0,0.225775,0.295023,0.325872,0.000000
2,splA,CGCATGAT-AAGCCTGA,2293141,2246343,46798,0,0.200246,0.262725,0.242100,0.000000
2,splB,ACGGAACA-ACGAGAAC,2318638,2270158,48480,0,0.202473,0.265510,0.250802,0.000000""")
        with open(os.path.join(self.tmp_in_folder, "Reports", "Top_Unknown_Barcodes.csv"), "w") as writer:
            writer.write("""Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads
1,GGGGGGGG,AGATCTCG,2610292,0.963872,0.227941
1,TGGACTCT,TCTTACGG,37763,0.013944,0.003298
2,TTGCATTC,AAGGCGTA,31132,0.011496,0.002719
2,TGGACTCT,TCTTACGG,28943,0.010687,0.002527""")
        stat = DemultStatFactory.get(self.tmp_in_folder)
        self.assertEqual(stat.__class__, DemultStatBclConvert)

    def testGetException(self):
        os.makedirs(self.tmp_in_folder)
        os.makedirs(os.path.join(self.tmp_in_folder, "Other"))
        with self.assertRaises(IOError):
            DemultStatFactory.get(self.tmp_in_folder)

    def tearDown(self):
        if os.path.exists(self.tmp_in_folder):
            shutil.rmtree(self.tmp_in_folder)


if __name__ == "__main__":
    unittest.main()
