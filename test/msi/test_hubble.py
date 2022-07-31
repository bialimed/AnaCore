#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU-Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
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

from anacore.msi.base import Status
from anacore.msi.hubble import HubbleDiff, HubbleDist, nbRepeatFromStart, parseHubbleResults
from anacore.sequenceIO import IdxFastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
FASTA_CONTENT = """>chr4
GTAGGCTTCCAGTTAGCTTTAAATGTTTTTTTTTTCCAGCTCAAAAAATT
>chr4b
ACATACTATTTACAGCTGCCTTTTTTTTTTTTTTTCGCTA
>chr11
TCGCACTGAGCACATCACATTTAGGACGTTGAAGAtttttttttttttttttttttAATATGCAGTTTGT
>chr11b
gttgttgttgttgttttgtgtttttcggtgtttttttttttaaagaaaaatttacctttg"""

FASTA_INDEX_CONTENT = """chr4	50	6	50	51
chr4b	40	64	40	41
chr11	70	112	70	71
chr11b	60	191	60	61"""


class TestHubbleDiff(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fa = os.path.join(tmp_folder, unique_id + "_ref.fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + "_ref.fa.fai")
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")

        with open(self.tmp_in, "w") as writer:
            writer.write("""#Chromosome	Start	RepeatUnit	Assessed	Distance	PValue
chr4	25	T	True	0.0166435953523248	0.0931187026743588
chr4b	20	T	False	0	0
chr11	35	T	False		
chr11b	30	T	True	0.17974233508488	2.74351359014847E-12
""")
        with open(self.tmp_fa, "w") as writer:
            writer.write(FASTA_CONTENT)
        with open(self.tmp_fai, "w") as writer:
            writer.write(FASTA_INDEX_CONTENT)

    def tearDown(self):
        for curr_file in [self.tmp_fa, self.tmp_fai, self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testInitException(self):
        with self.assertRaises(Exception):
            HubbleDiff(self.tmp_in)

    def testRead(self):
        # Expected
        expected = [
            {
                "position": "chr4:25-35",
                "status": None,
                "distance": 0.0166435953523248,
                "p_value": 0.0931187026743588
            },
            {
                "position": "chr4b:20-35",
                "status": "Undetermined",
                "distance": None,
                "p_value": None
            },
            {
                "position": "chr11:35-56",
                "status": "Undetermined",
                "distance": None,
                "p_value": None
            },
            {
                "position": "chr11b:30-41",
                "status": None,
                "distance": 0.17974233508488,
                "p_value": 2.74351359014847E-12
            }
        ]
        # Observed
        observed = []
        with HubbleDiff(self.tmp_in, ref_path=self.tmp_fa) as reader:
            for record in reader:
                observed.append({
                    "position": record.position,
                    "status": record.results["Hubble"].status,
                    "distance": record.results["Hubble"].data["distance"],
                    "p_value": record.results["Hubble"].data["p_value"]
                })
        # Assert
        self.assertEqual(observed, expected)

    def testWrite(self):
        # Read
        loci = []
        with HubbleDiff(self.tmp_in, ref_path=self.tmp_fa) as reader:
            loci = [record for record in reader]
        # Write
        with HubbleDiff(self.tmp_out, "w") as writer:
            for locus in loci:
                writer.write(locus)
        # Assert
        expected = None
        with open(self.tmp_in) as reader:
            expected = [
                elt.replace("chr4b	20	T	False	0	0", "chr4b	20	T	False		")
                for elt in reader.readlines()
            ]
        observed = None
        with open(self.tmp_out) as reader:
            observed = [
                elt.replace("2.74351359014847e-12", "2.74351359014847E-12")
                for elt in reader.readlines()
            ]
        self.assertEqual(observed, expected)


class TestHubbleDist(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fa = os.path.join(tmp_folder, unique_id + "_ref.fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + "_ref.fa.fai")
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")

        with open(self.tmp_in, "w") as writer:
            writer.write("""#chromosome	location	repeat_unit_bases	covered	length_distribution
chr4	25	T	True	0,0,0,0,0,0,0,0,0,104,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr4b	20	T	False	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr11	35	T	True	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,7,24,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr11b	30	T	True	0,0,0,0,0,0,0,0,0,1,15,124,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        with open(self.tmp_fa, "w") as writer:
            writer.write(FASTA_CONTENT)
        with open(self.tmp_fai, "w") as writer:
            writer.write(FASTA_INDEX_CONTENT)

    def tearDown(self):
        for curr_file in [self.tmp_fa, self.tmp_fai, self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testInitException(self):
        with self.assertRaises(Exception):
            HubbleDist(self.tmp_in)

    def testRead(self):
        # Expected
        expected = [
            {
                "position": "chr4:25-35",
                "status": Status.none,
                "nb_by_length": {10: 104, 11: 1}
            },
            {
                "position": "chr4b:20-35",
                "status": Status.undetermined,
                "nb_by_length": {}
            },
            {
                "position": "chr11:35-56",
                "status": Status.none,
                "nb_by_length": {17: 1, 19: 1, 20: 7, 21: 24, 22: 10}
            },
            {
                "position": "chr11b:30-41",
                "status": Status.none,
                "nb_by_length": {10: 1, 11: 15, 12: 124, 13: 9}
            }
        ]
        # Observed
        observed = []
        with HubbleDist(self.tmp_in, ref_path=self.tmp_fa) as reader:
            for record in reader:
                observed.append({
                    "position": record.position,
                    "status": record.results["Hubble"].status,
                    "nb_by_length": record.results["Hubble"].data["lengths"].ct_by_len
                })
        # Assert
        self.assertEqual(observed, expected)

    def testWrite(self):
        # Read
        loci = []
        with HubbleDist(self.tmp_in, ref_path=self.tmp_fa) as reader:
            loci = [record for record in reader]
        # Write
        with HubbleDist(self.tmp_out, "w") as writer:
            for locus in loci:
                writer.write(locus)
        # Assert
        expected = None
        with open(self.tmp_in) as reader:
            expected = reader.readlines()
        observed = None
        with open(self.tmp_out) as reader:
            observed = reader.readlines()
        self.assertEqual(observed, expected)


class TestNbRepeat(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_fa = os.path.join(tmp_folder, unique_id + "_ref.fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + "_ref.fa.fai")
        with open(self.tmp_fa, "w") as writer:
            writer.write(FASTA_CONTENT)
        with open(self.tmp_fai, "w") as writer:
            writer.write(FASTA_INDEX_CONTENT)

    def tearDown(self):
        for curr_file in [self.tmp_fa, self.tmp_fai]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        data = [
            {"chrom": "chr4", "start": 25, "motif": "T", "expected": 10},
            {"chrom": "chr4", "start": 25, "motif": "TT", "expected": 5},  # Motif with length > 1
            {"chrom": "chr4", "start": 25, "motif": "TTTTT", "expected": 2},  # Motif with length > 1
            {"chrom": "chr4", "start": 48, "motif": "T", "expected": 2},  # End of chrom
            {"chrom": "chr4b", "start": 20, "motif": "T", "expected": 15},
            {"chrom": "chr11", "start": 35, "motif": "T", "expected": 21},  # Repeat lower case
            {"chrom": "chr11", "start": 35, "motif": "TTT", "expected": 7},  # Repeat lower case and motif with length > 1
            {"chrom": "chr11b", "start": 30, "motif": "T", "expected": 11}  # All lower case
        ]
        expected = list()
        observed = list()
        with IdxFastaIO(self.tmp_fa) as ref_reader:
            for curr in data:
                expected.append(curr["expected"])
                observed.append(
                    nbRepeatFromStart(ref_reader, curr["chrom"], curr["start"], curr["motif"])
                )
        self.assertEqual(observed, expected)


class TestParseHubbleResults(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_fa = os.path.join(tmp_folder, unique_id + "_ref.fa")
        self.tmp_fai = os.path.join(tmp_folder, unique_id + "_ref.fa.fai")
        self.tmp_summary = os.path.join(tmp_folder, unique_id + "_sum.json")
        self.tmp_differential = os.path.join(tmp_folder, unique_id + "_diff.tsv")
        self.tmp_distributions = os.path.join(tmp_folder, unique_id + "_dist.tsv")

        with open(self.tmp_summary, "w") as writer:
            writer.write("""{"Settings":{"DetectionOptions":{"QueryBam":"/opt/illumina/analysis-folder/Logs_Intermediates/StitchedRealigned/splA-DNA/splA-DNA.bam","NormalBam":null,"MicrosatelliteFilePath":"/opt/illumina/resources/hubble/References/microsatellites-sites-concordant-500-hg19.list"},"TumorOnlyComparisonOptions":{"PValueThreshold":0.01,"DistanceThreshold":0.1,"NormalReferenceDirectory":"/opt/illumina/resources/hubble/data/normal","SpanningCoverageThreshold":60},"TumorNormalComparisonOptions":null,"AppOptions":{"OutputDirectory":"/opt/illumina/analysis-folder/Logs_Intermediates/Msi/splA-DNA","MaxNumThreads":6,"Debug":false}},"TotalMicrosatelliteSitesAssessed":2,"TotalMicrosatelliteSitesUnstable":1,"PercentageUnstableSites":50.0,"ResultMessage":null,"ResultIsValid":true}""")
        with open(self.tmp_differential, "w") as writer:
            writer.write("""#Chromosome	Start	RepeatUnit	Assessed	Distance	PValue
chr4	25	T	True	0.0166435953523248	0.0931187026743588
chr4b	20	T	False	0	0
chr11	35	T	False		
chr11b	30	T	True	0.17974233508488	2.74351359014847E-12""")
        with open(self.tmp_distributions, "w") as writer:
            writer.write("""#chromosome	location	repeat_unit_bases	covered	length_distribution
chr4	25	T	True	0,0,0,0,0,0,0,0,0,104,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr4b	20	T	False	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr11	35	T	True	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,7,24,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
chr11b	30	T	True	0,0,0,0,0,0,0,0,0,1,15,124,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0""")
        with open(self.tmp_fa, "w") as writer:
            writer.write(FASTA_CONTENT)
        with open(self.tmp_fai, "w") as writer:
            writer.write(FASTA_INDEX_CONTENT)

    def tearDown(self):
        for curr_file in [self.tmp_differential, self.tmp_distributions, self.tmp_fa, self.tmp_fai, self.tmp_summary]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testParser(self):
        msi_spl = parseHubbleResults(
            self.tmp_fa, self.tmp_summary, self.tmp_differential, self.tmp_distributions, 0.2
        )
        # Sample status
        self.assertEqual(msi_spl.results["Hubble"].status, Status.unstable)
        self.assertEqual(msi_spl.results["Hubble"].score, None)
        # Loci status
        self.assertEqual(msi_spl.getNbLoci(), 4)
        self.assertEqual(msi_spl.getNbUnstable("Hubble"), 1)
        self.assertEqual(msi_spl.getNbStable("Hubble"), 1)
        self.assertEqual(msi_spl.getNbDetermined("Hubble"), 2)
        self.assertEqual(msi_spl.getNbUndetermined("Hubble"), 2)
        self.assertEqual(msi_spl.getNbProcessed("Hubble"), 4)
        # Info by locus
        expected = [
            {
                "position": "chr11:35-56",
                "status": Status.undetermined,
                "nb_by_length": {17: 1, 19: 1, 20: 7, 21: 24, 22: 10},
                "nb": 43,
                "distance": None,
                "p_value": None,
                "score": None
            },
            {
                "position": "chr11b:30-41",
                "status": Status.unstable,
                "nb_by_length": {10: 1, 11: 15, 12: 124, 13: 9},
                "nb": 149,
                "distance": 0.17974233508488,
                "p_value": 2.74351359014847E-12,
                "score": 1 - 2.74351359014847E-12,
            },
            {
                "position": "chr4:25-35",
                "status": Status.stable,
                "nb_by_length": {10: 104, 11: 1},
                "nb": 105,
                "distance": 0.0166435953523248,
                "p_value": 0.0931187026743588,
                "score": 1 - 0.0931187026743588
            },
            {
                "position": "chr4b:20-35",
                "status": Status.undetermined,
                "nb_by_length": {},
                "nb": 0,
                "distance": None,
                "p_value": None,
                "score": None
            }
        ]
        observed = []
        for locus_id in sorted(msi_spl.loci):
            locus = msi_spl.loci[locus_id]
            observed.append({
                "position": locus.position,
                "status": locus.results["Hubble"].status,
                "nb_by_length": locus.results["Hubble"].data["lengths"].ct_by_len,
                "nb": locus.results["Hubble"].data["lengths"].getCount(),
                "distance": locus.results["Hubble"].data["distance"],
                "p_value": locus.results["Hubble"].data["p_value"],
                "score": locus.results["Hubble"].score
            })
        self.assertEqual(observed, expected)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
