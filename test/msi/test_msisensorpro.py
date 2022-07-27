#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU-Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.msi.base import Status
from anacore.msi.locus import Locus, LocusDataDistrib
from anacore.msi.msisensorpro import BaselineIO, BaselineRecord, binaryToString, DistIO, parseProResults, ProEval, ProIO, stringToBinary
from anacore.msi.sample import MSISample
from anacore.region import Region
from anacore.sequenceIO import IdxFastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestBaselineIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")

        with open(self.tmp_in, "w") as writer:
            writer.write("""chromosome	location	repeat_unit_length	repeat_unit_binary	repeat_times	left_flank_binary 	right_flank_binary	repeat_unit_bases	left_flank_bases	right_flank_bases	threshold	supportSamples
2	47414420	1	0	27	299	687	A	CAGGT	GGGTT	0.203254	93
2	95183613	1	3	23	860	738	T	TCCTA	GTGAG	0.1071	72
""")

        self.records = [
            BaselineRecord("2", 47414420, 27, "A", "CAGGT", "GGGTT", 0.203254, 93),
            BaselineRecord("2", 95183613, 23, "T", "TCCTA", "GTGAG", 0.1071, 72),
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testRead(self):
        with BaselineIO(self.tmp_in) as reader:
            self.assertEqual(
                [rec.__dict__ for rec in reader],
                [rec.__dict__ for rec in self.records]
            )

    def testWrite(self):
        with BaselineIO(self.tmp_out, "w") as writer:
            for rec in self.records:
                writer.write(rec)
        with open(self.tmp_in) as reader:
            expected = reader.readlines()
        expected[0] = expected[0].replace("left_flank_binary ", "left_flank_binary")  # Bug in MSIsensor file
        with open(self.tmp_out) as reader:
            observed = reader.readlines()
        self.assertEqual(observed, expected)


class TestBaselineRecord(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_fasta = os.path.join(tmp_folder, unique_id + ".fasta")
        with open(self.tmp_fasta, "w") as FH_out:
            FH_out.write(">chr1\nATGGCATAGCATAATCATTTTTTTTTTAGTTTCACA")
        self.tmp_fasta_idx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        with open(self.tmp_fasta_idx, "w") as FH_out:
            FH_out.write("chr1	36	6	100	101")
        self.record = BaselineRecord("chr1", 17, 10, "T", "ATCA", "AGTT", 0.13, 3)

    def tearDown(self):
        for curr_file in [self.tmp_fasta, self.tmp_fasta_idx]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testFromModel(self):
        microsat = Region(18, 27, None, "chr1")
        models = [
            MSISample.fromDict({
                "name": "splA",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "score": 1,
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.09
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splB",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.1
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splC",
                "results": {
                    "train": {
                        "status": "MSI",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSI",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.2
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splD",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.11
                                    }
                                }
                            }
                        }
                    }
                }
            })
        ]
        with IdxFastaIO(self.tmp_fasta) as ref_fh:
            observed = BaselineRecord.fromModel(ref_fh, microsat, models, "train", 4)
            self.assertEqual(
                observed.__dict__,
                self.record.__dict__
            )


class TestBinary(unittest.TestCase):
    def setUp(self):
        self.data = [
            {"bin": 0, "str": "A"},
            {"bin": 1003, "str": "TTGGT"},
            {"bin": 1004, "str": "TTGTA"},
            {"bin": 1005, "str": "TTGTC"},
            {"bin": 1006, "str": "TTGTG"},
            {"bin": 1007, "str": "TTGTT"},
            {"bin": 1008, "str": "TTTAA"},
            {"bin": 1009, "str": "TTTAC"},
            {"bin": 100, "str": "CGCA"},
            {"bin": 1010, "str": "TTTAG"},
            {"bin": 1011, "str": "TTTAT"},
            {"bin": 1012, "str": "TTTCA"},
            {"bin": 1013, "str": "TTTCC"},
            {"bin": 1015, "str": "TTTCT"},
            {"bin": 1016, "str": "TTTGA"},
            {"bin": 1017, "str": "TTTGC"},
            {"bin": 1018, "str": "TTTGG"},
            {"bin": 1019, "str": "TTTGT"},
            {"bin": 101, "str": "CGCC"},
            {"bin": 1020, "str": "TTTTA"},
            {"bin": 1021, "str": "TTTTC"},
            {"bin": 1022, "str": "TTTTG"},
            {"bin": 103, "str": "CGCT"},
            {"bin": 106, "str": "CGGG"},
            {"bin": 108, "str": "ACGTA"},
            {"bin": 109, "str": "CGTC"},
            {"bin": 10, "str": "AAAGG"},
            {"bin": 10, "str": "AAGG"},
            {"bin": 10, "str": "AGG"},
            {"bin": 110, "str": "ACGTG"},
            {"bin": 110, "str": "CGTG"},
            {"bin": 111, "str": "ACGTT"},
            {"bin": 112, "str": "ACTAA"},
            {"bin": 112, "str": "CTAA"},
            {"bin": 113, "str": "ACTAC"},
            {"bin": 113, "str": "CTAC"},
            {"bin": 114, "str": "ACTAG"},
            {"bin": 114, "str": "CTAG"},
            {"bin": 115, "str": "ACTAT"},
            {"bin": 115, "str": "CTAT"},
            {"bin": 116, "str": "CTCA"},
            {"bin": 117, "str": "ACTCC"},
            {"bin": 117, "str": "CTCC"},
            {"bin": 118, "str": "CTCG"},
            {"bin": 119, "str": "ACTCT"},
            {"bin": 11, "str": "AAAGT"},
            {"bin": 11, "str": "AAGT"},
            {"bin": 11, "str": "AGT"},
            {"bin": 11, "str": "GT"},
            {"bin": 120, "str": "CTGA"},
            {"bin": 121, "str": "CTGC"},
            {"bin": 122, "str": "CTGG"},
            {"bin": 123, "str": "CTGT"},
            {"bin": 124, "str": "CTTA"},
            {"bin": 125, "str": "CTTC"},
            {"bin": 126, "str": "CTTG"},
            {"bin": 127, "str": "CTTT"},
            {"bin": 128, "str": "AGAAA"},
            {"bin": 128, "str": "GAAA"},
            {"bin": 129, "str": "AGAAC"},
            {"bin": 129, "str": "GAAC"},
            {"bin": 12, "str": "AAATA"},
            {"bin": 12, "str": "AATA"},
            {"bin": 12, "str": "ATA"},
            {"bin": 12, "str": "TA"},
            {"bin": 130, "str": "AGAAG"},
            {"bin": 130, "str": "GAAG"},
            {"bin": 131, "str": "AGAAT"},
            {"bin": 131, "str": "GAAT"},
            {"bin": 132, "str": "GACA"},
            {"bin": 134, "str": "GACG"},
            {"bin": 135, "str": "GACT"},
            {"bin": 136, "str": "AGAGA"},
            {"bin": 137, "str": "GAGC"},
            {"bin": 138, "str": "AGAGG"},
            {"bin": 138, "str": "GAGG"},
            {"bin": 139, "str": "GAGT"},
            {"bin": 13, "str": "AAATC"},
            {"bin": 13, "str": "AATC"},
            {"bin": 13, "str": "ATC"},
            {"bin": 13, "str": "TC"},
            {"bin": 140, "str": "AGATA"},
            {"bin": 140, "str": "GATA"},
            {"bin": 141, "str": "GATC"},
            {"bin": 142, "str": "AGATG"},
            {"bin": 142, "str": "GATG"},
            {"bin": 143, "str": "GATT"},
            {"bin": 144, "str": "GCAA"},
            {"bin": 145, "str": "AGCAC"},
            {"bin": 145, "str": "GCAC"},
            {"bin": 146, "str": "AGCAG"},
            {"bin": 146, "str": "GCAG"},
            {"bin": 147, "str": "AGCAT"},
            {"bin": 147, "str": "GCAT"},
            {"bin": 148, "str": "GCCA"},
            {"bin": 149, "str": "AGCCC"},
            {"bin": 149, "str": "GCCC"},
            {"bin": 14, "str": "AAATG"},
            {"bin": 14, "str": "AATG"},
            {"bin": 14, "str": "ATG"},
            {"bin": 14, "str": "TG"},
            {"bin": 150, "str": "GCCG"},
            {"bin": 151, "str": "GCCT"},
            {"bin": 154, "str": "GCGG"},
            {"bin": 155, "str": "GCGT"},
            {"bin": 156, "str": "GCTA"},
            {"bin": 157, "str": "GCTC"},
            {"bin": 158, "str": "GCTG"},
            {"bin": 159, "str": "GCTT"},
            {"bin": 15, "str": "AAATT"},
            {"bin": 15, "str": "AATT"},
            {"bin": 15, "str": "ATT"},
            {"bin": 160, "str": "AGGAA"},
            {"bin": 160, "str": "GGAA"},
            {"bin": 161, "str": "AGGAC"},
            {"bin": 161, "str": "GGAC"},
            {"bin": 162, "str": "AGGAG"},
            {"bin": 162, "str": "GGAG"},
            {"bin": 163, "str": "GGAT"},
            {"bin": 164, "str": "GGCA"},
            {"bin": 166, "str": "GGCG"},
            {"bin": 167, "str": "GGCT"},
            {"bin": 168, "str": "AGGGA"},
            {"bin": 168, "str": "GGGA"},
            {"bin": 169, "str": "AGGGC"},
            {"bin": 169, "str": "GGGC"},
            {"bin": 16, "str": "AACAA"},
            {"bin": 16, "str": "ACAA"},
            {"bin": 16, "str": "CAA"},
            {"bin": 170, "str": "AGGGG"},
            {"bin": 171, "str": "GGGT"},
            {"bin": 172, "str": "GGTA"},
            {"bin": 173, "str": "GGTC"},
            {"bin": 174, "str": "AGGTG"},
            {"bin": 174, "str": "GGTG"},
            {"bin": 175, "str": "GGTT"},
            {"bin": 176, "str": "GTAA"},
            {"bin": 177, "str": "AGTAC"},
            {"bin": 177, "str": "GTAC"},
            {"bin": 178, "str": "AGTAG"},
            {"bin": 178, "str": "GTAG"},
            {"bin": 179, "str": "AGTAT"},
            {"bin": 179, "str": "GTAT"},
            {"bin": 17, "str": "AACAC"},
            {"bin": 17, "str": "CAC"},
            {"bin": 180, "str": "GTCA"},
            {"bin": 181, "str": "GTCC"},
            {"bin": 183, "str": "GTCT"},
            {"bin": 184, "str": "AGTGA"},
            {"bin": 184, "str": "GTGA"},
            {"bin": 185, "str": "GTGC"},
            {"bin": 186, "str": "AGTGG"},
            {"bin": 186, "str": "GTGG"},
            {"bin": 187, "str": "AGTGT"},
            {"bin": 188, "str": "GTTA"},
            {"bin": 189, "str": "GTTC"},
            {"bin": 18, "str": "AACAG"},
            {"bin": 18, "str": "ACAG"},
            {"bin": 18, "str": "CAG"},
            {"bin": 190, "str": "GTTG"},
            {"bin": 191, "str": "AGTTT"},
            {"bin": 191, "str": "GTTT"},
            {"bin": 192, "str": "ATAAA"},
            {"bin": 192, "str": "TAAA"},
            {"bin": 193, "str": "ATAAC"},
            {"bin": 193, "str": "TAAC"},
            {"bin": 194, "str": "ATAAG"},
            {"bin": 194, "str": "TAAG"},
            {"bin": 195, "str": "ATAAT"},
            {"bin": 195, "str": "TAAT"},
            {"bin": 196, "str": "ATACA"},
            {"bin": 196, "str": "TACA"},
            {"bin": 197, "str": "ATACC"},
            {"bin": 197, "str": "TACC"},
            {"bin": 198, "str": "ATACG"},
            {"bin": 199, "str": "ATACT"},
            {"bin": 199, "str": "TACT"},
            {"bin": 19, "str": "AACAT"},
            {"bin": 19, "str": "ACAT"},
            {"bin": 19, "str": "CAT"},
            {"bin": 1, "str": "AAAAC"},
            {"bin": 1, "str": "AAAC"},
            {"bin": 1, "str": "AAC"},
            {"bin": 1, "str": "AC"},
            {"bin": 1, "str": "C"},
            {"bin": 200, "str": "ATAGA"},
            {"bin": 200, "str": "TAGA"},
            {"bin": 201, "str": "ATAGC"},
            {"bin": 201, "str": "TAGC"},
            {"bin": 202, "str": "ATAGG"},
            {"bin": 202, "str": "TAGG"},
            {"bin": 203, "str": "ATAGT"},
            {"bin": 203, "str": "TAGT"},
            {"bin": 204, "str": "ATATA"},
            {"bin": 205, "str": "TATC"},
            {"bin": 206, "str": "ATATG"},
            {"bin": 206, "str": "TATG"},
            {"bin": 207, "str": "ATATT"},
            {"bin": 207, "str": "TATT"},
            {"bin": 208, "str": "ATCAA"},
            {"bin": 208, "str": "TCAA"},
            {"bin": 209, "str": "ATCAC"},
            {"bin": 209, "str": "TCAC"},
            {"bin": 20, "str": "AACCA"},
            {"bin": 20, "str": "ACCA"},
            {"bin": 20, "str": "CCA"},
            {"bin": 210, "str": "TCAG"},
            {"bin": 211, "str": "ATCAT"},
            {"bin": 211, "str": "TCAT"},
            {"bin": 212, "str": "TCCA"},
            {"bin": 213, "str": "ATCCC"},
            {"bin": 213, "str": "TCCC"},
            {"bin": 214, "str": "TCCG"},
            {"bin": 215, "str": "ATCCT"},
            {"bin": 215, "str": "TCCT"},
            {"bin": 216, "str": "TCGA"},
            {"bin": 217, "str": "TCGC"},
            {"bin": 21, "str": "AACCC"},
            {"bin": 21, "str": "ACCC"},
            {"bin": 220, "str": "ATCTA"},
            {"bin": 220, "str": "TCTA"},
            {"bin": 221, "str": "ATCTC"},
            {"bin": 222, "str": "TCTG"},
            {"bin": 223, "str": "ATCTT"},
            {"bin": 223, "str": "TCTT"},
            {"bin": 224, "str": "ATGAA"},
            {"bin": 224, "str": "TGAA"},
            {"bin": 225, "str": "ATGAC"},
            {"bin": 225, "str": "TGAC"},
            {"bin": 226, "str": "TGAG"},
            {"bin": 227, "str": "ATGAT"},
            {"bin": 227, "str": "TGAT"},
            {"bin": 228, "str": "TGCA"},
            {"bin": 229, "str": "ATGCC"},
            {"bin": 229, "str": "TGCC"},
            {"bin": 22, "str": "ACCG"},
            {"bin": 22, "str": "CCG"},
            {"bin": 230, "str": "TGCG"},
            {"bin": 231, "str": "ATGCT"},
            {"bin": 231, "str": "TGCT"},
            {"bin": 232, "str": "ATGGA"},
            {"bin": 232, "str": "TGGA"},
            {"bin": 233, "str": "TGGC"},
            {"bin": 234, "str": "ATGGG"},
            {"bin": 234, "str": "TGGG"},
            {"bin": 235, "str": "ATGGT"},
            {"bin": 235, "str": "TGGT"},
            {"bin": 236, "str": "ATGTA"},
            {"bin": 236, "str": "TGTA"},
            {"bin": 237, "str": "ATGTC"},
            {"bin": 237, "str": "TGTC"},
            {"bin": 238, "str": "ATGTG"},
            {"bin": 239, "str": "ATGTT"},
            {"bin": 239, "str": "TGTT"},
            {"bin": 23, "str": "ACCT"},
            {"bin": 23, "str": "CCT"},
            {"bin": 240, "str": "ATTAA"},
            {"bin": 240, "str": "TTAA"},
            {"bin": 241, "str": "ATTAC"},
            {"bin": 241, "str": "TTAC"},
            {"bin": 242, "str": "ATTAG"},
            {"bin": 242, "str": "TTAG"},
            {"bin": 243, "str": "ATTAT"},
            {"bin": 243, "str": "TTAT"},
            {"bin": 244, "str": "ATTCA"},
            {"bin": 244, "str": "TTCA"},
            {"bin": 245, "str": "ATTCC"},
            {"bin": 245, "str": "TTCC"},
            {"bin": 246, "str": "TTCG"},
            {"bin": 247, "str": "ATTCT"},
            {"bin": 247, "str": "TTCT"},
            {"bin": 248, "str": "ATTGA"},
            {"bin": 248, "str": "TTGA"},
            {"bin": 249, "str": "ATTGC"},
            {"bin": 249, "str": "TTGC"},
            {"bin": 24, "str": "CGA"},
            {"bin": 250, "str": "TTGG"},
            {"bin": 251, "str": "ATTGT"},
            {"bin": 251, "str": "TTGT"},
            {"bin": 252, "str": "ATTTA"},
            {"bin": 252, "str": "TTTA"},
            {"bin": 253, "str": "ATTTC"},
            {"bin": 253, "str": "TTTC"},
            {"bin": 254, "str": "ATTTG"},
            {"bin": 254, "str": "TTTG"},
            {"bin": 255, "str": "ATTTT"},
            {"bin": 256, "str": "CAAAA"},
            {"bin": 257, "str": "CAAAC"},
            {"bin": 258, "str": "CAAAG"},
            {"bin": 25, "str": "ACGC"},
            {"bin": 25, "str": "CGC"},
            {"bin": 260, "str": "CAACA"},
            {"bin": 266, "str": "CAAGG"},
            {"bin": 268, "str": "CAATA"},
            {"bin": 26, "str": "AACGG"},
            {"bin": 26, "str": "CGG"},
            {"bin": 270, "str": "CAATG"},
            {"bin": 271, "str": "CAATT"},
            {"bin": 272, "str": "CACAA"},
            {"bin": 273, "str": "CACAC"},
            {"bin": 274, "str": "CACAG"},
            {"bin": 277, "str": "CACCC"},
            {"bin": 278, "str": "CACCG"},
            {"bin": 27, "str": "CGT"},
            {"bin": 284, "str": "CACTA"},
            {"bin": 285, "str": "CACTC"},
            {"bin": 286, "str": "CACTG"},
            {"bin": 288, "str": "CAGAA"},
            {"bin": 289, "str": "CAGAC"},
            {"bin": 28, "str": "AACTA"},
            {"bin": 28, "str": "ACTA"},
            {"bin": 28, "str": "CTA"},
            {"bin": 292, "str": "CAGCA"},
            {"bin": 293, "str": "CAGCC"},
            {"bin": 295, "str": "CAGCT"},
            {"bin": 297, "str": "CAGGC"},
            {"bin": 298, "str": "CAGGG"},
            {"bin": 29, "str": "ACTC"},
            {"bin": 29, "str": "CTC"},
            {"bin": 2, "str": "AAAAG"},
            {"bin": 2, "str": "AAAG"},
            {"bin": 2, "str": "AAG"},
            {"bin": 2, "str": "AG"},
            {"bin": 2, "str": "G"},
            {"bin": 300, "str": "CAGTA"},
            {"bin": 302, "str": "CAGTG"},
            {"bin": 304, "str": "CATAA"},
            {"bin": 306, "str": "CATAG"},
            {"bin": 309, "str": "CATCC"},
            {"bin": 30, "str": "AACTG"},
            {"bin": 30, "str": "ACTG"},
            {"bin": 30, "str": "CTG"},
            {"bin": 311, "str": "CATCT"},
            {"bin": 315, "str": "CATGT"},
            {"bin": 316, "str": "CATTA"},
            {"bin": 317, "str": "CATTC"},
            {"bin": 319, "str": "CATTT"},
            {"bin": 31, "str": "ACTT"},
            {"bin": 31, "str": "CTT"},
            {"bin": 320, "str": "CCAAA"},
            {"bin": 323, "str": "CCAAT"},
            {"bin": 324, "str": "CCACA"},
            {"bin": 325, "str": "CCACC"},
            {"bin": 327, "str": "CCACT"},
            {"bin": 329, "str": "CCAGC"},
            {"bin": 32, "str": "AAGAA"},
            {"bin": 32, "str": "AGAA"},
            {"bin": 32, "str": "GAA"},
            {"bin": 332, "str": "CCATA"},
            {"bin": 333, "str": "CCATC"},
            {"bin": 335, "str": "CCATT"},
            {"bin": 337, "str": "CCCAC"},
            {"bin": 338, "str": "CCCAG"},
            {"bin": 339, "str": "CCCAT"},
            {"bin": 33, "str": "AAGAC"},
            {"bin": 33, "str": "AGAC"},
            {"bin": 33, "str": "GAC"},
            {"bin": 340, "str": "CCCCA"},
            {"bin": 342, "str": "CCCCG"},
            {"bin": 343, "str": "CCCCT"},
            {"bin": 345, "str": "CCCGC"},
            {"bin": 347, "str": "CCCGT"},
            {"bin": 348, "str": "CCCTA"},
            {"bin": 349, "str": "CCCTC"},
            {"bin": 34, "str": "AAGAG"},
            {"bin": 34, "str": "GAG"},
            {"bin": 350, "str": "CCCTG"},
            {"bin": 351, "str": "CCCTT"},
            {"bin": 354, "str": "CCGAG"},
            {"bin": 357, "str": "CCGCC"},
            {"bin": 358, "str": "CCGCG"},
            {"bin": 35, "str": "AAGAT"},
            {"bin": 35, "str": "AGAT"},
            {"bin": 35, "str": "GAT"},
            {"bin": 369, "str": "CCTAC"},
            {"bin": 36, "str": "AAGCA"},
            {"bin": 36, "str": "AGCA"},
            {"bin": 36, "str": "GCA"},
            {"bin": 371, "str": "CCTAT"},
            {"bin": 372, "str": "CCTCA"},
            {"bin": 373, "str": "CCTCC"},
            {"bin": 375, "str": "CCTCT"},
            {"bin": 377, "str": "CCTGC"},
            {"bin": 379, "str": "CCTGT"},
            {"bin": 37, "str": "AGCC"},
            {"bin": 37, "str": "GCC"},
            {"bin": 380, "str": "CCTTA"},
            {"bin": 381, "str": "CCTTC"},
            {"bin": 382, "str": "CCTTG"},
            {"bin": 383, "str": "CCTTT"},
            {"bin": 387, "str": "CGAAT"},
            {"bin": 38, "str": "AGCG"},
            {"bin": 38, "str": "GCG"},
            {"bin": 393, "str": "CGAGC"},
            {"bin": 394, "str": "CGAGG"},
            {"bin": 39, "str": "AGCT"},
            {"bin": 39, "str": "GCT"},
            {"bin": 3, "str": "AAAAT"},
            {"bin": 3, "str": "AAAT"},
            {"bin": 3, "str": "AAT"},
            {"bin": 3, "str": "AT"},
            {"bin": 3, "str": "T"},
            {"bin": 405, "str": "CGCCC"},
            {"bin": 406, "str": "CGCCG"},
            {"bin": 407, "str": "CGCCT"},
            {"bin": 40, "str": "AAGGA"},
            {"bin": 40, "str": "AGGA"},
            {"bin": 40, "str": "GGA"},
            {"bin": 414, "str": "CGCTG"},
            {"bin": 418, "str": "CGGAG"},
            {"bin": 41, "str": "AAGGC"},
            {"bin": 41, "str": "AGGC"},
            {"bin": 41, "str": "GGC"},
            {"bin": 423, "str": "CGGCT"},
            {"bin": 426, "str": "CGGGG"},
            {"bin": 42, "str": "AAGGG"},
            {"bin": 42, "str": "AGGG"},
            {"bin": 432, "str": "CGTAA"},
            {"bin": 437, "str": "CGTCC"},
            {"bin": 43, "str": "AAGGT"},
            {"bin": 43, "str": "AGGT"},
            {"bin": 43, "str": "GGT"},
            {"bin": 448, "str": "CTAAA"},
            {"bin": 449, "str": "CTAAC"},
            {"bin": 44, "str": "AAGTA"},
            {"bin": 44, "str": "AGTA"},
            {"bin": 44, "str": "GTA"},
            {"bin": 453, "str": "CTACC"},
            {"bin": 455, "str": "CTACT"},
            {"bin": 45, "str": "AGTC"},
            {"bin": 45, "str": "GTC"},
            {"bin": 460, "str": "CTATA"},
            {"bin": 461, "str": "CTATC"},
            {"bin": 462, "str": "CTATG"},
            {"bin": 463, "str": "CTATT"},
            {"bin": 464, "str": "CTCAA"},
            {"bin": 467, "str": "CTCAT"},
            {"bin": 468, "str": "CTCCA"},
            {"bin": 469, "str": "CTCCC"},
            {"bin": 46, "str": "AAGTG"},
            {"bin": 46, "str": "AGTG"},
            {"bin": 46, "str": "GTG"},
            {"bin": 471, "str": "CTCCT"},
            {"bin": 474, "str": "CTCGG"},
            {"bin": 477, "str": "CTCTC"},
            {"bin": 478, "str": "CTCTG"},
            {"bin": 479, "str": "CTCTT"},
            {"bin": 47, "str": "AGTT"},
            {"bin": 47, "str": "GTT"},
            {"bin": 485, "str": "CTGCC"},
            {"bin": 487, "str": "CTGCT"},
            {"bin": 489, "str": "CTGGC"},
            {"bin": 48, "str": "AATAA"},
            {"bin": 48, "str": "ATAA"},
            {"bin": 48, "str": "TAA"},
            {"bin": 490, "str": "CTGGG"},
            {"bin": 492, "str": "CTGTA"},
            {"bin": 493, "str": "CTGTC"},
            {"bin": 494, "str": "CTGTG"},
            {"bin": 495, "str": "CTGTT"},
            {"bin": 49, "str": "AATAC"},
            {"bin": 49, "str": "ATAC"},
            {"bin": 49, "str": "TAC"},
            {"bin": 4, "str": "AAACA"},
            {"bin": 4, "str": "AACA"},
            {"bin": 4, "str": "ACA"},
            {"bin": 4, "str": "CA"},
            {"bin": 501, "str": "CTTCC"},
            {"bin": 503, "str": "CTTCT"},
            {"bin": 506, "str": "CTTGG"},
            {"bin": 509, "str": "CTTTC"},
            {"bin": 50, "str": "AATAG"},
            {"bin": 50, "str": "ATAG"},
            {"bin": 50, "str": "TAG"},
            {"bin": 511, "str": "CTTTT"},
            {"bin": 512, "str": "GAAAA"},
            {"bin": 514, "str": "GAAAG"},
            {"bin": 515, "str": "GAAAT"},
            {"bin": 516, "str": "GAACA"},
            {"bin": 517, "str": "GAACC"},
            {"bin": 518, "str": "GAACG"},
            {"bin": 51, "str": "AATAT"},
            {"bin": 51, "str": "TAT"},
            {"bin": 520, "str": "GAAGA"},
            {"bin": 522, "str": "GAAGG"},
            {"bin": 524, "str": "GAATA"},
            {"bin": 526, "str": "GAATG"},
            {"bin": 52, "str": "AATCA"},
            {"bin": 52, "str": "ATCA"},
            {"bin": 52, "str": "TCA"},
            {"bin": 530, "str": "GACAG"},
            {"bin": 538, "str": "GACGG"},
            {"bin": 53, "str": "AATCC"},
            {"bin": 53, "str": "ATCC"},
            {"bin": 53, "str": "TCC"},
            {"bin": 540, "str": "GACTA"},
            {"bin": 544, "str": "GAGAA"},
            {"bin": 546, "str": "GAGAG"},
            {"bin": 547, "str": "GAGAT"},
            {"bin": 54, "str": "TCG"},
            {"bin": 552, "str": "GAGGA"},
            {"bin": 553, "str": "GAGGC"},
            {"bin": 554, "str": "GAGGG"},
            {"bin": 555, "str": "GAGGT"},
            {"bin": 556, "str": "GAGTA"},
            {"bin": 558, "str": "GAGTG"},
            {"bin": 55, "str": "AATCT"},
            {"bin": 55, "str": "ATCT"},
            {"bin": 55, "str": "TCT"},
            {"bin": 560, "str": "GATAA"},
            {"bin": 561, "str": "GATAC"},
            {"bin": 563, "str": "GATAT"},
            {"bin": 568, "str": "GATGA"},
            {"bin": 569, "str": "GATGC"},
            {"bin": 56, "str": "AATGA"},
            {"bin": 56, "str": "ATGA"},
            {"bin": 56, "str": "TGA"},
            {"bin": 570, "str": "GATGG"},
            {"bin": 572, "str": "GATTA"},
            {"bin": 574, "str": "GATTG"},
            {"bin": 575, "str": "GATTT"},
            {"bin": 576, "str": "GCAAA"},
            {"bin": 577, "str": "GCAAC"},
            {"bin": 579, "str": "GCAAT"},
            {"bin": 57, "str": "AATGC"},
            {"bin": 57, "str": "ATGC"},
            {"bin": 57, "str": "TGC"},
            {"bin": 586, "str": "GCAGG"},
            {"bin": 587, "str": "GCAGT"},
            {"bin": 588, "str": "GCATA"},
            {"bin": 589, "str": "GCATC"},
            {"bin": 58, "str": "AATGG"},
            {"bin": 58, "str": "ATGG"},
            {"bin": 58, "str": "TGG"},
            {"bin": 594, "str": "GCCAG"},
            {"bin": 595, "str": "GCCAT"},
            {"bin": 596, "str": "GCCCA"},
            {"bin": 597, "str": "GCCCC"},
            {"bin": 598, "str": "GCCCG"},
            {"bin": 599, "str": "GCCCT"},
            {"bin": 59, "str": "AATGT"},
            {"bin": 59, "str": "ATGT"},
            {"bin": 59, "str": "TGT"},
            {"bin": 5, "str": "AAACC"},
            {"bin": 5, "str": "AACC"},
            {"bin": 5, "str": "ACC"},
            {"bin": 602, "str": "GCCGG"},
            {"bin": 607, "str": "GCCTT"},
            {"bin": 60, "str": "AATTA"},
            {"bin": 60, "str": "ATTA"},
            {"bin": 60, "str": "TTA"},
            {"bin": 616, "str": "GCGGA"},
            {"bin": 617, "str": "GCGGC"},
            {"bin": 618, "str": "GCGGG"},
            {"bin": 61, "str": "AATTC"},
            {"bin": 61, "str": "ATTC"},
            {"bin": 61, "str": "TTC"},
            {"bin": 622, "str": "GCGTG"},
            {"bin": 627, "str": "GCTAT"},
            {"bin": 628, "str": "GCTCA"},
            {"bin": 62, "str": "AATTG"},
            {"bin": 62, "str": "ATTG"},
            {"bin": 62, "str": "TTG"},
            {"bin": 630, "str": "GCTCG"},
            {"bin": 631, "str": "GCTCT"},
            {"bin": 634, "str": "GCTGG"},
            {"bin": 635, "str": "GCTGT"},
            {"bin": 639, "str": "GCTTT"},
            {"bin": 63, "str": "AATTT"},
            {"bin": 63, "str": "ATTT"},
            {"bin": 640, "str": "GGAAA"},
            {"bin": 641, "str": "GGAAC"},
            {"bin": 642, "str": "GGAAG"},
            {"bin": 643, "str": "GGAAT"},
            {"bin": 646, "str": "GGACG"},
            {"bin": 647, "str": "GGACT"},
            {"bin": 648, "str": "GGAGA"},
            {"bin": 64, "str": "ACAAA"},
            {"bin": 64, "str": "CAAA"},
            {"bin": 650, "str": "GGAGG"},
            {"bin": 651, "str": "GGAGT"},
            {"bin": 654, "str": "GGATG"},
            {"bin": 656, "str": "GGCAA"},
            {"bin": 658, "str": "GGCAG"},
            {"bin": 65, "str": "ACAAC"},
            {"bin": 65, "str": "CAAC"},
            {"bin": 661, "str": "GGCCC"},
            {"bin": 662, "str": "GGCCG"},
            {"bin": 665, "str": "GGCGC"},
            {"bin": 666, "str": "GGCGG"},
            {"bin": 667, "str": "GGCGT"},
            {"bin": 668, "str": "GGCTA"},
            {"bin": 66, "str": "ACAAG"},
            {"bin": 66, "str": "CAAG"},
            {"bin": 670, "str": "GGCTG"},
            {"bin": 672, "str": "GGGAA"},
            {"bin": 674, "str": "GGGAG"},
            {"bin": 675, "str": "GGGAT"},
            {"bin": 676, "str": "GGGCA"},
            {"bin": 678, "str": "GGGCG"},
            {"bin": 679, "str": "GGGCT"},
            {"bin": 67, "str": "ACAAT"},
            {"bin": 67, "str": "CAAT"},
            {"bin": 680, "str": "GGGGA"},
            {"bin": 681, "str": "GGGGC"},
            {"bin": 683, "str": "GGGGT"},
            {"bin": 686, "str": "GGGTG"},
            {"bin": 68, "str": "ACACA"},
            {"bin": 695, "str": "GGTCT"},
            {"bin": 698, "str": "GGTGG"},
            {"bin": 69, "str": "CACC"},
            {"bin": 6, "str": "AAACG"},
            {"bin": 6, "str": "ACG"},
            {"bin": 6, "str": "CG"},
            {"bin": 700, "str": "GGTTA"},
            {"bin": 701, "str": "GGTTC"},
            {"bin": 703, "str": "GGTTT"},
            {"bin": 705, "str": "GTAAC"},
            {"bin": 709, "str": "GTACC"},
            {"bin": 70, "str": "CACG"},
            {"bin": 714, "str": "GTAGG"},
            {"bin": 716, "str": "GTATA"},
            {"bin": 717, "str": "GTATC"},
            {"bin": 719, "str": "GTATT"},
            {"bin": 71, "str": "ACACT"},
            {"bin": 71, "str": "CACT"},
            {"bin": 727, "str": "GTCCT"},
            {"bin": 72, "str": "CAGA"},
            {"bin": 733, "str": "GTCTC"},
            {"bin": 735, "str": "GTCTT"},
            {"bin": 737, "str": "GTGAC"},
            {"bin": 73, "str": "ACAGC"},
            {"bin": 73, "str": "CAGC"},
            {"bin": 743, "str": "GTGCT"},
            {"bin": 744, "str": "GTGGA"},
            {"bin": 746, "str": "GTGGG"},
            {"bin": 74, "str": "ACAGG"},
            {"bin": 74, "str": "CAGG"},
            {"bin": 750, "str": "GTGTG"},
            {"bin": 751, "str": "GTGTT"},
            {"bin": 753, "str": "GTTAC"},
            {"bin": 755, "str": "GTTAT"},
            {"bin": 757, "str": "GTTCC"},
            {"bin": 75, "str": "CAGT"},
            {"bin": 764, "str": "GTTTA"},
            {"bin": 766, "str": "GTTTG"},
            {"bin": 767, "str": "GTTTT"},
            {"bin": 768, "str": "TAAAA"},
            {"bin": 769, "str": "TAAAC"},
            {"bin": 76, "str": "ACATA"},
            {"bin": 76, "str": "CATA"},
            {"bin": 770, "str": "TAAAG"},
            {"bin": 771, "str": "TAAAT"},
            {"bin": 772, "str": "TAACA"},
            {"bin": 77, "str": "ACATC"},
            {"bin": 77, "str": "CATC"},
            {"bin": 780, "str": "TAATA"},
            {"bin": 781, "str": "TAATC"},
            {"bin": 782, "str": "TAATG"},
            {"bin": 783, "str": "TAATT"},
            {"bin": 784, "str": "TACAA"},
            {"bin": 785, "str": "TACAC"},
            {"bin": 788, "str": "TACCA"},
            {"bin": 789, "str": "TACCC"},
            {"bin": 78, "str": "ACATG"},
            {"bin": 78, "str": "CATG"},
            {"bin": 796, "str": "TACTA"},
            {"bin": 797, "str": "TACTC"},
            {"bin": 798, "str": "TACTG"},
            {"bin": 79, "str": "ACATT"},
            {"bin": 79, "str": "CATT"},
            {"bin": 7, "str": "AAACT"},
            {"bin": 7, "str": "AACT"},
            {"bin": 7, "str": "ACT"},
            {"bin": 7, "str": "CT"},
            {"bin": 800, "str": "TAGAA"},
            {"bin": 801, "str": "TAGAC"},
            {"bin": 802, "str": "TAGAG"},
            {"bin": 803, "str": "TAGAT"},
            {"bin": 804, "str": "TAGCA"},
            {"bin": 809, "str": "TAGGC"},
            {"bin": 80, "str": "ACCAA"},
            {"bin": 80, "str": "CCAA"},
            {"bin": 810, "str": "TAGGG"},
            {"bin": 812, "str": "TAGTA"},
            {"bin": 813, "str": "TAGTC"},
            {"bin": 815, "str": "TAGTT"},
            {"bin": 816, "str": "TATAA"},
            {"bin": 817, "str": "TATAC"},
            {"bin": 818, "str": "TATAG"},
            {"bin": 819, "str": "TATAT"},
            {"bin": 81, "str": "CCAC"},
            {"bin": 820, "str": "TATCA"},
            {"bin": 821, "str": "TATCC"},
            {"bin": 823, "str": "TATCT"},
            {"bin": 824, "str": "TATGA"},
            {"bin": 825, "str": "TATGC"},
            {"bin": 826, "str": "TATGG"},
            {"bin": 827, "str": "TATGT"},
            {"bin": 828, "str": "TATTA"},
            {"bin": 829, "str": "TATTC"},
            {"bin": 82, "str": "ACCAG"},
            {"bin": 82, "str": "CCAG"},
            {"bin": 830, "str": "TATTG"},
            {"bin": 831, "str": "TATTT"},
            {"bin": 835, "str": "TCAAT"},
            {"bin": 836, "str": "TCACA"},
            {"bin": 83, "str": "ACCAT"},
            {"bin": 83, "str": "CCAT"},
            {"bin": 841, "str": "TCAGC"},
            {"bin": 844, "str": "TCATA"},
            {"bin": 845, "str": "TCATC"},
            {"bin": 846, "str": "TCATG"},
            {"bin": 847, "str": "TCATT"},
            {"bin": 849, "str": "TCCAC"},
            {"bin": 84, "str": "CCCA"},
            {"bin": 851, "str": "TCCAT"},
            {"bin": 852, "str": "TCCCA"},
            {"bin": 853, "str": "TCCCC"},
            {"bin": 854, "str": "TCCCG"},
            {"bin": 855, "str": "TCCCT"},
            {"bin": 859, "str": "TCCGT"},
            {"bin": 85, "str": "ACCCC"},
            {"bin": 860, "str": "TCCTA"},
            {"bin": 861, "str": "TCCTC"},
            {"bin": 862, "str": "TCCTG"},
            {"bin": 863, "str": "TCCTT"},
            {"bin": 86, "str": "CCCG"},
            {"bin": 87, "str": "ACCCT"},
            {"bin": 87, "str": "CCCT"},
            {"bin": 882, "str": "TCTAG"},
            {"bin": 883, "str": "TCTAT"},
            {"bin": 884, "str": "TCTCA"},
            {"bin": 885, "str": "TCTCC"},
            {"bin": 886, "str": "TCTCG"},
            {"bin": 887, "str": "TCTCT"},
            {"bin": 890, "str": "TCTGG"},
            {"bin": 891, "str": "TCTGT"},
            {"bin": 892, "str": "TCTTA"},
            {"bin": 893, "str": "TCTTC"},
            {"bin": 894, "str": "TCTTG"},
            {"bin": 895, "str": "TCTTT"},
            {"bin": 896, "str": "TGAAA"},
            {"bin": 89, "str": "CCGC"},
            {"bin": 8, "str": "AAAGA"},
            {"bin": 8, "str": "AAGA"},
            {"bin": 8, "str": "AGA"},
            {"bin": 8, "str": "GA"},
            {"bin": 904, "str": "TGAGA"},
            {"bin": 905, "str": "TGAGC"},
            {"bin": 908, "str": "TGATA"},
            {"bin": 913, "str": "TGCAC"},
            {"bin": 914, "str": "TGCAG"},
            {"bin": 915, "str": "TGCAT"},
            {"bin": 917, "str": "TGCCC"},
            {"bin": 919, "str": "TGCCT"},
            {"bin": 91, "str": "CCGT"},
            {"bin": 924, "str": "TGCTA"},
            {"bin": 925, "str": "TGCTC"},
            {"bin": 926, "str": "TGCTG"},
            {"bin": 928, "str": "TGGAA"},
            {"bin": 929, "str": "TGGAC"},
            {"bin": 92, "str": "CCTA"},
            {"bin": 930, "str": "TGGAG"},
            {"bin": 931, "str": "TGGAT"},
            {"bin": 933, "str": "TGGCC"},
            {"bin": 936, "str": "TGGGA"},
            {"bin": 937, "str": "TGGGC"},
            {"bin": 938, "str": "TGGGG"},
            {"bin": 93, "str": "CCTC"},
            {"bin": 941, "str": "TGGTC"},
            {"bin": 942, "str": "TGGTG"},
            {"bin": 943, "str": "TGGTT"},
            {"bin": 944, "str": "TGTAA"},
            {"bin": 945, "str": "TGTAC"},
            {"bin": 946, "str": "TGTAG"},
            {"bin": 947, "str": "TGTAT"},
            {"bin": 949, "str": "TGTCC"},
            {"bin": 94, "str": "CCTG"},
            {"bin": 951, "str": "TGTCT"},
            {"bin": 952, "str": "TGTGA"},
            {"bin": 953, "str": "TGTGC"},
            {"bin": 954, "str": "TGTGG"},
            {"bin": 955, "str": "TGTGT"},
            {"bin": 956, "str": "TGTTA"},
            {"bin": 957, "str": "TGTTC"},
            {"bin": 958, "str": "TGTTG"},
            {"bin": 959, "str": "TGTTT"},
            {"bin": 95, "str": "ACCTT"},
            {"bin": 95, "str": "CCTT"},
            {"bin": 960, "str": "TTAAA"},
            {"bin": 963, "str": "TTAAT"},
            {"bin": 964, "str": "TTACA"},
            {"bin": 965, "str": "TTACC"},
            {"bin": 966, "str": "TTACG"},
            {"bin": 968, "str": "TTAGA"},
            {"bin": 96, "str": "CGAA"},
            {"bin": 970, "str": "TTAGG"},
            {"bin": 971, "str": "TTAGT"},
            {"bin": 972, "str": "TTATA"},
            {"bin": 973, "str": "TTATC"},
            {"bin": 974, "str": "TTATG"},
            {"bin": 975, "str": "TTATT"},
            {"bin": 976, "str": "TTCAA"},
            {"bin": 979, "str": "TTCAT"},
            {"bin": 980, "str": "TTCCA"},
            {"bin": 981, "str": "TTCCC"},
            {"bin": 982, "str": "TTCCG"},
            {"bin": 983, "str": "TTCCT"},
            {"bin": 988, "str": "TTCTA"},
            {"bin": 989, "str": "TTCTC"},
            {"bin": 98, "str": "CGAG"},
            {"bin": 990, "str": "TTCTG"},
            {"bin": 991, "str": "TTCTT"},
            {"bin": 992, "str": "TTGAA"},
            {"bin": 996, "str": "TTGCA"},
            {"bin": 997, "str": "TTGCC"},
            {"bin": 999, "str": "TTGCT"},
            {"bin": 9, "str": "AAAGC"},
            {"bin": 9, "str": "AAGC"},
            {"bin": 9, "str": "AGC"},
            {"bin": 9, "str": "GC"}
        ]

    def testBinaryToString(self):
        expected = [elt["str"] for elt in self.data]
        observed = [binaryToString(elt["bin"], len(elt["str"])) for elt in self.data]
        self.assertEqual(observed, expected)

    def testStringToBinary(self):
        expected = [elt["bin"] for elt in self.data]
        observed = [stringToBinary(elt["str"]) for elt in self.data]
        self.assertEqual(observed, expected)


class TestDistIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.txt")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.txt")

        with open(self.tmp_in, "w") as writer:
            writer.write("""2 47414420 CAGGT 27[A] GGGTT
T: 0 0 0 0 0 0 0 0 0 0 0 0 0 7 2 2 9 15 16 22 29 63 67 62 61 33 12 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
2 95183613 TCCTA 23[T] GTGAG
T: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 2 6 8 11 23 17 20 7 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
""")
        self.records = [
            (
                Locus.fromDict({
                    "position": "2:47414420",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": Status.none,
                            "data": {
                                "lengths": {
                                    "ct_by_len": {
                                        14: 7,
                                        15: 2,
                                        16: 2,
                                        17: 9,
                                        18: 15,
                                        19: 16,
                                        20: 22,
                                        21: 29,
                                        22: 63,
                                        23: 67,
                                        24: 62,
                                        25: 61,
                                        26: 33,
                                        27: 12,
                                        28: 2,
                                        29: 1,
                                        30: 1
                                    }
                                }
                            }
                        }
                    }
                }),
                {
                    "chromosome": "2",
                    "location": 47414420,
                    "left_flank_bases": "CAGGT",
                    "repeat_times": 27,
                    "repeat_unit_bases": "A",
                    "right_flank_bases": "GGGTT"
                }
            ), (
                Locus.fromDict({
                    "position": "2:95183613",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": Status.none,
                            "data": {
                                "lengths": {
                                    "ct_by_len": {
                                        16: 3,
                                        17: 2,
                                        18: 6,
                                        19: 8,
                                        20: 11,
                                        21: 23,
                                        22: 17,
                                        23: 20,
                                        24: 7,
                                        25: 1,
                                        26: 1
                                    }
                                }
                            }
                        }
                    }
                }),
                {
                    "chromosome": "2",
                    "location": 95183613,
                    "left_flank_bases": "TCCTA",
                    "repeat_times": 23,
                    "repeat_unit_bases": "T",
                    "right_flank_bases": "GTGAG"
                }
            )
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testRead(self):
        with DistIO(self.tmp_in) as reader:
            self.assertEqual(
                [(toDict(locus), meta) for locus, meta in reader],
                [(toDict(locus), meta) for locus, meta in self.records]
            )

    def testWrite(self):
        with DistIO(self.tmp_out, "w") as writer:
            for locus, meta in self.records:
                writer.write(locus, meta)
        with open(self.tmp_in) as reader:
            expected = [line.replace(" \n", "\n") for line in reader.readlines()]  # Line ends with ' \n'
        with open(self.tmp_out) as reader:
            observed = reader.readlines()
        self.assertEqual(observed, expected)


class TestParseMSIsensorProResults(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_all = os.path.join(tmp_folder, unique_id + "_all.tsv")
        self.tmp_distrib = os.path.join(tmp_folder, unique_id + "_dis.txt")
        with open(self.tmp_all, "w") as writer:
            writer.write("""chromosome	location	left_flank_bases	repeat_times	repeat_unit_bases	right_flank_bases	pro_p	pro_q	CovReads	threshold
11	125620870	GAAGA	21	T	AATAT	0.094124	0.002322	1064	0.082020
13	31148483	CACTT	17	A	TCATG	0.088336	0.000664	4070	0.118326""")
        with open(self.tmp_distrib, "w") as writer:
            writer.write("""11 125620870 GAAGA 21[T] AATAT
T: 0 0 0 0 0 0 0 0 0 0 0 0 3 8 18 40 113 179 238 255 168 32 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
13 31148483 CACTT 17[A] TCATG
T: 0 0 0 0 0 0 0 0 0 1 2 21 84 423 1210 1967 320 38 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 """)

    def tearDown(self):
        for curr_file in [self.tmp_all, self.tmp_distrib]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        expected = MSISample.fromDict({
            "name": "splA",
            "results": {
                "MSIsensor-pro_pro": {
                    "status": "MSI",
                    "score": 0.5,
                    "method": "MSIsensor-pro_pro",
                    "param": {
                        "aggregation_method": "instability ratio",
                        "instability_threshold": 0.1,
                        "min_voting_loci": 0.5

                    }
                }
            },
            "loci": {
                "11:125620870": {
                    "position": "11:125620870",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": "MSI",
                            "data": {
                                "pro_p": 0.094124,
                                "pro_q": 0.002322,
                                "depth": 1064,
                                "lengths": {
                                    "ct_by_len": {
                                        13: 3,
                                        14: 8,
                                        15: 18,
                                        16: 40,
                                        17: 113,
                                        18: 179,
                                        19: 238,
                                        20: 255,
                                        21: 168,
                                        22: 32,
                                        23: 10
                                    }
                                }
                            }
                        }
                    }
                },
                "13:31148483": {
                    "position": "13:31148483",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": "MSS",
                            "data": {
                                "pro_p": 0.088336,
                                "pro_q": 0.000664,
                                "depth": 4070,
                                "lengths": {
                                    "ct_by_len": {
                                        10: 1,
                                        11: 2,
                                        12: 21,
                                        13: 84,
                                        14: 423,
                                        15: 1210,
                                        16: 1967,
                                        17: 320,
                                        18: 38,
                                        19: 4
                                    }
                                }
                            }
                        }
                    }
                }
            }
        })
        observed = parseProResults(
            "splA",
            self.tmp_all,
            self.tmp_distrib,
            unstable_rate_threshold=0.1,
            determined_rate_threshold=0.5,
            min_support=20
        )
        self.assertEqual(toDict(observed), toDict(expected))


class TestProEval(unittest.TestCase):
    def testGetSlippageScores(self):
        data = [
            {
                "ref_len": 27,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 7 2 2 9 15 16 22 29 63 67 62 61 33 12 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.159597, 0.000641)
            },
            {
                "ref_len": 23,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 2 6 8 11 23 17 20 7 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.083443, 0.005242)
            },
            {
                "ref_len": 25,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 1 8 0 10 7 18 16 23 35 23 11 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.099361, 0.003831)
            },
            {
                "ref_len": 26,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 4 9 9 19 26 42 49 47 54 21 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.121951, 0.001049)
            },
            {
                "ref_len": 21,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 11 28 67 100 150 152 101 23 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.054563, 0.011692)
            },
            {
                "ref_len": 17,
                "distrib_str": "0 0 0 0 0 0 0 0 0 1 2 21 84 423 1210 1967 320 38 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.088336, 0.000664)
            },
            {
                "ref_len": 21,
                "distrib_str": "0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 1 3 7 11 33 50 68 74 43 12 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "expected": (0.015962, 0.060828)
            }
        ]
        for curr in data:
            locus_distrib = LocusDataDistrib.fromDense(
                [int(elt) for elt in curr["distrib_str"].split(" ")]
            )
            self.assertEqual(
                tuple(round(elt, 6) for elt in ProEval.getSlippageScores(locus_distrib, curr["ref_len"])),
                curr["expected"]
            )

    def testGetThreshold(self):
        models = [
            MSISample.fromDict({
                "name": "splA",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "score": 1,
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.09
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splB",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.1
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splC",
                "results": {
                    "train": {
                        "status": "MSI",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSI",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.2
                                    }
                                }
                            }
                        }
                    }
                }
            }),
            MSISample.fromDict({
                "name": "splD",
                "results": {
                    "train": {
                        "status": "MSS",
                        "score": 1
                    }
                },
                "loci": {
                    "chr1:17": {
                        "position": "chr1:17",
                        "results": {
                            "train": {
                                "status": "MSS",
                                "data": {
                                    "MSIsensor-pro": {
                                        "pro_p": 0.11
                                    }
                                }
                            }
                        }
                    }
                }
            })
        ]
        self.assertEqual(
            ProEval.getThreshold(models, "chr1:17", "train"),
            0.13
        )

    def testGetThresholdFromScores(self):
        data = [
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 3 6 9 23 22 32 45 46 49 23 4 3 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 2 0 3 0 2 12 27 39 50 81 79 90 76 23 6 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 2 0 0 0 1 8 11 10 28 49 72 91 128 174 147 72 21 3 0 0 0 0 0 0 1 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 1 6 9 23 33 54 68 110 107 80 45 10 2 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 4 15 33 49 56 88 132 132 130 60 24 4 2 0 1 0 0 0 0 0 0"
                ],
                "ref_len": 27,
                "expected": 0.171801
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 1 3 5 9 17 20 10 3 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 1 1 7 12 21 42 61 47 20 3 1 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 3 9 18 31 51 70 83 59 18 4 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 5 9 14 24 49 56 64 37 7 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5 8 21 29 57 95 103 36 12 4 2 0 0 0 0 0 0 0"
                ],
                "ref_len": 23,
                "expected": 0.076841
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 6 8 23 18 36 42 56 18 6 0 1 0 0 0 0 0 0 0 0 1 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 4 5 21 10 23 53 69 65 29 10 0 1 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 2 2 5 1 2 2 13 16 20 40 54 100 123 50 31 8 3 1 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 7 22 27 29 64 70 122 87 39 8 2 1 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 3 4 11 11 32 37 52 62 35 9 1 0 0 0 0 0 0 0 0 0 0 0 0"
                ],
                "ref_len": 25,
                "expected": 0.0870342
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 1 8 14 32 31 45 37 47 11 3 1 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 4 9 20 25 32 54 58 70 48 19 4 0 1 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 7 15 32 51 69 86 104 90 59 15 1 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 3 3 14 23 45 46 51 44 19 7 0 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 4 5 8 23 28 52 61 79 56 41 9 4 0 0 0 0 0 0 0 0 0 0 0 0 0"
                ],
                "ref_len": 26,
                "expected": 0.166947
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 6 2 9 44 75 173 154 102 17 7 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 13 26 84 173 245 275 161 26 6 1 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 8 42 110 232 337 408 232 36 4 1 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 18 49 139 229 196 32 1 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 3 4 19 56 128 213 325 251 73 8 2 0 0 0 0 0 0 0 0"
                ],
                "ref_len": 21,
                "expected": 0.0786864
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 1 4 35 136 360 900 463 65 8 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 7 21 66 267 682 1284 832 138 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 1 1 27 145 465 1307 2662 1433 214 26 2 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 1 1 1 45 159 448 926 1406 180 18 5 0 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 1 10 41 152 609 1294 1840 305 44 9 1 0 0 0 0 0 0 0 0 0 0 0"
                ],
                "ref_len": 17,
                "expected": 0.0964866
            },
            {
                "distrib": [
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 5 7 11 37 51 66 87 30 4 0 0 0 0 0 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 3 10 23 33 64 110 151 130 45 11 3 1 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 2 0 6 2 19 27 64 127 200 278 278 105 23 3 0 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 7 10 18 41 57 103 142 138 76 27 9 1 1 0 0 0 0 0 0",
                    "0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 1 6 17 32 44 71 121 131 76 25 4 1 0 1 0 0 0 0 0 0"
                ],
                "ref_len": 21,
                "expected": 0.0231436
            }
        ]
        for curr in data:
            pro_p = []
            for distrib in curr["distrib"]:
                locus_distrib = LocusDataDistrib.fromDense([int(elt) for elt in distrib.split(" ")])
                pro_p.append(
                    ProEval.getSlippageScores(locus_distrib, curr["ref_len"])[0]
                )
            self.assertEqual(
                round(ProEval.getThresholdFromScores(pro_p), 6),
                round(curr["expected"], 6)
            )
        # self.assertEqual(
        #     round(ProEval.getThresholdFromScores([0.147127, 0.158391, 0.174515, 0.154106, 0.150265, 0.152119, 0.154577, 0.149493, 0.151734, 0.133869, 0.147504, 0.162312, 0.142707, 0.149554, 0.14993, 0.133616, 0.146007, 0.144878, 0.147363, 0.148772, 0.201536, 0.141008, 0.160981, 0.192323, 0.129846, 0.133303, 0.140621, 0.152649, 0.214573, 0.14024, 0.142125, 0.15114, 0.156582, 0.145829, 0.136654, 0.166352, 0.152788, 0.147627, 0.15286, 0.149677, 0.17488, 0.1588, 0.149415, 0.144182, 0.161492, 0.153196, 0.138065, 0.191469, 0.153171, 0.135274, 0.144125, 0.152325, 0.147295, 0.152189, 0.148673, 0.148799, 0.144999, 0.140025, 0.151651, 0.149045, 0.146288, 0.166864, 0.159557, 0.161954, 0.196679, 0.162283, 0.146167, 0.168555, 0.150034, 0.159113, 0.148918, 0.124418, 0.151375, 0.154594, 0.150862, 0.165631, 0.145705, 0.135396, 0.153937, 0.150396, 0.148585, 0.183152, 0.180425, 0.163239, 0.211433, 0.155945, 0.175925, 0.145602, 0.162011, 0.12887, 0.144698, 0.143911, 0.159597], 0), 6),
        #     0.203375
        # ) # Problem with floating point


class TestProIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.txt")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.txt")

        with open(self.tmp_in, "w") as writer:
            writer.write("""chromosome	location	left_flank_bases	repeat_times	repeat_unit_bases	right_flank_bases	pro_p	pro_q	CovReads	threshold
2	47414420	CAGGT	27	A	GGGTT	0.159597	0.000641	404	0.203254
2	95183613	TCCTA	23	T	GTGAG	0.083443	0.005242	99	0.107100
""")
        self.records = [
            (
                Locus.fromDict({
                    "position": "2:47414420",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": Status.none,
                            "data": {
                                "pro_p": 0.159597,
                                "pro_q": 0.000641,
                                "depth": 404
                            }
                        }
                    }
                }),
                {
                    "chromosome": "2",
                    "location": 47414420,
                    "left_flank_bases": "CAGGT",
                    "MSIsensor-pro": {"pro_p_threshold": 0.203254},
                    "repeat_times": 27,
                    "repeat_unit_bases": "A",
                    "right_flank_bases": "GGGTT"
                }
            ), (
                Locus.fromDict({
                    "position": "2:95183613",
                    "results": {
                        "MSIsensor-pro_pro": {
                            "status": Status.none,
                            "data": {
                                "pro_p": 0.083443,
                                "pro_q": 0.005242,
                                "depth": 99
                            }
                        }
                    }
                }),
                {
                    "chromosome": "2",
                    "location": 95183613,
                    "left_flank_bases": "TCCTA",
                    "MSIsensor-pro": {"pro_p_threshold": 0.107100},
                    "repeat_times": 23,
                    "repeat_unit_bases": "T",
                    "right_flank_bases": "GTGAG"
                }
            )
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testLoci(self):
        self.assertEqual(
            ProIO.loci(self.tmp_in),
            [loci.position for loci, meta in self.records]
        )

    def testRead(self):
        with ProIO(self.tmp_in) as reader:
            self.assertEqual(
                [(toDict(locus), meta) for locus, meta in reader],
                [(toDict(locus), meta) for locus, meta in self.records]
            )

    def testWrite(self):
        with ProIO(self.tmp_out, "w") as writer:
            for locus, meta in self.records:
                writer.write(locus, meta)
        with open(self.tmp_in) as reader:
            expected = reader.readlines()
        with open(self.tmp_out) as reader:
            observed = reader.readlines()
        self.assertEqual(observed, expected)


def toDict(obj):
    return json.loads(
        json.dumps(obj, default=lambda o: getattr(o, '__dict__', str(o)))
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
