#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.7.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import gzip
import uuid
import filecmp
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.sequence import RNAAlphabet, Sequence
from anacore.sequenceIO import IdxFastaIO, FastaIO, FastqIO, getStrandedSeqFromPos


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFastaIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_multi_line = os.path.join(tmp_folder, unique_id + "_multi.fasta")
        self.tmp_multi_line_gz = os.path.join(tmp_folder, unique_id + "_multi.fasta.gz")
        self.tmp_mono_line = os.path.join(tmp_folder, unique_id + "_mono.fasta")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.fasta")

        # Expected records
        self.expected_rec = [
            Sequence("seq1", "ATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATAC", "test description"),
            Sequence("seq2", ""),
            Sequence("seq3", "ATGAAAAAAAAAAAAANTGATGAAAAAAAAAAAAANTG", "test description 2"),
            Sequence("seq4", "", "trimmed")
        ]

        # Create multi line
        content = """>seq1 test description
ATAGATAGCATCCCCCCNATAC
ATAGATAGCATCCCCCCNATAC
ATAGATAGCATCCCCCCNATAC
>seq2

>seq3 test description 2
ATGAAAAAAAAAAAAANTG
ATGAAAAAAAAAAAAANTG
>seq4 trimmed

"""
        with open(self.tmp_multi_line, "w") as FH_out:
            FH_out.write(content)

        # Create mono line
        content = """>seq1 test description
ATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATACATAGATAGCATCCCCCCNATAC
>seq2

>seq3 test description 2
ATGAAAAAAAAAAAAANTGATGAAAAAAAAAAAAANTG
>seq4 trimmed

"""
        with open(self.tmp_mono_line, "w") as FH_out:
            FH_out.write(content)
        with gzip.open(self.tmp_multi_line_gz, "wt") as FH_out:
            FH_out.write(content)

    def testNbSeq(self):
        nb_seq = FastaIO.nbSeq(self.tmp_mono_line)
        self.assertEqual(4, nb_seq)
        nb_seq = FastaIO.nbSeq(self.tmp_multi_line)
        self.assertEqual(4, nb_seq)
        nb_seq = FastaIO.nbSeq(self.tmp_multi_line_gz)
        self.assertEqual(4, nb_seq)

    def testNbSeqAndNt(self):
        nb_seq, nb_nt = FastaIO.nbSeqAndNt(self.tmp_mono_line)
        self.assertEqual(nb_seq, 4)
        self.assertEqual(nb_nt, 104)
        nb_seq, nb_nt = FastaIO.nbSeqAndNt(self.tmp_multi_line)
        self.assertEqual(nb_seq, 4)
        self.assertEqual(nb_nt, 104)
        nb_seq, nb_nt = FastaIO.nbSeqAndNt(self.tmp_multi_line_gz)
        self.assertEqual(nb_seq, 4)
        self.assertEqual(nb_nt, 104)

    def testWrite(self):
        with FastaIO(self.tmp_out, "w") as FH_out:
            for curr_rec in self.expected_rec:
                FH_out.write(curr_rec)
        self.assertTrue(FastaIO.isValid(self.tmp_out))
        self.assertTrue(filecmp.cmp(self.tmp_out, self.tmp_mono_line))

    def testIsValid(self):
        # Valid
        self.assertTrue(FastaIO.isValid(self.tmp_mono_line))
        self.assertTrue(FastaIO.isValid(self.tmp_multi_line))
        self.assertTrue(FastaIO.isValid(self.tmp_multi_line_gz))
        # Valid long file
        content = ">seq1\nATGC\n>seq2\nATGC\n>seq3\nATGC\n>seq4\nATGC\n>seq5\nATGC\n>seq6\nATGC\n>seq7\nATGC\n>seq8\nATGC\n>seq9\nATGC\n>seq10\nATGC\n>seq11\nATGC\n>seq12\nATGC"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastaIO.isValid(self.tmp_out))
        # Valid empty file
        content = ""
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastaIO.isValid(self.tmp_out))
        # Valid empty sequence
        content = ">seq1\n"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastaIO.isValid(self.tmp_out))
        # Invalid file (two headers)
        content = ">seq1\nATGC\n>seq2\n>seq3\nATGC"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(not FastaIO.isValid(self.tmp_out))
        # Invalid file (no header at the first line)
        content = "seq1\nATGC\n>seq2\nATGC"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(not FastaIO.isValid(self.tmp_out))
        # Invalid file (fastq)
        content = "@seq1\nATGC\n+\n####"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(not FastaIO.isValid(self.tmp_out))

    def testIter(self):
        with FastaIO(self.tmp_mono_line) as FH_in:
            for idx, record in enumerate(FH_in):
                self.assertTrue(cmpSequences(record, self.expected_rec[idx]))
            self.assertEqual(idx + 1, 4)
        with FastaIO(self.tmp_multi_line) as FH_in:
            for idx, record in enumerate(FH_in):
                self.assertTrue(cmpSequences(record, self.expected_rec[idx]))
            self.assertEqual(idx + 1, 4)
        with FastaIO(self.tmp_multi_line_gz) as FH_in:
            for idx, record in enumerate(FH_in):
                self.assertTrue(cmpSequences(record, self.expected_rec[idx]))
            self.assertEqual(idx + 1, 4)
        with FastaIO(self.tmp_mono_line) as FH_in_mono:
            with FastaIO(self.tmp_multi_line) as FH_in_multi:
                idx = 0
                for rec_expected, rec_mono_line, rec_multi_line in zip(self.expected_rec, FH_in_mono, FH_in_multi):
                    self.assertTrue(cmpSequences(rec_mono_line, rec_expected))
                    self.assertTrue(cmpSequences(rec_multi_line, rec_expected))
                    idx += 1
                self.assertEqual(idx, 4)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_mono_line, self.tmp_multi_line, self.tmp_multi_line_gz, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


class TestFastqIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_seq = os.path.join(tmp_folder, unique_id + "_ref.fastq")
        self.tmp_seq_gz = os.path.join(tmp_folder, unique_id + "_ref.fastq.gz")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.fastq")

        # Expected
        self.expected_rec = [
            Sequence("seq1", "ATAGATAGCATCCCCCCNATAC", "test description", "@?>=<;:9876543210##-,+"),
            Sequence("seq2", "", None, ""),
            Sequence("seq3", "ATGAAAAAAAAAAAAANTG", "test description 2", "@?>=<;:987654321#/."),
            Sequence("seq4", "", "trimmed", "")
        ]

        # Create file
        content = """@seq1 test description
ATAGATAGCATCCCCCCNATAC
+
@?>=<;:9876543210##-,+
@seq2

+

@seq3 test description 2
ATGAAAAAAAAAAAAANTG
+
@?>=<;:987654321#/.
@seq4 trimmed

+

"""
        with open(self.tmp_seq, "w") as FH_out:
            FH_out.write(content)
        with gzip.open(self.tmp_seq_gz, "wt") as FH_out:
            FH_out.write(content)

    def testQualOffset(self):
        # Illumina 1.8 with under 59
        content = """@M70265:234:000000000-CCC3N:1:1101:21165:1697 1:N:0:13
ATGTCCTTGTGCACAATGCCCTGGCTATGCAGGTACTCCAGGCCGTCAATCAGCTGACAGAAGTACCTGCGGGCAGCACACACCCGTCCTGGGGCCGAGGCCTCCCTGCCCCTCTCAGGGGCGAATTTCGACGATCGTTGCATTAACTCGC
+
-A-A@EFF9E,C9C,,C,CEEF,,6C9E,@,,,C<EEEE,,,:B7@:,,CC,,CE,,;,,,,,<9@CE,C+++@+,@,C,C,B@>+BBFE,,,+87++++8ABA=FE,B?BFDC==,,,,+6+++@BFD,+8++>>7@D,<,@@,,@7>*>
@M70265:234:000000000-CCC3N:1:1101:14142:1764 1:N:0:13
TGTCAATCAATATCAGGACAAGCAGTGTGTCCTCACGGAAAGGAGCCTGCCCTGCCTGGCCCCCGGCCCCCGCCCCACCCTGGCCCCTGCCCCGCGCACCCACCCGTTGGCCTTGCCCCCTCGGAAACGCTTCTCCCGCACCCTTGCGAAT
+
B<CB9-CF9,F9FDC,,,C8<,C8,C,C,,<CEE9C,,+,6,<,<CBF@<@EE,CFD,,@ADCF7::C@B@+@>CC,:FFC,,4CDC<,,CED+@+>+6+8+?CC+8,+,4:>,,:,:+83++++++8+83>:,33,+3+5*68,,,**1*
@M70265:234:000000000-CCC3N:1:1101:9715:1775 1:N:0:13
TCCAGGGCTTTTGTCTTCTTCCCTTTAGATTCTCTTCTTCTGTACTGCCTGTGCTTTTGCATTCTCTACACTCATCTGTGCCACCGTTTGGAAAGCTAGTGGTTCAGAGTTCTATATATTCTCGAATTTCGCCGATCGTTTCATTAACTCT
+
-8A----8FFGG,E@E@EEF<@6CFF9,,<,;C6C,6CE@C,C,CF,@C,,;,,,<,;,;,,,6,;C,,,6;;,<<E,,,6C,C<+CBA,,,,,,6C,,:,,CC@,,,,<E@F,C,,,<EA,C,,,9?E,,,8++8>+BE+559E,,5=E,"""
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertEqual(FastqIO.qualOffset(self.tmp_out), 33)
        # Illumina 1.8
        content = """@M70265:234:000000000-CCC3N:1:1102:19767:8584 1:N:0:35
TCATGACTGATATGGTAGACAGAGCCTAAACATCCCCTTAAATTGGATTAAAAAGAAATATACCTTTGTTGTTACCTTTAAATGCAAAGTTAAAATAGGCAGAAGTCTTGCCCACATCGTTGTAGGCCTTACATTCAACCGGCGAATTTCG
+
CCCCCFGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFFGGGGGGGGGGGGFGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGEGGGGFGF
@M70265:234:000000000-CCC3N:1:1102:17014:8587 1:N:0:35
TCCATAACTTCTTGCTAAGTCCTGAGCCTGTTTTGTGTCTACTGTTCTAGAAGGCAAATCACATTTATTTCCTACTAGGACCACAGGTACATCTTCAGAGTCCTTAACTCTTTTAATTTGTTCTCTGGGAAAGAGCGAATTTCGACGATCG
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@M70265:234:000000000-CCC3N:1:1102:16174:8588 1:N:0:35
CTTGAGTGAAGGACTGAGAAAATCCCTGTTCCCACTCATACAGGACTTGGGAGGTATCCACATCCTCTTCCTCAGGATTGCCTTTACCACTCTGAGAAGGAGCTGTGGTAGTGGCACCAGAATGGATTCCAGAGTCCAGGTAAGACTGCGC
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"""
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertEqual(FastqIO.qualOffset(self.tmp_out), 33)
        # Solexa
        content = """@SRR1296011.1 1 length=107
CGGCAAGTTAACAAAAAGAAAAATGGTGAATGATACCCGGTGCTGGCAATCTCGTTTAAACTACATGCAGGAACAGCAAAGGAAATCCGGCAAATTTGCGCAGTCAT
+SRR1296011.1 1 length=107
dddddaa]aafffc`c_ccc`cccf_^cddf_fcddd`ddWdd^a]daadf[fdcffaafcffcfcff]fcfffW^I^a^^KZdaffc_cWbc[cN[[X^]`a``ca
@SRR1296011.2 2 length=107
AAATTTGCCGGATTTCCTTTGCTGTTCCTGCATGTAGTTTAAACGAGATTGCCAGCACCGGGTATCATTCACCATTTTTCTTTTTGTTAACTTGCCGTCAGCCTTTT
+SRR1296011.2 2 length=107
gggggggaggggggggggfgfgggggggcgggggggc_geggfggggggggggaggggggggggdggffggfgggaggeegcgggeggggeffgac]dbcaggeab_
@SRR1296011.3 3 length=107
CTTTCTGTTCATGTGTATCTGCTGTCTCTTAGCCCAGACTTCCCGTGTCCTTTCCACTGGGCCTTTGGGAGGTCACAGGGTCTTGATGCTGTGGTCTTGATCTGCAG
+SRR1296011.3 3 length=107
fffdffgggggc_aggaggggfe_afffffgggggfgggggggggddgge_aWdaggggg]]cfffffedfeUeaacff_Wcfcc`bb]d__b^Zacaa[]\```_b"""
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertEqual(FastqIO.qualOffset(self.tmp_out), 64)

    def testWrite(self):
        with FastqIO(self.tmp_out, "w") as FH_out:
            for curr_rec in self.expected_rec:
                FH_out.write(curr_rec)
        self.assertTrue(FastqIO.isValid(self.tmp_out))
        self.assertTrue(filecmp.cmp(self.tmp_out, self.tmp_seq))

    def testNbSeq(self):
        nb_seq = FastqIO.nbSeq(self.tmp_seq)
        self.assertEqual(nb_seq, 4)

    def testNbSeqAndNt(self):
        nb_seq, nb_nt = FastqIO.nbSeqAndNt(self.tmp_seq)
        self.assertEqual(nb_seq, 4)
        self.assertEqual(nb_nt, 41)
        nb_seq, nb_nt = FastqIO.nbSeqAndNt(self.tmp_seq_gz)
        self.assertEqual(nb_seq, 4)
        self.assertEqual(nb_nt, 41)

    def testIsValid(self):
        # Valid
        self.assertTrue(FastqIO.isValid(self.tmp_seq))
        self.assertTrue(FastqIO.isValid(self.tmp_seq_gz))
        # Valid long file
        content = "@seq1\nATGC\n+\n####\n@seq2\nATGC\n+\n####\n@seq3\nATGC\n+\n####\n@seq4\nATGC\n+\n####\n@seq5\nATGC\n+\n####\n@seq6\nATGC\n+\n####\n@seq7\nATGC\n+\n####\n@seq8\nATGC\n+\n####\n@seq9\nATGC\n+\n####\n@seq10\nATGC\n+\n####\n@seq11\nATGC\n+\n####\n"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastqIO.isValid(self.tmp_out))
        # Valid empty file
        content = ""
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastqIO.isValid(self.tmp_out))
        # Valid empty sequence
        content = "@seq1\n\n+\n\n"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(FastqIO.isValid(self.tmp_out))
        # Invalid file (fasta)
        content = ">seq1\nATGC\n>seq2\nATGC"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(not FastqIO.isValid(self.tmp_out))
        # Invalid file (not a sequence)
        content = "@seq1\nAT1GC2\n+\n######"
        with open(self.tmp_out, "w") as FH_out:
            FH_out.write(content)
        self.assertTrue(not FastqIO.isValid(self.tmp_out))

    def testIter(self):
        with FastqIO(self.tmp_seq) as FH_in:
            for idx, record in enumerate(FH_in):
                self.assertTrue(cmpSequences(record, self.expected_rec[idx]))
            self.assertEqual(idx + 1, 4)
        with FastqIO(self.tmp_seq_gz) as FH_in:
            for idx, record in enumerate(FH_in):
                self.assertTrue(cmpSequences(record, self.expected_rec[idx]))
            self.assertEqual(idx + 1, 4)
        with FastqIO(self.tmp_seq) as FH_in:
            idx = 0
            for rec_expected, rec_observed in zip(self.expected_rec, FH_in):
                self.assertTrue(cmpSequences(rec_observed, rec_expected))
                idx += 1
            self.assertEqual(idx, 4)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_seq, self.tmp_seq_gz, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


def cmpSequences(rec_1, rec_2):
    are_equal = False
    if rec_1.id == rec_2.id:
        if rec_1.description == rec_2.description:
            if rec_1.string == rec_2.string:
                if rec_1.quality == rec_2.quality:
                    are_equal = True
    return are_equal


class TestIdxFastaIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fasta_idx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        self.tmp_fasta = os.path.join(tmp_folder, unique_id + ".fasta")

        # Expected
        self.expected_rec = {
            "one": Sequence("one", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"),
            "two": Sequence("two", "ATGCATGCATGCATGCATGCATGCATGC", "another chromosome")
        }

        # Create sequence file
        content_fasta = """>one
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
>two another chromosome
ATGCATGCATGCAT
GCATGCATGCATGC"""
        with open(self.tmp_fasta, "w") as FH_out:
            FH_out.write(content_fasta)

        # Create index
        content_fasta_idx = """one	66	5	30	31
two	28	98	14	15"""
        with open(self.tmp_fasta_idx, "w") as FH_out:
            FH_out.write(content_fasta_idx)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta, self.tmp_fasta_idx]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGet(self):
        # Without cache
        with IdxFastaIO(self.tmp_fasta) as FH:
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )
            self.assertEqual(
                self.expected_rec["two"].string,
                FH.get("two").string
            )
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )
        # With cache
        with IdxFastaIO(self.tmp_fasta, "r", True) as FH:
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )
            self.assertEqual(
                self.expected_rec["two"].string,
                FH.get("two").string
            )
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.get("one").string
            )


    def testGetSub(self):
        # From one
        one_len = len(self.expected_rec["one"].string)
        with IdxFastaIO(self.tmp_fasta) as FH:
            # Full length
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.getSub("one", 1, one_len)
            )
            # Full length without end parameter
            self.assertEqual(
                self.expected_rec["one"].string,
                FH.getSub("one", 1)
            )
            # First nt
            self.assertEqual(
                self.expected_rec["one"].string[0],
                FH.getSub("one", 1, 1)
            )
            # 2 firsts nt
            self.assertEqual(
                self.expected_rec["one"].string[0:2],
                FH.getSub("one", 1, 2)
            )
            # 3 firsts nt
            self.assertEqual(
                self.expected_rec["one"].string[0:3],
                FH.getSub("one", 1, 3)
            )
            # Last nt
            self.assertEqual(
                self.expected_rec["one"].string[-1],
                FH.getSub("one", one_len, one_len)
            )
            # 2 lasts nt
            self.assertEqual(
                self.expected_rec["one"].string[-2:],
                FH.getSub("one", one_len - 1)
            )
            # 3 lasts nt
            self.assertEqual(
                self.expected_rec["one"].string[-3:],
                FH.getSub("one", one_len - 2)
            )
            # 1 nt on start of second line
            self.assertEqual(
                self.expected_rec["one"].string[30],
                FH.getSub("one", 31, 31)
            )
            # 1 nt in second line
            self.assertEqual(
                self.expected_rec["one"].string[31],
                FH.getSub("one", 32, 32)
            )
            # 11 nt in second line
            self.assertEqual(
                self.expected_rec["one"].string[31:42],
                FH.getSub("one", 32, 42)
            )
            # nt on first and second line
            self.assertEqual(
                self.expected_rec["one"].string[29:35],
                FH.getSub("one", 30, 35)
            )
            # nt on second and third line
            self.assertEqual(
                self.expected_rec["one"].string[57:63],
                FH.getSub("one", 58, 63)
            )
            # nt on first, second and third line
            self.assertEqual(
                self.expected_rec["one"].string[29:63],
                FH.getSub("one", 30, 63)
            )
            # End of first
            self.assertEqual(
                self.expected_rec["one"].string[63:66],
                FH.getSub("one", 64, 66)
            )
            # Ends out of first
            self.assertEqual(
                self.expected_rec["one"].string[63:66],
                FH.getSub("one", 64, 80)
            )
            # Starts out of first
            self.assertEqual(
                "",
                FH.getSub("one", 67, 80)
            )


class TestUtils(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_fasta_idx = os.path.join(tmp_folder, unique_id + ".fasta.fai")
        self.tmp_fasta = os.path.join(tmp_folder, unique_id + ".fasta")

        # Create sequence file
        content_fasta = """>one
ATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGCATGCATGCATGC
ATGCAT
>two another chromosome
ATGCATGCATGCAT
GCATGCATGCATGC"""
        with open(self.tmp_fasta, "w") as FH_out:
            FH_out.write(content_fasta)

        # Create index
        content_fasta_idx = """one	66	5	30	31
two	28	98	14	15"""
        with open(self.tmp_fasta_idx, "w") as FH_out:
            FH_out.write(content_fasta_idx)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta, self.tmp_fasta_idx]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testGetStrandedSeqFromPos(self):
        with IdxFastaIO(self.tmp_fasta) as reader:
            self.assertEqual(
                "ATGCCAT",
                getStrandedSeqFromPos("one", [5,6,7,8,20,21,34], "+", reader)
            )
            self.assertEqual(
                "ATGGCAT",
                getStrandedSeqFromPos("one", [5,6,7,8,20,21,34], "-", reader)
            )
            self.assertEqual(
                "AUGGCAU",
                getStrandedSeqFromPos("one", [5,6,7,8,20,21,34], "-", reader, RNAAlphabet)
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
