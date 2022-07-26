#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.sequence import AA3LettersAlphabet, CodonAlphabet, DNAAlphabet, getShortestRepeatUnit, RNAAlphabet, Sequence


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestAA3LettersAlphabet(unittest.TestCase):
    def testIsValid(self):
        self.assertTrue(AA3LettersAlphabet.isValid("AlaGlyHisTer"))
        self.assertFalse(AA3LettersAlphabet.isValid("AlaGlaHisTer"))  # Invalid aa
        self.assertFalse(AA3LettersAlphabet.isValid("AlaHiTer"))  # Invalid length

    def testToOneLetter(self):
        self.assertEqual(
            AA3LettersAlphabet.toOneLetter("AlaGlyHisTer"),
            "AGH*"
        )
        with self.assertRaises(Exception):
            AA3LettersAlphabet.toOneLetter("AlaGlaHisTer")  # Invalid aa
        with self.assertRaises(Exception):
            AA3LettersAlphabet.toOneLetter("AlaHiTer")  # Invalid length


class TestCodonAlphabet(unittest.TestCase):
    def testIsValid(self):
        self.assertTrue(CodonAlphabet.isValid("GCCGGGATCTGA"))
        self.assertFalse(CodonAlphabet.isValid("GCCGIGATCTGA"))  # Invalid nt
        self.assertFalse(CodonAlphabet.isValid("GCCGGGACTGA"))  # Invalid length

    def testTranslate(self):
        self.assertEqual(
            CodonAlphabet.translate("GCCGGGATCTGA"),
            "AlaGlyIleTer"
        )
        with self.assertRaises(Exception):
            CodonAlphabet.translate("GCCGIGATCTGA")  # Invalid nt
        with self.assertRaises(Exception):
            CodonAlphabet.translate("GCCGGGACTGA")  # Invalid length


class TestDNAAlphabet(unittest.TestCase):
    def testIsValid(self):
        self.assertTrue(DNAAlphabet.isValid("GCCGGGATCTGM"))
        self.assertFalse(DNAAlphabet.isValid("GCCGUGAICTGA"))  # Invalid nt

    def testRevCom(self):
        self.assertEqual(
            DNAAlphabet.revCom("GCCGGGATCTGM"),
            "KCAGATCCCGGC"
        )
        with self.assertRaises(Exception):
            DNAAlphabet.revCom("GCCGUGAICTGA")  # Invalid nt


class TestGetShortestRepeatUnit(unittest.TestCase):
    def test(self):
        self.assertEqual(
            getShortestRepeatUnit("A"),
            None
        )
        self.assertEqual(
            getShortestRepeatUnit("AA"),
            "A"
        )
        self.assertEqual(
            getShortestRepeatUnit("ATATATAT"),
            "AT"
        )
        self.assertEqual(
            getShortestRepeatUnit("AAAAA"),
            "A"
        )
        self.assertEqual(
            getShortestRepeatUnit("ATCGAATCGA"),
            "ATCGA"
        )
        self.assertEqual(
            getShortestRepeatUnit("ATATATAG"),
            None
        )


class TestRNAAlphabet(unittest.TestCase):
    def testIsValid(self):
        self.assertTrue(RNAAlphabet.isValid("GCCGGGAUCUGM"))
        self.assertFalse(RNAAlphabet.isValid("GCCGIGAUCUGA"))  # Invalid nt

    def testRevCom(self):
        self.assertEqual(
            RNAAlphabet.revCom("GCCGGGAUCUGM"),
            "KCAGAUCCCGGC"
        )
        with self.assertRaises(Exception):
            RNAAlphabet.revCom("GCCGIGAUCUGA")  # Invalid nt


class TestSequence(unittest.TestCase):
    def testDnaRevCom(self):
        # Without quality
        test_record = Sequence("id1", "AtMGCUN").dnaRevCom()
        expected_revcom = Sequence("id1", "NAGCKaT")
        self.assertTrue(cmpSequences(test_record, expected_revcom))
        # With quality
        test_record = Sequence("id1", "AtMGCUN", "", "18AAEGH").dnaRevCom()
        expected_revcom = Sequence("id1", "NAGCKaT", "", "HGEAA81")
        self.assertTrue(cmpSequences(test_record, expected_revcom))

    def testRnaRevCom(self):
        # Without quality
        test_record = Sequence("id1", "AuMGCUN").rnaRevCom()
        expected_revcom = Sequence("id1", "NAGCKaU")
        self.assertTrue(cmpSequences(test_record, expected_revcom))
        # With quality
        test_record = Sequence("id1", "AuMGCUN", "", "18AAEGH").rnaRevCom()
        expected_revcom = Sequence("id1", "NAGCKaU", "", "HGEAA81")
        self.assertTrue(cmpSequences(test_record, expected_revcom))


def cmpSequences(rec_1, rec_2):
    are_equal = False
    if rec_1.id == rec_2.id:
        if rec_1.description == rec_2.description:
            if rec_1.string == rec_2.string:
                if rec_1.quality == rec_2.quality:
                    are_equal = True
    return are_equal


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
