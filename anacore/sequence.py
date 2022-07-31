# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing sequences."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from textwrap import wrap


class Alphabet:
    pass


class AA3LettersAlphabet(Alphabet):
    """Static methods and attributes to manage amino acid alphabet in three letters."""
    codons_by_aa = {
        "Ter": ["TAA", "TGA", "TAG"],
        "Ala": ["GCT", "GCC", "GCA", "GCG"],
        "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "Asn": ["AAT", "AAC"],
        "Asp": ["GAT", "GAC"],
        "Cys": ["TGT", "TGC"],
        "Gln": ["CAA", "CAG"],
        "Glu": ["GAA", "GAG"],
        "Gly": ["GGT", "GGC", "GGA", "GGG"],
        "His": ["CAT", "CAC"],
        "Ile": ["ATT", "ATC", "ATA"],
        "Leu": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
        "Lys": ["AAA", "AAG"],
        "Met": ["ATG"],
        "Phe": ["TTT", "TTC"],
        "Pro": ["CCT", "CCC", "CCA", "CCG"],
        "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "Thr": ["ACT", "ACC", "ACA", "ACG"],
        "Trp": ["TGG"],
        "Tyr": ["TAT", "TAC"],
        "Val": ["GTT", "GTC", "GTA", "GTG"]
    }

    one_by_three = {
        "Ter": "*",
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Asx": "B",
        "Glx": "Z",
        "Xle": "J",
        "Sec": "U",
        "Pyl": "O",
        "Xaa": "X"
    }

    words = set(codons_by_aa.keys())

    @staticmethod
    def isValid(seq):
        """
        Return True if the sequence contains only elements of the amino acid alphabet.

        :param seq: Sequence to validate.
        :type seq: str
        :return: True if the sequence contains only elements of the amino acid alphabet.
        :rtype: bool
        """
        is_valid = True
        try:
            for aa in wrap(seq, 3):
                if aa.capitalize() not in __class__.words:
                    is_valid = False
                    break
        except Exception:
            pass
        return is_valid

    @staticmethod
    def toOneLetter(seq):
        """
        Return sequence of amino acids in on letter alphabet.

        :param seq: Sequence to convert.
        :type seq: str
        :return: Sequence of amino acids in on letter alphabet.
        :rtype: str
        """
        new_seq = ""
        for aa in wrap(seq, 3):
            new_seq += __class__.one_by_three[aa]
        return new_seq


class CodonAlphabet(Alphabet):
    """Static methods and attributes to manage codon alphabet."""
    aa_by_codon = {curr_codon: aa for aa, codons in AA3LettersAlphabet.codons_by_aa.items() for curr_codon in codons}

    words = set(aa_by_codon.keys())

    @staticmethod
    def isValid(seq):
        """
        Return True if the sequence contains only elements of the codon alphabet.

        :param seq: Sequence to validate.
        :type seq: str
        :return: True if the sequence contains only elements of the codon alphabet.
        :rtype: bool
        """
        is_valid = True
        try:
            for codon in wrap(seq, 3):
                if codon.upper() not in __class__.words:
                    is_valid = False
                    break
        except Exception:
            pass
        return is_valid

    @staticmethod
    def translate(seq):
        """
        Return amino acids sequence from codon sequence.

        :param seq: Sequence to translate.
        :type seq: str
        :return: Translated sequence.
        :rtype: str
        """
        return "".join([CodonAlphabet.aa_by_codon[codon.upper()] for codon in wrap(seq, 3)])


class DNAAlphabet(Alphabet):
    """Static methods and attributes to manage DNA alphabet."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a', 'n': 'n',
        'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
        'w': 'w', 's': 's', 'm': 'k', 'k': 'm', 'r': 'y', 'y': 'r', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'
    }

    words = set(complement.keys())

    @staticmethod
    def isValid(seq):
        """
        Return True if the sequence contains only elements of the DNA alphabet.

        :param seq: Sequence to validate.
        :type seq: str
        :return: True if the sequence contains only elements of the DNA alphabet.
        :rtype: bool
        """
        is_valid = True
        for nt in seq:
            if nt not in __class__.words:
                is_valid = False
                break
        return is_valid

    @staticmethod
    def revCom(seq):
        """
        Return the reverse complement.

        :param seq: Sequence to reverse complement.
        :type seq: str
        :return: The reverse complement.
        :rtype: str
        """
        return "".join([__class__.complement[base] for base in seq[::-1]])


def getShortestRepeatUnit(seq):
    """
    Return the shortest repeat unit that can be used to reproduce the entire sequence or None if it does not exist.

    :param str: Complete sequence.
    :type seq: str
    :return: The shortest repeat unit that can be used to reproduce the entire sequence or None if it does not exist.
    :rtype: str
    """
    uc_seq = seq.upper()
    seq_len = len(seq)
    repeat_unit = None
    for unit_size in range(1, int(seq_len / 2) + 1):
        if repeat_unit is None and seq_len % unit_size == 0:
            eval_unit = seq[:unit_size]
            chunks = {chunk for chunk in wrap(uc_seq, unit_size)}
            if len(chunks - {eval_unit}) == 0:
                repeat_unit = eval_unit
    return repeat_unit


class RNAAlphabet(Alphabet):
    """Static methods and attributes to manage RNA alphabet."""
    complement = {
        'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N',
        'a': 'u', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a', 'n': 'n',
        'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
        'w': 'w', 's': 's', 'm': 'k', 'k': 'm', 'r': 'y', 'y': 'r', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'
    }

    words = set(complement.keys())

    @staticmethod
    def isValid(seq):
        """
        Return True if the sequence contains only elements of the RNA alphabet.

        :param seq: Sequence to validate.
        :type seq: str
        :return: True if the sequence contains only elements of the RNA alphabet.
        :rtype: bool
        """
        is_valid = True
        for nt in seq:
            if nt not in __class__.words:
                is_valid = False
                break
        return is_valid

    @staticmethod
    def revCom(seq):
        """
        Return the reverse complement.

        :param seq: Sequence to reverse complement.
        :type seq: str
        :return: The reverse complement.
        :rtype: str
        """
        return "".join([__class__.complement[base] for base in seq[::-1]])


class Sequence:
    def __init__(self, id, string, description=None, quality=None):
        """
        Build and return an instance of Sequence.

        :param id: Id of the sequence.
        :type id: str
        :param string: Sequence of the sequence.
        :type string: str
        :param description: The sequence description.
        :type description: str
        :param quality: The quality of the sequence (same length as string).
        :type quality: str
        :return: The new instance.
        :rtype: Sequence
        """
        self.id = id
        self.description = description
        self.string = string
        self.quality = quality

    def dnaRevCom(self):
        """
        Return the sequence corresponding to the DNA reverse complement.

        :return: The reverse complement.
        :rtype: Sequence
        """
        return Sequence(
            self.id,
            DNAAlphabet.revCom(self.string),
            self.description,
            (None if self.quality is None else self.quality[::-1])
        )

    def rnaRevCom(self):
        """
        Return the sequence corresponding to the RNA reverse complement.

        :return: The reverse complement.
        :rtype: Sequence
        """
        return Sequence(
            self.id,
            RNAAlphabet.revCom(self.string),
            self.description,
            (None if self.quality is None else self.quality[::-1])
        )
