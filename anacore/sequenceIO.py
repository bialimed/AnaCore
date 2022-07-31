# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing sequence files and their indices."""

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '2.6.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.abstractFile import isGzip
from anacore.sequence import Sequence, DNAAlphabet
from anacore.sv import SVIO
import gzip
import re


class FaidxRecord:
    """Record for faidx (see http://www.htslib.org/doc/faidx.html)."""

    def __init__(self, name, length, offset, line_bases, line_width, qual_offset=None):
        """
        Build and return an instance of FaidxRecord.

        :param name: Name of this reference sequence.
        :type name: str
        :param length: Total length of this reference sequence, in bases.
        :type length: int
        :param offset: Offset in the FASTA/FASTQ file of this sequence's first base.
        :type offset: int
        :param line_bases: The number of bases on each line.
        :type line_bases: int
        :param line_width: The number of bytes in each line, including the newline.
        :type line_width: int
        :param qual_offset: Offset of sequence's first quality within the FASTQ file.
        :type qual_offset: int
        :return: The new instance.
        :rtype: FaidxRecord
        """
        self.name = name
        self.length = length
        self.offset = offset
        self.line_bases = line_bases
        self.line_width = line_width
        self.qual_offset = qual_offset


class Faidx(SVIO):
    """Class to read and write a faidx file (see http://www.htslib.org/doc/faidx.html)."""

    def __init__(self, filepath, mode="r"):
        """
        Build and return an instance of Faidx.

        :param filepath: Path to the file.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: Faidx
        """
        super().__init__(filepath, mode=mode, separator="\t", title_starter=None, has_title=False)

    def _parseLine(self):
        """
        Return a structured record from the current line.

        :return: The record defined by the current line.
        :rtype: FaidxItem
        """
        fields = []
        for idx, elt in enumerate(self.current_line.split('\t')):
            elt = elt.strip()
            if idx != 0:
                elt = int(elt)
            fields.append(elt)
        return FaidxRecord(*fields)

    def readById(self):
        """
        Return all the record from the file.

        :return: By sequence ID his FaiItem.
        :rtype: dict
        """
        return {record.name: record for record in self}


class SequenceFileReader(object):
    @staticmethod
    def factory(filepath):
        if FastqIO.isValid(filepath):
            return FastqIO(filepath)
        elif FastaIO.isValid(filepath):
            return FastaIO(filepath)
        else:
            raise IOError("The file " + filepath + " does not have a valid format for 'SequenceFileReader'.")


class FastqIO:
    """Class to read and write in fastq file (https://en.wikipedia.org/wiki/FASTQ_format)."""

    def __init__(self, filepath, mode="r"):
        """
        Build and return an instance of FastqIO.

        :param filepath: Path to the file.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: FastqIO
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and isGzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 1

    def __del__(self):
        self.close()

    def close(self):
        """Close file handle it is opened."""
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        is_end = False
        while not is_end:
            record = self.nextSeq()
            if record is not None:
                yield(record)
            else:
                is_end = True

    def nextSeq(self):
        """
        Return the next sequence.

        :return: The next sequence or None if it is the end of file.
        :rtype: anacore.sequence.Sequence
        """
        seq_record = None
        try:
            prev_file_pos = self.file_handle.tell()
            header = self.file_handle.readline().strip()
            new_file_pos = self.file_handle.tell()
            if prev_file_pos != new_file_pos:
                # Header
                fields = header[1:].split(None, 1)
                seq_id = fields[0]
                seq_desc = fields[1] if len(fields) == 2 else None
                self.current_line_nb += 1
                # Sequence
                seq_str = self.file_handle.readline().strip()
                self.current_line_nb += 1
                # Separator
                self.file_handle.readline()
                self.current_line_nb += 1
                # Quality
                seq_qual = self.file_handle.readline().strip()
                self.current_line_nb += 1
                # Record
                seq_record = Sequence(seq_id, seq_str, seq_desc, seq_qual)
        except Exception:
            raise IOError(
                'The line {} in "{}" cannot be parsed by {}.'.format(
                    self.current_line_nb,
                    self.filepath,
                    self.__class__.__name__
                )
            )
        return seq_record

    @staticmethod
    def qualOffset(filepath):
        """
        Return the offset used to encode the quality in the file.

        :param filepath: The file path.
        :type filepath: str
        :return: The offset 33 for sanger and Illumina >=1.8, 64 for Solexa and Illumina <1.8 or None if the offset cannot be determined.
        :rtype: int
        """
        offset = None
        nb_qual = 0
        count_by_qual = {elt: 0 for elt in range(-5, 127)}
        with FastqIO(filepath) as FH_in:
            record = FH_in.nextSeq()
            while record and offset is None:
                for curr_nt, curr_ascii in zip(record.string, record.quality):
                    num_ascii = ord(curr_ascii)
                    if num_ascii < 59:
                        offset = 33
                        break
                    if curr_nt != "N":
                        nb_qual += 1
                        count_by_qual[num_ascii] += 1
                record = FH_in.nextSeq()
        if offset is None:
            checked_idx = int((nb_qual / 100) * 30)  # Index of the 30 percentile in sorted qualities
            curr_count_sum = 0
            for curr_ascii, curr_count in sorted(count_by_qual.items()):
                curr_count_sum += curr_count
                if offset is None and curr_count_sum >= checked_idx:  # 30% of qualities under this point
                    if curr_ascii > 84:  # 70% of qualities are superior than Q61 with offset 33 and Q20 with offset 64
                        offset = 64
                    else:
                        offset = 33
        return offset

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in fastq format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in fastq format.
        :rtype: bool
        """
        is_valid = False
        FH_in = FastqIO(filepath)
        try:
            seq_idx = 0
            end_of_file = False
            while seq_idx < 10 and not end_of_file:
                curr_line = FH_in.file_handle.readline()
                if not curr_line:
                    end_of_file = True
                else:
                    seq_idx += 1
                    # Record ID
                    header = curr_line
                    if not header.startswith("@"):
                        raise IOError('The record {} in "{}" has an invalid header.'.format(seq_idx, filepath))
                    # Record sequence
                    unstriped_sequence = FH_in.file_handle.readline()
                    if not unstriped_sequence:  # No line
                        raise IOError('The record {} in "{}" is truncated.'.format(seq_idx, filepath))
                    if not re.search("^[A-Za-z]*$", unstriped_sequence.strip()):
                        raise IOError('The sequence {} in "{}" contains invalid characters.'.format(seq_idx, filepath))
                    # Record separator
                    unstriped_separator = FH_in.file_handle.readline()
                    if not unstriped_separator:  # No line
                        raise IOError('The record {} in "{}" is truncated.'.format(seq_idx, filepath))
                    # Record quality
                    unstriped_quality = FH_in.file_handle.readline()
                    # Sequence and quality length
                    if len(unstriped_sequence.strip()) != len(unstriped_quality.strip()):
                        raise IOError('The record {} in "{}" contains a sequence and a quality with different length.'.format(seq_idx, filepath))
            is_valid = True
        except Exception:
            pass
        finally:
            FH_in.close()
        return is_valid

    @staticmethod
    def nbSeq(filepath):
        """
        Return the number of sequences in file.

        :param filepath: Path to the file.
        :type filepath: str
        :return: The number of sequences.
        :rtype: int
        """
        handler = open
        handler_options = "r"
        if isGzip(filepath):
            handler = gzip.open
            handler_options = "rt"
        nb_lines = 0
        with handler(filepath, handler_options) as reader:
            for line in reader:
                nb_lines += 1
        return int(nb_lines / 4)

    @staticmethod
    def nbSeqAndNt(filepath):
        """
        Return the number of sequences and nucleotids in file.

        :param filepath: Path to the file.
        :type filepath: str
        :return: The number of sequences and the number of nucleotids.
        :rtype: int, int
        """
        nb_seq = 0
        nb_nt = 0
        handler = open
        handler_options = "r"
        if isGzip(filepath):
            handler = gzip.open
            handler_options = "rt"
        with handler(filepath, handler_options) as reader:
            for line in reader:
                nb_seq += 1
                nb_nt += len(reader.readline().rstrip())
                reader.readline()
                reader.readline()
        return nb_seq, nb_nt

    def write(self, sequence_record):
        """
        Write record lines in file.

        :param sequence_record: The record to write.
        :type sequence_record: anacore.sequence.Sequence
        """
        self.file_handle.write(self.seqToFastqLine(sequence_record) + "\n")
        self.current_line_nb += 1

    def seqToFastqLine(self, sequence):
        """
        Return the sequence in fastq format.

        :param sequence: The sequence to process.
        :type sequence: anacore.sequence.Sequence
        :return: The sequence.
        :rtype: str
        """
        seq = "@" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        seq += "\n" + sequence.string
        seq += "\n+"
        seq += "\n" + sequence.quality
        return seq


class FastaIO:
    """Class to read and write in fasta file (https://en.wikipedia.org/wiki/FASTA_format)."""

    def __init__(self, filepath, mode="r"):
        """
        Build and return an instance of FastaIO.

        :param filepath: Path to the file.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :return: The new instance.
        :rtype: FastaIO
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and isGzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 1
        self._next_id = None
        self._end_of_file = False

    def __del__(self):
        self.close()

    def close(self):
        """Close file handle it is opened."""
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        is_end = False
        while not is_end:
            record = self.nextSeq()
            if record is not None:
                yield(record)
            else:
                is_end = True

    def nextSeq(self):
        """
        Return the next sequence.

        :return: The next sequence.
        :rtype: anacore.sequence.Sequence
        """
        seq_record = None
        if not self._end_of_file:
            line = ""
            try:
                # First line in file
                if self.current_line_nb == 1:
                    self._next_id = self.file_handle.readline().strip()
                    self.current_line_nb += 1
                # Sequence
                seq_str = ""
                while not line.startswith('>'):
                    seq_str += line.strip()
                    line = self.file_handle.readline()
                    if not line:
                        line = None
                        self._end_of_file = True
                        break
                    self.current_line_nb += 1
                fields = self._next_id[1:].split(None, 1)
                seq_id = fields[0]
                seq_desc = fields[1].strip() if len(fields) == 2 else None
                seq_record = Sequence(seq_id, seq_str, seq_desc)
                self._next_id = line  # next seq_id
            except Exception:
                raise IOError(
                    'The line {} in "{}" cannot be parsed by {}.\ncontent: {}'.format(
                        self.current_line_nb,
                        self.filepath,
                        self.__class__.__name__,
                        line
                    )
                )
        return seq_record

    @staticmethod
    def isValid(filepath):
        """
        Return true if the file is in fasta format.

        :param filepath: The file path.
        :type filepath: str
        :return: True if the file is in fasta format.
        :rtype: bool
        """
        is_valid = False
        FH_in = FastaIO(filepath)
        try:
            seq_idx = 0
            prev_is_header = None
            end_of_file = False
            while seq_idx < 10 and not end_of_file:
                curr_line = FH_in.file_handle.readline()
                if not curr_line:
                    end_of_file = True
                else:
                    if curr_line.startswith(">"):
                        if prev_is_header:
                            raise IOError('The fasta file "{}" cotains an header without sequence.'.format(filepath))
                        prev_is_header = True
                        seq_idx += 1
                    else:
                        if seq_idx == 0:  # The fasta file do not starts with ">"
                            raise IOError('The fasta file "{}" does not start with ">".'.format(filepath))
                        prev_is_header = False
            is_valid = True
        except Exception:
            pass
        finally:
            FH_in.close()
        return is_valid

    @staticmethod
    def nbSeq(filepath):
        """
        Return the number of sequences in file.

        :param filepath: Path to the file.
        :type filepath: str
        :return: The number of sequences.
        :rtype: int
        """
        nb_seq = 0
        handler = open
        handler_options = "r"
        if isGzip(filepath):
            handler = gzip.open
            handler_options = "rt"
        with handler(filepath, handler_options) as reader:
            for line in reader:
                if line.startswith(">"):
                    nb_seq += 1
        return nb_seq

    @staticmethod
    def nbSeqAndNt(filepath):
        """
        Return the number of sequences and nucleotids in file.

        :param filepath: Path to the file.
        :type filepath: str
        :return: The number of sequences and the number of nucleotids.
        :rtype: int, int
        """
        nb_seq = 0
        nb_nt = 0
        handler = open
        handler_options = "r"
        if isGzip(filepath):
            handler = gzip.open
            handler_options = "rt"
        with handler(filepath, handler_options) as reader:
            for line in reader:
                if line.startswith(">"):
                    nb_seq += 1
                else:
                    nb_nt += len(line.rstrip())
        return nb_seq, nb_nt

    def write(self, sequence_record):
        """
        Write record lines in file.

        :param sequence_record: The record to write.
        :type sequence_record: anacore.sequence.Sequence
        """
        self.file_handle.write(self.seqToFastaLine(sequence_record) + "\n")
        self.current_line_nb += 1

    def seqToFastaLine(self, sequence):
        """
        Return the sequence in fasta format.

        :param sequence: The sequence to process.
        :type sequence: anacore.sequence.Sequence
        :return: The sequence.
        :rtype: str
        """
        header = ">" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        return header + "\n" + sequence.string


class IdxFastaIO(FastaIO):
    """Class to read and write a fasta file with an faidx index file."""

    def __init__(self, filepath, mode="r", use_cache=False, fai_path=None):
        """
        Build and return an instance of IdxFastaIO.

        :param filepath: Path to the file.
        :type filepath: str
        :param mode: Mode to open the file ('r', 'w', 'a').
        :type mode: str
        :param use_cache: If True the last sequence retrieved with get is kept in memory to speed-up successive call to get() for the same sequence.
        :type use_cache: bool
        :param fai_path: Path to the fai file.
        :type fai_path: str
        :return: The new instance.
        :rtype: IdxFastaIO
        """
        super().__init__(filepath, mode=mode)
        if mode != "r":
            raise NotImplementedError("Write and Append mode are not currently implemented for IdxFastaIO.")
        self.fai_path = self.filepath + ".fai" if fai_path is None else fai_path
        self.use_cache = use_cache
        self.cached = None  # cached sequence
        self._setIndex()

    def _setIndex(self):
        """Set self.index with content of self.fai_path."""
        self.index = Faidx(self.fai_path).readById()

    def getSub(self, id, start, end=None):
        """
        Return the selected sub of the sequence from file. Return empty sequence if start position is out. Return fitted sequence if end is out.

        :param id: The sequence ID.
        :type id: str
        :param start: The start position of the selected sub-sequence (1-based).
        :type start: int
        :param end: The end position of the selected sub-sequence (1-based). Default: The end of the sequence.
        :type end: int
        :return: The sequence selected.
        :rtype: str
        """
        if start > self.index[id].length + 1:
            return ""
        endline_marker_len = self.index[id].line_width - self.index[id].line_bases
        # Start position
        start_line_idx = int((start - 1) / self.index[id].line_bases)
        read_start = self.index[id].offset + start - 1 + start_line_idx * endline_marker_len
        # End position
        if end is None or end > self.index[id].length + 1:
             end = self.index[id].length
        end_line_idx = int((end - 1) / self.index[id].line_bases)
        read_end = self.index[id].offset + end - 1 + end_line_idx * endline_marker_len
        # Get sequence
        self.file_handle.seek(read_start)
        seq = self.file_handle.read(read_end - read_start + 1).replace("\n", "").replace("\r", "")
        return seq

    def get(self, id):
        """
        Return the sequence from file.

        :param id: The sequence ID.
        :type id: str
        :return: The sequence selected from the file.
        :rtype: anacore.sequence.Sequence
        """
        if self.cached is not None:
            if self.cached.id == id:
                return self.cached
        # The sequence is not already cached
        seq = self.getSub(id, 1, None)
        selected = Sequence(id, seq)
        if self.use_cache:
            self.cached = selected
        return selected

    # def write(self, sequence_record):
    #     """
    #     Write record lines in file.
    #
    #     :param sequence_record: The record to write.
    #     :type sequence_record: anacore.sequence.Sequence
    #     """
    #     self.file_handle.write(self.seqToFastaLine(sequence_record) + "\n")
    #     self.faidx.write(sequence_record)
    #     self.current_line_nb += 1


def getStrandedSeqFromPos(reference, pos_on_ref, strand, sequence_reader, complement_alphabet=None):
    """
    Return the stranded sequence for the positions (continuous or not).

    :param sequence_reader: File handle to the sequences file.
    :type sequence_reader: anacore.sequenceIO.IdxFastaIO
    :param pos_list: List of reference positions. These positions can be sorted/unsorted and continuous/discontinuous.
    :type pos_list: list
    :param complement_alphabet: The alphabet used to complement the sequence (default: DNA alphabet).
    :type complement_alphabet: anacore.sequence.Alphabet
    :return: Stranded sequence for the positions.
    :rtype: str
    """
    if complement_alphabet is None:
        complement_alphabet = DNAAlphabet
    seq = ""
    # Get sequence
    pos_on_ref = sorted(pos_on_ref)
    stack_start = pos_on_ref[0]
    stack_end = pos_on_ref[0]
    for curr_pos in pos_on_ref[1:]:
        if curr_pos == stack_end + 1:  # Continuous
            stack_end = curr_pos
        else:  # Dicontinuous
            seq += sequence_reader.getSub(reference, stack_start, stack_end)
            stack_start = curr_pos
            stack_end = curr_pos
    if stack_start is not None:
        seq += sequence_reader.getSub(reference, stack_start, stack_end)
    # Reverse sequence to be stranded +
    if strand == "-":
        seq = complement_alphabet.revCom(seq)
    return seq
