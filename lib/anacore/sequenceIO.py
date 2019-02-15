#
# Copyright (C) 2014 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import gzip


def is_gzip(file):
    """
    @return: [bool] True if the file is gziped.
    @param file : [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open(file)
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip


class Sequence:
    def __init__(self, id, string, description=None, quality=None):
        """
        @param id: [str] Id of the sequence.
        @param string: [str] Sequence of the sequence.
        @param description: [str] The sequence description.
        @param quality: [str] The quality of the sequence (same length as string).
        """
        self.id = id
        self.description = description
        self.string = string
        self.quality = quality


class SequenceFileReader(object):
    @staticmethod
    def factory(filepath):
        if FastqIO.is_valid(filepath):
            return FastqIO(filepath)
        elif FastaIO.is_valid(filepath):
            return FastaIO(filepath)
        else:
            raise IOError("The file " + filepath + " does not have a valid format for 'SequenceFileReader'.")


class FastqIO:
    def __init__(self, filepath, mode="r"):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 1
        self.current_line = None

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        seq_id = None
        seq_desc = None
        seq_str = None
        seq_qual = None
        try:
            for line in self.file_handle:
                line = line.rstrip()
                if (self.current_line_nb % 4) == 1:
                    fields = line[1:].split(None, 1)
                    seq_id = fields[0]
                    seq_desc = fields[1] if len(fields) == 2 else None
                elif (self.current_line_nb % 4) == 2:
                    seq_str = line
                elif (self.current_line_nb % 4) == 0:
                    seq_qual = line
                    yield Sequence(seq_id, seq_str, seq_desc, seq_qual)
                self.current_line_nb += 1
        except:
            raise IOError("The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".")

    def next_seq(self):
        """
        @summary: Returns the next sequence.
        @return: [Sequence] The next sequence or None if it is the end of file.
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
                separator = self.file_handle.readline()
                self.current_line_nb += 1
                # Quality
                seq_qual = self.file_handle.readline().strip()
                self.current_line_nb += 1
                # Record
                seq_record = Sequence(seq_id, seq_str, seq_desc, seq_qual)
        except:
            raise IOError("The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".")
        return seq_record

    @staticmethod
    def qual_offset(filepath):
        """
        @summary: Returns the offset used to encode the quality in the file.
        @filepath: [str] The file path
        @return: [int] The offset: 33 for sanger and Illumina >=1.8, 64 for Solexa and Illumina <1.8 or None if the offset cannot be determined.
        """
        offset = None
        nb_qual = 0
        count_by_qual = {elt: 0 for elt in range(-5, 127)}
        with FastqIO(filepath) as FH_in:
            record = FH_in.next_seq()
            while record and offset is None:
                for curr_ascii in record.quality:
                    num_ascii = ord(curr_ascii)
                    nb_qual += 1
                    count_by_qual[num_ascii] += 1
                    if num_ascii < 59:
                        offset = 33
                        break
                record = FH_in.next_seq()
        if offset is None:
            checked_idx = int((float(nb_qual)/100)*20)  # Index of the 20 percentile in sorted qualities
            curr_idx = 0
            for curr_ascii in count_by_qual:
                curr_idx += count_by_qual[curr_ascii]
                if curr_idx >= checked_idx:  # 20% of qualities before this point
                    if curr_ascii > 84:  # 80% of qualities are superior than Q61 with offset 33 and Q20 with offset 64
                        offset = 64
                        break
        return offset

    @staticmethod
    def is_valid(filepath):
        is_valid = False
        FH_in = FastqIO(filepath)
        try:
            seq_idx = 0
            header = FH_in.file_handle.readline()
            while seq_idx < 10 and header:
                if not header.startswith("@"):
                    raise IOError("The line '" + str(header) + "' in '" + filepath + "' is not a fastq header.")
                unstriped_sequence = FH_in.file_handle.readline()
                if unstriped_sequence == "":  # No line
                    raise IOError(filepath + "' is not a fastq.")
                unstriped_separator = FH_in.file_handle.readline()
                if unstriped_separator == "":  # No line
                    raise IOError(filepath + "' is not a fastq.")
                unstriped_quality = FH_in.file_handle.readline()
                if unstriped_quality == "":  # No line
                    raise IOError(filepath + "' is not a fastq.")
                if len(unstriped_sequence.strip()) != len(unstriped_quality.strip()):
                    raise IOError(filepath + "' is not a fastq.")
                header = FH_in.file_handle.readline()
                seq_idx += 1
            is_valid = True
        except:
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
        nb_line = 0
        if is_gzip(filepath):
            with gzip.open(filepath, "rt") as FH_in:
                for line in FH_in:
                    nb_line += 1
        else:
            with open(filepath, "r") as FH_in:
                for line in FH_in:
                    nb_line += 1
        return nb_line / 4

    def write(self, sequence_record):
        self.file_handle.write(self.seqToFastqLine(sequence_record) + "\n")

    def seqToFastqLine(self, sequence):
        """
        @summary: Returns the sequence in fastq format.
        @param sequence: [Sequence] The sequence to process.
        @return: [str] The sequence.
        """
        seq = "@" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        seq += "\n" + sequence.string
        seq += "\n+"
        seq += "\n" + sequence.quality
        return seq


class FastaIO:
    def __init__(self, filepath, mode="r"):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 1
        self.current_line = None

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        seq_id = None
        seq_desc = None
        seq_str = None
        try:
            for line in self.file_handle:
                line = line.rstrip()
                self.current_line_nb += 1
                if line.startswith('>'):
                    if seq_id is not None:
                        seq_record = Sequence(seq_id, seq_str, seq_desc)
                        yield seq_record
                    # New seq
                    fields = line[1:].split(None, 1)
                    seq_id = fields[0]
                    seq_desc = fields[1] if len(fields) == 2 else None
                    seq_str = ""
                else:
                    seq_str += line
            if seq_id is not None:
                seq_record = Sequence(seq_id, seq_str, seq_desc)
                yield seq_record
        except:
            raise IOError("The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".")

    def next_seq(self):
        """
        @summary: Returns the next sequence.
        @return: [Sequence] The next sequence.
        """
        seq_record = None
        line = ""
        try:
            # First line in file
            if self.current_line_nb == 1:
                self.next_id = self.file_handle.readline().strip()
                self.current_line_nb += 1
            # Sequence
            seq_str = ""
            while not line.startswith('>'):
                seq_str += line.strip()
                line = self.file_handle.readline()
                if not line:
                    break
                self.current_line_nb += 1
            fields = self.next_id[1:].split(None, 1)
            seq_id = fields[0]
            seq_desc = fields[1] if len(fields) == 2 else None
            seq_record = Sequence(seq_id, seq_str, seq_desc)
            self.next_id = line  # next seq_id
        except:
            raise IOError("The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\ncontent : " + line)
        return seq_record

    @staticmethod
    def is_valid(filepath):
        is_valid = False
        FH_in = FastaIO(filepath)
        try:
            seq_idx = 0
            header = FH_in.file_handle.readline()
            while seq_idx < 10 and header:
                if not header.startswith(">"):
                    raise IOError("The line '" + str(header) + "' in '" + filepath + "' is not a fasta header.")
                previous_is_header = True
                current_line = FH_in.file_handle.readline()
                while not current_line.startswith(">") and current_line:
                    previous_is_header = False
                    current_line = FH_in.file_handle.readline()
                if previous_is_header:
                    raise IOError(filepath + "' is not a fasta.")
                header = current_line
                seq_idx += 1
            is_valid = True
        except:
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
        if is_gzip(filepath):
            with gzip.open(filepath, "rt") as FH_in:
                for line in FH_in:
                    if line.startswith(">"):
                        nb_seq += 1
        else:
            with open(filepath, "r") as FH_in:
                for line in FH_in:
                    if line.startswith(">"):
                        nb_seq += 1
        return nb_seq

    def write(self, sequence_record):
        self.file_handle.write(self.seqToFastaLine(sequence_record) + "\n")

    def seqToFastaLine(self, sequence):
        """
        @summary: Returns the sequence in fasta format.
        @param sequence: [Sequence] The sequence to process.
        @return: [str] The sequence.
        """
        header = ">" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        return header + "\n" + sequence.string
