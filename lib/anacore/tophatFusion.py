#
# Copyright (C) 2017 IUCT-O
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

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
from anacore.abstractFile import isGzip
from anacore.bioStructures import Reference, Region


class TophatFusionIO:
    def __init__(self, filepath, mode="r"):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and isGzip(filepath)):
            self.file_handle = gzip.open(filepath, mode + "t")
        else:
            self.file_handle = open(filepath, mode)
        self.current_line_nb = 0
        self.current_line = None

    def __del__(self):
        self.close()

    def __iter__(self):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            if self.current_line.startswith('#'):
                continue
            try:
                fusion_record = self._parse_line()
            except:
                raise IOError("The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                              "Line content : " + self.current_line)
            else:
                yield fusion_record

    def close( self ):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _parse_line(self):
        """
        @summary: Returns a structured record from the TophatFusion current line.
        @return: [dict] The fusion described by the current line.
        """
        fusion, trash_1, contig_a, contig_b, depth_a, depth_b, mate_distances = [elt.strip() for elt in self.current_line.split('@')]
        chrom, break_a, break_b, orientation, nb_splitted_reads, nb_splitted_pairs, nb_pairs_splitted_reads, nb_contradict, base_cover_left, base_cover_right, trash_1 = [field.strip() for field in fusion.split("\t")]
        chrom_a, chrom_b = chrom.split("-")
        break_a = int(break_a)
        break_b = int(break_b)
        strand_a, strand_b = [("+" if elt == "f" else "-") for elt in orientation]

        return {
            "partner_a": Region(break_a, break_a, strand_a, Reference(chrom_a, chrom_a, None)),
            "partner_b": Region(break_b, break_b, strand_b, Reference(chrom_b, chrom_b, None)),
            "nb_splitted_reads": int(nb_splitted_reads),
            "nb_splitted_pairs": int(nb_splitted_pairs),
            "nb_pairs_splitted_reads": int(nb_pairs_splitted_reads),
            "nb_contradict": int(nb_contradict),
            "base_cover_left": int(base_cover_left),
            "base_cover_right": int(base_cover_right)
        }
