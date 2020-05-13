# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing TopHatFusion outputs."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.abstractFile import AbstractFile
from anacore.region import Region


class TopHatFusionIO(AbstractFile):
    """Class to manage read and write in TopHatFusionIO file (see: https://ccb.jhu.edu/software/tophat/fusion_manual.shtml)."""

    def _parseLine(self):
        """
        Return a structured record from the TopHatFusionIO current line.

        :return: The record.
        :rtype: dict
        """
        fusion, trash_1, contig_a, contig_b, depth_a, depth_b, mate_distances = [elt.strip() for elt in self.current_line.split('@')]
        chrom, break_a, break_b, orientation, nb_splitted_reads, nb_splitted_pairs, nb_pairs_splitted_reads, nb_contradict, base_cover_left, base_cover_right, trash_1 = [field.strip() for field in fusion.split("\t")]
        chrom_a, chrom_b = chrom.split("-")
        break_a = int(break_a)
        break_b = int(break_b)
        strand_a, strand_b = [("+" if elt == "f" else "-") for elt in orientation]
        return {
            "partner_a": Region(break_a, break_a, strand_a, chrom_a),
            "partner_b": Region(break_b, break_b, strand_b, chrom_b),
            "nb_splitted_reads": int(nb_splitted_reads),
            "nb_splitted_pairs": int(nb_splitted_pairs),
            "nb_pairs_splitted_reads": int(nb_pairs_splitted_reads),
            "nb_contradict": int(nb_contradict),
            "base_cover_left": int(base_cover_left),
            "base_cover_right": int(base_cover_right)
        }
