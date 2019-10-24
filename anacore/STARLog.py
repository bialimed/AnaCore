# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing STAR log."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


class STARLog:
    """
    Reader of log file produced by STAR.

    :Example of content:
        .. highlight:: python
        .. code-block:: python

                                         Started job on |       Jan 06 13:04:49
                                     Started mapping on |       Jan 06 13:06:46
                                            Finished on |       Jan 06 13:10:10
               Mapping speed, Million of reads per hour |       237.83

                                  Number of input reads |       13477055
                              Average input read length |       150
                                            UNIQUE READS:
                           Uniquely mapped reads number |       7988755
                                Uniquely mapped reads % |       59.28%
                                  Average mapped length |       150.47
                               Number of splices: Total |       2900295
                    Number of splices: Annotated (sjdb) |       2866020
                               Number of splices: GT/AG |       2870546
                               Number of splices: GC/AG |       22895
                               Number of splices: AT/AC |       2759
                       Number of splices: Non-canonical |       4095
                              Mismatch rate per base, % |       0.26%
                                 Deletion rate per base |       0.01%
                                Deletion average length |       1.36
                                Insertion rate per base |       0.00%
                               Insertion average length |       1.41
                                     MULTI-MAPPING READS:
                Number of reads mapped to multiple loci |       4957483
                     % of reads mapped to multiple loci |       36.78%
                Number of reads mapped to too many loci |       23143
                     % of reads mapped to too many loci |       0.17%
                                          UNMAPPED READS:
               % of reads unmapped: too many mismatches |       0.00%
                         % of reads unmapped: too short |       3.72%
                             % of reads unmapped: other |       0.05%
    """

    def __init__(self, path):
        """
        Build and return an instance of STARLog.

        :param path: Path to STAR log file.
        :type path: str
        :return: The new instance.
        :rtype: STARLog
        """
        self.filepath = path
        self.input = {"nb": 0, "average_length": 0}
        self.unique = {
            "nb": 0,  # Uniquely mapped reads number |       7988755
            "average_mapped_length": 0,  # Average mapped length |       150.47
            "nb_spliced": 0,  # Number of splices: Total |       2900295
            "nb_spliced_annotated": 0,  # Number of splices: Annotated (sjdb) |       2866020
            "nb_spliced_non_canonical": 0,  # Number of splices: Non-canonical |       4095
        }
        self.multi = {
            "nb": 0,  # Number of reads mapped to multiple loci |       4957483       + Number of reads mapped to too many loci |       23143
            "nb_too_many": 0  # Number of reads mapped to too many loci |       23143
        }
        self.unmapped = {
            "nb": 0
        }
        self._parse()

    def _parse(self):
        with open(self.filepath) as FH_log:
            for line in FH_log:
                fields = [field.strip() for field in line.split("|")]
                if len(fields) == 2:
                    key = fields[0]
                    value = fields[1]
                    if key.startswith("Number of input reads"):
                        self.input["nb"] = int(value)
                    elif key.startswith("Average input read length"):
                        self.input["average_length"] = int(value)
                    elif key.startswith("Uniquely mapped reads number"):
                        self.unique["nb"] = int(value)
                    elif key.startswith("Number of reads mapped to multiple loci"):
                        self.multi["nb"] += int(value)
                    elif key.startswith("Number of reads mapped to too many loci"):
                        self.multi["nb"] += int(value)
        self.unmapped["nb"] = self.input["nb"] - self.unique["nb"] - self.multi["nb"]
