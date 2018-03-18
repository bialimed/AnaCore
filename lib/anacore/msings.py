#
# Copyright (C) 2018 IUCT-O
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
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

import gzip
from anacore.abstractFile import isGzip


class MsingsSample:
    def __init__(self, name, score, is_stable, loci=None):
        """
        @param name: [str] Name of the sample.
        @param score: [float] Score of the MSI.
        @param is_stable: [bool] True if the sample is stable.
        @param loci: [dict] The stability status by locus (with True the locus is stable).
        """
        self.name = name
        self.score = score
        self.is_stable = is_stable
        self.loci = dict() if loci is None else loci


    def getNbUnstable(self):
        """
        @summary: Returns the number of instable loci for the sample.
        @return: [int] The number of instable loci.
        """
        nb_unstable = 0
        for locus, is_stable in self.loci.items():
            if is_stable is not None and not is_stable:
                nb_unstable += 1
        return nb_unstable


    def getNbStable(self):
        """
        @summary: Returns the number of stable loci for the sample.
        @return: [int] The number of stable loci.
        """
        nb_stable = 0
        for locus, is_stable in self.loci.items():
            if is_stable is not None and is_stable:
                nb_stable += 1
        return nb_stable


    def getNbUsable(self):
        """
        @summary: Returns the number of loci usable in MSI evaluation for the sample.
        @return: [int] The number of instable loci.
        """
        nb_usable = 0
        for locus, is_stable in self.loci.items():
            if is_stable is not None:
                nb_usable += 1
        return nb_usable


    def getNbLoci(self):
        """
        @summary: Returns the number of loci in the sample.
        @return: [int] The number of loci.
        """
        return len(self.loci)


class CountMSI(object):
    """
    @summary: Manage output file produced by the command "msi count_msi_samples"
    of mSINGS (https://bitbucket.org/uwlabmed/msings).
    """

    def __init__(self, filepath):
        """
        @param filepath: [str] The filepath.
        """
        self.filepath = filepath
        self.samples = dict()
        self.loci = list()
        self.parse()


    def _parseFileHandle(self, FH):
        # Parse general information
        samples = [elt.strip() for elt in FH.readline().split("\t")][1:]
        nb_unstable = [elt.strip() for elt in FH.readline().split("\t")][1:]
        nb_evaluated = [elt.strip() for elt in FH.readline().split("\t")][1:]
        scores = [elt.strip() for elt in FH.readline().split("\t")][1:]
        for idx, elt in enumerate(scores):
            curr_score = None
            if elt != "":
                curr_score = float(elt)
            scores[idx] = curr_score
        status = [elt.strip() for elt in FH.readline().split("\t")][1:]
        for idx, elt in enumerate(status):
            is_stable = None
            if elt == "NEG":
                is_stable = True
            elif elt == "POS":
                is_stable = False
            status[idx] = is_stable
        for spl_idx, curr_spl in enumerate(samples):
            self.samples[curr_spl] = MsingsSample(curr_spl, scores[spl_idx], status[spl_idx])
        # Parse loci information
        for curr_line in FH:
            fields = [elt.strip() for elt in curr_line.split("\t")]
            curr_locus = fields[0]
            self.loci.append(curr_locus)
            for idx, curr_val in enumerate(fields[1:]):
                curr_spl = samples[idx]
                if curr_val == "":
                    self.samples[curr_spl].loci[curr_locus] = None
                elif curr_val == "1":
                    self.samples[curr_spl].loci[curr_locus] = False
                else:
                    self.samples[curr_spl].loci[curr_locus] = True


    def parse(self):
        if isGzip(self.filepath):
            with gzip.open(self.filepath, "rt") as FH:
                self._parseFileHandle(FH)
        else:
            with open(self.filepath, "r") as FH:
                self._parseFileHandle(FH)
