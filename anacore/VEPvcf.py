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
__version__ = '1.2.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
from copy import deepcopy
from anacore.vcf import *



class VEPVCFIO(VCFIO):
    def __init__(self, *args, **kwargs):
        self.CSQ_titles = list()
        super().__init__(*args, **kwargs)

    def _parseHeader(self):
        super()._parseHeader()
        if "CSQ" in self.info:
            # Get CSQ fields
            match = re.search("Format: ([^ ]+)", self.info["CSQ"]["description"])
            if match is None:
                raise Exception("The CSQ description cannot be parsed in file " + self.filepath)
            self.CSQ_titles = match.group(1).split("|")

    def _parseLine(self):
        record = super()._parseLine()
        self._parseCSQ(record)
        return record

    def _parseCSQ(self, record):
        if "CSQ" in self.info:
            if "CSQ" not in record.info or record.info["CSQ"] is None:  # CSQ is empty
                record.info["CSQ"] = list()
            else:  # The variant has consequence(s)
                for annot_idx, annot in enumerate(record.info["CSQ"]):
                    annot_val = dict()
                    for field_idx, field_value in enumerate(annot.split("|")):
                        csq_key = self.CSQ_titles[field_idx]
                        if field_value.strip() != "":
                            annot_val[csq_key] = field_value.strip()
                        else:
                            annot_val[csq_key] = None
                    record.info["CSQ"][annot_idx] = annot_val

    def copyHeader(self, model):
        """
        @summary: Copy header fields from the specified VCF.
        @param model: [VEPVCFIO] The VCF source.
        """
        super().copyHeader(model)
        self.CSQ_titles = model.CSQ_titles

    def recToVCFLine(self, record):
        std_record = deepcopy(record)
        if "CSQ" in record.info:
            if record.info["CSQ"] is None or len(record.info["CSQ"]) == 0:  # CSQ is empty
                std_record.info.pop('CSQ', None)
            else:  # The variant has consequence(s)
                csq_fields = list()
                for annot in record.info["CSQ"]:
                    annot_fields = list()
                    for csq_key in self.CSQ_titles:
                        if annot[csq_key] is None:
                            annot_fields.append("")
                        else:
                            annot_fields.append(str(annot[csq_key]))
                    csq_fields.append("|".join(annot_fields))
                std_record.info["CSQ"] = csq_fields
        return super().recToVCFLine(std_record)
