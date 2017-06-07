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
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'

import re
from copy import deepcopy
from vcf import *


    
class VEPVCFIO(VCFIO):
    def __init__(self, *args, **kwargs):
        self.CSQ_titles = list()
        super().__init__(*args, **kwargs)

    def _parse_header( self ):
        super()._parse_header()
        if "CSQ" in self.info:
            # Get CSQ fields
            match = re.search( "Format: ([^ ]+)", self.info["CSQ"]["description"] )
            if match is None:
                raise Exception("The CSQ description cannot be parsed in file " + self.filepath)
            self.CSQ_titles = match.group(1).split("|")            

    def _parse_line(self):
        record = super()._parse_line()
        self._parse_CSQ( record )
        return record

    def _parse_CSQ(self, record):
        for annot_idx, annot in enumerate(record.info["CSQ"]):
            annot_val = dict()
            for field_idx, field_value in enumerate(annot.split("|")):
                csq_key = self.CSQ_titles[field_idx]
                if field_value.strip() != "":
                    annot_val[csq_key] = field_value.strip()
                else:
                    annot_val[csq_key] = None
            record.info["CSQ"][annot_idx] = annot_val

    def recToVCFLine(self, record):
        std_record = deepcopy( record )
        if "CSQ" in record.info:
            csq_fields = list()
            for annot in record.info["CSQ"]:
                annot_fields = list()
                for csq_key in self.CSQ_titles:
                    if annot[csq_key] is None:
                        annot_fields.append("")
                    else:
                        annot_fields.append( str(annot[csq_key]) )
                csq_fields.append( "|".join(annot_fields) )
            std_record.info["CSQ"] = csq_fields
        return super().recToVCFLine(std_record)
