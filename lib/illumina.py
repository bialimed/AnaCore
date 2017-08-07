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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


class SampleSheetIO(object):
    def __init__(self, path ):
        self.filepath = path
        self.run = None
        self.samples = None
        self.header = None
        self.manifests = None
        self._parse()

    def _parse(self):
        # Retrieve information by section
        with open( self.filepath ) as FH_sheet:
            sections_by_title = dict()
            section_title = None
            for line in FH_sheet:
                if line.strip() in ["[Header]", "[Manifests]", "[Reads]", "[Settings]", "[Data]"]:
                    section_title = line.strip()[1:-1]
                    sections_by_title[section_title] = list()
                elif line.strip() == "":
                    section_title = None
                else:
                    sections_by_title[section_title].append( line.strip() )
        # Process information
        self.samples = self._getSamplesFromData( sections_by_title["Data"] )
        self.header = self._getInfoFromSection( sections_by_title["Header"] )
        self.manifests = self._getInfoFromSection(sections_by_title["Manifests"]) if "Manifests" in sections_by_title else dict()

    def _getSamplesFromData(self, data_section):
        samples = list()
        data_titles = [field.strip() for field in data_section[0].split(",")]
        for line in data_section[1:]:
            samples.append( {data_titles[idx]:field.strip() for idx, field in enumerate(line.split(","))} )
        return( samples )

    def _getInfoFromSection(self, section):
        info = dict()
        for line in section:
            key, value = [field.strip() for field in line.split(",", 1)]
            info[key] = value
        return( info )
