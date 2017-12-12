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

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class SamtoolsFlagStat (Component):

    def define_parameters(self, aln_files, search_dupl=False, optical_duplicate_pixel_distance=10,
                          read_name_regex=None, validation_stringency="LENIENT",
                          max_file_handles=40, sorting_collection_size_ratio=0.15):
        # Parameters
        self.add_parameter("search_dupl", "search_dupl", default=search_dupl, type=bool)
        self.add_parameter("optical_duplicate_pixel_distance", "optical_duplicate_pixel_distance", default=optical_duplicate_pixel_distance, type=int)
        self.add_parameter("read_name_regex", "read_name_regex", default=read_name_regex)
        self.add_parameter("validation_stringency", "validation_stringency", default=validation_stringency)
        self.add_parameter("max_file_handles", "max_file_handles", default=max_file_handles, type=int)
        self.add_parameter("sorting_collection_size_ratio", "sorting_collection_size_ratio", default=sorting_collection_size_ratio, type=float)

        # Files
        self.add_input_file_list("aln_files", "Path to the alignments files (format: BAM).", default=aln_files, required=True)
        self.add_output_file_list("flagstat_files", "flagstat stdout", pattern='{basename_woext}_flagstat.txt', items=self.aln_files)
        self.add_output_file_list("flagstat_stderr", "flagstat stderr", pattern='{basename_woext}_flagstat.stderr', items=self.aln_files)
        self.duplication_options = ""
        if self.search_dupl:
            self.add_output_file_list("duplication_files", "duplication metrics", pattern='{basename_woext}_dupl.txt', items=self.aln_files)
            self.add_output_file_list("duplication_stderr", "duplication stderr", pattern='{basename_woext}_dupl.stderr', items=self.aln_files)
            self.duplication_options += " ASSUME_SORTED=true"
            self.duplication_options += " MAX_FILE_HANDLES=" + str(self.max_file_handles)
            self.duplication_options += " SORTING_COLLECTION_SIZE_RATIO=" + str(self.sorting_collection_size_ratio)
            self.duplication_options += " VALIDATION_STRINGENCY=" + self.validation_stringency
            if self.read_name_regex != None:
                self.duplication_options += " READ_NAME_REGEX=" + self.read_name_regex
            if self.optical_duplicate_pixel_distance != None:
                self.duplication_options += " OPTICAL_DUPLICATE_PIXEL_DISTANCE=" + str(self.optical_duplicate_pixel_distance)


    def process(self):
        tmp_bam = self.aln_files

        # Duplication stats
        if self.search_dupl:
            tmp_bam = self.get_outputs('{basename_woext}_noDupl.bam', self.aln_files)
            cmd_mark = 'java -Xmx6g -jar ' + self.get_exec_path("picardtools") + ' MarkDuplicates INPUT=$1 METRICS_FILE=$2 OUTPUT=$3' + self.duplication_options + ' 2> $4'
            markduplicates = ShellFunction(cmd_mark, cmd_format='{EXE} {IN} {OUT}')
            MultiMap(markduplicates, inputs=self.aln_files, outputs=[self.duplication_files, tmp_bam, self.duplication_stderr])

        # Alignment summary
        cmd_flagstat = self.get_exec_path("samtools") + " view -F0x0100 -bh $1 | " + self.get_exec_path("samtools") + " flagstat - > $2 2> $3"
        flagstat = ShellFunction(cmd_flagstat, cmd_format='{EXE} {IN} {OUT}')
        MultiMap(flagstat, inputs=tmp_bam, outputs=[self.flagstat_files, self.flagstat_stderr])
