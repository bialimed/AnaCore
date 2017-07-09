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
__copyright__ = 'Copyright (C) 2017 IUCT'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction

class ReadsStat (Component):
    
    def define_parameters(self, input_files, no_group=False, nb_threads=1):
        self.add_input_file_list( "input_files", "Paths to the input files (format: FASTQ).", default=input_files, required=True, file_format='fastq' )
        self.add_parameter( "no_group", "True disables grouping of bases for reads >50bp.", type='bool', default=no_group )
        self.add_parameter( "nb_threads", "Number of threads for fastQC", type='int', default=nb_threads )
        self.add_output_file_list( "stdouts", "Fastqc stdout files", pattern='{basename_woext}.stdout', items=input_files )
        self.add_output_file_list( "stderrs", "Fastqc stderr files", pattern='{basename_woext}.stderr', items=input_files )

    #~ def _parse_summary_file (self, summary_file):
        #~ """
        #~ @summary: Parse the summary file
        #~ @param summary_file: [str] the fastqc summary file
        #~ @return: [dict] {"basicstat" : "PASS", ...}
        #~ """
        #~ stats = {}
        #~ for line in open(summary_file, 'r').readlines():
            #~ parts = line.strip().split("\t")
            #~ if parts[1] == "Basic Statistics": stats["basicstat"] = parts[0]
            #~ if parts[1] == "Per base sequence quality": stats["perbasequal"] = parts[0]
            #~ if parts[1] == "Per sequence quality scores": stats["perseqqual"] = parts[0] 
            #~ if parts[1] == "Per base sequence content": stats["seqcontent"] = parts[0]
            #~ if parts[1] == "Per base GC content": stats["baseGC"] = parts[0]
            #~ if parts[1] == "Per sequence GC content": stats["seqGC"] = parts[0] 
            #~ if parts[1] == "Per base N content": stats["nperbase"] = parts[0]
            #~ if parts[1] == "Sequence Length Distribution": stats["seqlen"] = parts[0]
            #~ if parts[1] == "Sequence Duplication Levels": stats["seqduplication"] = parts[0] 
            #~ if parts[1] == "Overrepresented sequences": stats["seqoverrepresentation"] = parts[0]
            #~ if parts[1] == "Kmer Content": stats["kmers"] = parts[0]
            #~ if parts[1] == "Per tile sequence quality": stats["pertilequal"] = parts[0]
            #~ if parts[1] == "Adapter Content": stats["adaptercontent"] = parts[0]
        #~ return stats


    #~ def _parse_data_file (self, data_file):
        #~ """
        #~ Parse the data file
          #~ @param data_file : the fastqc data file
          #~ @return          : {"nbseq" : x, ...}
        #~ """
        #~ stats = {}
        #~ for line in open(data_file, 'r').readlines():
            #~ if line.startswith("Total Sequences"):
                #~ stats["nbseq"] = line.strip().split()[2]
            #~ if line.startswith("Sequence length"):
                #~ stats["seqlen"] = line.strip().split()[2]
            #~ if line.startswith("%GC"):
                #~ stats["gc"] = line.strip().split()[1]
        #~ return stats

    #~ def post_process(self):
        #~ results_files = []
        #~ metrics = []
        #~ for file in os.listdir(self.output_directory):
            #~ full_file_path = os.path.join(self.output_directory, file)
            #~ if file.endswith(".zip"):
                #~ results_files.append(full_file_path)
            #~ elif os.path.isdir(full_file_path):
                #~ sample = file[:-7]
                #~ summary_info = self.__parse_summary_file(os.path.join(full_file_path, "summary.txt"))
                #~ data_info = self.__parse_data_file(os.path.join(full_file_path, "fastqc_data.txt"))
                #~ self._add_result_element(sample, "nbseq", str(data_info["nbseq"]))
                #~ self._add_result_element(sample, "value", str(data_info["gc"]), "psgcpng")
                #~ self._add_result_element(sample, "value", str(data_info["seqlen"]), "sldpng")
                #~ self._add_result_element(sample, "seqoverrepresentation", str(summary_info["seqoverrepresentation"])) 
                #~ if "seqoverrepresentation" not in metrics : metrics.append("seqoverrepresentation")
                
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_tile_quality.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_tile_quality.png"), 
                                                                            #~ sample + ".per_tile_quality.png"), "ptqpng")
                    #~ self._add_result_element(sample, "result", str(summary_info["pertilequal"]), "ptqpng")
                    #~ if "pertilequal" not in metrics : metrics.append("pertilequal")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "adapter_content.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "adapter_content.png"), 
                                                                            #~ sample + ".adapter_content.png"), "acqpng")
                    #~ self._add_result_element(sample, "result", str(summary_info["adaptercontent"]), "acqpng")
                    #~ if "adaptercontent" not in metrics : metrics.append("adaptercontent")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_base_quality.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_base_quality.png"), 
                                                                            #~ sample + ".per_base_quality.png"), "pbqpng")
                    #~ self._add_result_element(sample, "result", str(summary_info["perbasequal"]), "pbqpng")
                    #~ if "perbasequal" not in metrics : metrics.append("perbasequal")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_base_gc_content.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_base_gc_content.png"), 
                                                                            #~ sample + ".per_base_gc_content.png"), "pbgcpng")
                    #~ self._add_result_element(sample, "result", str(summary_info["baseGC"]), "pbgcpng")
                    #~ if "baseGC" not in metrics : metrics.append("baseGC")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_base_n_content.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_base_n_content.png"), 
                                                                            #~ sample + ".per_base_n_content.png"), "pbnspng")  
                    #~ self._add_result_element(sample, "result", str(summary_info["nperbase"]), "pbnspng")
                    #~ if "nperbase" not in metrics : metrics.append("nperbase")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_base_sequence_content.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_base_sequence_content.png"), 
                                                                            #~ sample + ".per_base_sequence_content.png"), "pbspng")     
                    #~ self._add_result_element(sample, "result", str(summary_info["seqcontent"]), "pbspng")
                    #~ if "seqcontent" not in metrics : metrics.append("seqcontent")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_sequence_quality.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_sequence_quality.png"), 
                                                                            #~ sample + ".per_sequence_quality.png"), "psqpng")
                    #~ self._add_result_element(sample, "result", str(summary_info["perseqqual"]), "psqpng")
                    #~ if "perseqqual" not in metrics : metrics.append("perseqqual")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "per_sequence_gc_content.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "per_sequence_gc_content.png"), 
                                                                            #~ sample + ".per_sequence_gc_content.png"), "psgcpng") 
                    #~ self._add_result_element(sample, "result", str(summary_info["seqGC"]), "psgcpng")
                    #~ if "seqGC" not in metrics : metrics.append("seqGC")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "sequence_length_distribution.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "sequence_length_distribution.png"), 
                                                                            #~ sample + ".sequence_length_distribution.png"), "sldpng")     
                    #~ self._add_result_element(sample, "result", str(summary_info["seqlen"]), "sldpng")
                    #~ if "seqlen" not in metrics : metrics.append("seqlen")
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "duplication_levels.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "duplication_levels.png"), 
                                                                            #~ sample + ".duplication_levels.png"), "dlpng")     
                    #~ self._add_result_element(sample, "result", str(summary_info["seqduplication"]), "dlpng")  
                    #~ if "seqduplication" not in metrics : metrics.append("seqduplication")            
                #~ if os.path.isfile(os.path.join(full_file_path, "Images", "kmer_profiles.png")):
                    #~ self._add_result_element(sample, "img", self._save_file(os.path.join(full_file_path, "Images", "kmer_profiles.png"), 
                                                                            #~ sample + ".kmer_profiles.png"), "kppng")
                    #~ self._add_result_element(sample, "result", str(summary_info["kmers"]), "kppng")
                    #~ if "kmers" not in metrics : metrics.append("kmers")
        
        #~ for metric in metrics :
             #~ self._add_result_element("metrics", "metric", metric, metric)
        #~ # Finaly create and add the archive to the analysis
        #~ self._create_and_archive(results_files, self.archive_name)
    
    #~ def get_version(self):
        #~ cmd = [self.get_exec_path("fastqc"), "--version"]
        #~ p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        #~ stdout, stderr = p.communicate()
        #~ return stdout.split()[1]
    
    def process(self):
        cmd_fastqc = self.get_exec_path("fastqc") + \
                     ' --quiet' + \
                     ' --format fastq' + \
                     ' --threads ' + str(self.nb_threads) + \
                     ' --dir ' + self.output_directory + \
                     ' --outdir ' + self.output_directory + \
                     (' --nogroup' if self.no_group else '') + \
                     ' $1 > $2 2> $3'           ######################## --noextract
        fastqc = ShellFunction( cmd_fastqc, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( fastqc, inputs=self.input_files, outputs=[self.stdouts, self.stderrs] )
