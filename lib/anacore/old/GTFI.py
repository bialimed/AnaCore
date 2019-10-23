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
__copyright__ = 'Copyright (C) 2014 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'


import re
from anacore.bioStructures import *

class GTFI(object):
    def __init__( self, filepath ):
        """
        @param filepath : [str] The filepath.
        """
        self.filepath = filepath
        self.current_line = None
        self.current_line_nb = None
        self.file_handle = None
        self.open()

    def __del__( self ):
        self.close()

    def __iter__( self ):
        for line in self.file_handle:
            line = line.rstrip()
            self.current_line = line
            self.current_line_nb += 1
            if line.startswith('#') :
                continue
            try:
                record = self._parse_line()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                                "Line content : " + self.current_line )
            else:
                yield record

    def open(self):
        self.file_handle = open( self.filepath )
        self.current_line = ""
        self.current_line_nb = 0

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line = None
            self.current_line_nb = None

    def load_model( self ):
        """
        @summary: Returns the structure model from the GTF.
        @return: [list] Two dict : transcripts and proteins. You can use transcript or protein ID to 
                 access at the corresponding object and her gene, exon, ....
        @see: bioStructure.
        """
        transcripts = dict()
        proteins = dict()
        prev_chromosome_name = ""
        chromosome = None
        current_gene = None

        self.open()
        for record in self:
            if prev_chromosome_name != record["region"]:
                chromosome = Reference( "", record["region"], "", "" )
                prev_chromosome_name = record["region"]
            if record["feature"] == "exon":
                # Transcript
                transcript = None
                if record["attr"]["transcript_id"] not in transcripts:
                    if current_gene is None:
                        gene_name = record["attr"]["gene_name"] if "gene_name" in record["attr"] else "Gene_noName"
                        current_gene = BiologicalRegion( -1, -1, record["strand"], chromosome, record["attr"]["gene_id"], gene_name )
                    transcript_name = record["attr"]["transcript_name"] if "transcript_name" in record["attr"] else "Transcript_noName"
                    transcript = Transcript( record["attr"]["transcript_id"], transcript_name, current_gene, record["strand"], list() )
                    transcripts[record["attr"]["transcript_id"]] = transcript
                else:
                    transcript = transcripts[record["attr"]["transcript_id"]]
                # Exon
                exon_id = record["attr"]["exon_id"] if "exon_id" in record["attr"] else record["attr"]["transcript_id"] + "_" + str(len(transcript.regions))
                exon = Exon( record["start"], record["end"], record["strand"], chromosome, exon_id, str(len(transcript.regions)) )
                transcript.addRegion( exon )
            elif record["feature"] == "CDS":
                # Transcript
                transcript = None
                if record["attr"]["transcript_id"] not in transcripts:
                    if current_gene is None:
                        gene_name = record["attr"]["gene_name"] if "gene_name" in record["attr"] else "Gene_noName"
                        current_gene = BiologicalRegion( -1, -1, record["strand"], chromosome, record["attr"]["gene_id"], gene_name )
                    transcript_name = record["attr"]["transcript_name"] if "transcript_name" in record["attr"] else "Transcript_noName"
                    transcript = Transcript( record["attr"]["transcript_id"], transcript_name, current_gene, record["strand"], list() )
                    transcripts[record["attr"]["transcript_id"]] = transcript
                else:
                    transcript = transcripts[record["attr"]["transcript_id"]]
                # Protein
                protein_id = record["attr"]["protein_id"] if "protein_id" in record["attr"] else "protein_NoId_" + record["attr"]["transcript_id"]
                protein = None
                if protein_id not in proteins:
                    protein = Protein( protein_id, "unknown", transcript, record["strand"], list() )
                    proteins[protein_id] = protein
                else:
                    protein = proteins[protein_id]
                # CDS
                cds = CDS( record["start"], record["end"], record["strand"], chromosome, protein_id + "_" + str(len(protein.regions)), str(len(protein.regions)) )
                protein.addRegion( cds )
            elif record["feature"] == "gene":
                gene_name = record["attr"]["gene_name"] if "gene_name" in record["attr"] else "Gene_noName"
                current_gene = BiologicalRegion( record["start"], record["end"], record["strand"], chromosome, record["attr"]["gene_id"], gene_name )
        self.close()

        return transcripts, proteins

    def _parse_line(self):
        """
        @summary : Returns a structured record from the GTF current line.
        @Returns : [dict] The record.
        """
        record = dict()
        record_attr = dict()
        record["region"], record["source"], record["feature"], record["start"], record["end"], record["score"], record["strand"], record["frame"], attributes = self.current_line.split("\t")
        record["start"] = int(record["start"])
        record["end"] = int(record["end"])
        attributes = attributes.strip()
        for current_attr in attributes.split(";"):
            current_attr = current_attr.strip()
            if current_attr != "" :
                matches = re.search("(.+)\s+\"([^\"]+)\"", current_attr)
                record_attr[matches.groups(1)[0]] = matches.groups(1)[1]
        record["attr"] = record_attr
        return record
