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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
from anacore.abstractFile import AbstractFile
from anacore.region import RegionTree, RegionList
from anacore.genomicRegion import Gene, Transcript, Protein, Exon, CDS


class GTFIO(AbstractFile):
    def _parseLine(self):
        chrom, source, feature, start, end, score, strand, frame, attributes = [elt.strip() for elt in self.current_line.split('\t')]
        # Core
        record = RegionTree(
            start,
            end,
            strand,
            chrom
        )
        # Attributes
        record_attr = {
            "source": source,
            "feature": feature,
            "score": score,
            "frame": frame
        }
        attributes = attributes.strip()
        for current_attr in attributes.split(";"):
            current_attr = current_attr.strip()
            if current_attr != "":
                matches = re.search("(.+)\s+\"([^\"]+)\"", current_attr)
                record_attr[matches.groups(1)[0]] = matches.groups(1)[1]
        record.annot = record_attr
        return record

    def isRecordLine(self, line):
        is_record = True
        if line.startswith("#"):
            is_record = False
        return is_record

    #~ def write(self, record):
        #~ self.file_handle.write(self.RecordToGTFLine(bed_record) + "\n")
        #~ self.current_line_nb += 1

def _castedRegionTree(region, new_class):
    casted_obj = new_class(
        start=region.start,
        end=region.end,
        strand=region.strand,
        reference=region.reference,
        name=region.name,
        annot=region.annot,
        parent=region.parent,
        children=region.children
    )
    return casted_obj


def loadModel(gtf_path, rtype="exons"):
    if rtype not in ["genes", "exons", "transcripts", "proteins", "cds"]:
        raise Exception("The retrun type {} for loadModel is invalid.".format(rtype))
    genes = {}
    transcripts = {}
    proteins = {}
    exons = RegionList()
    cds = RegionList()

    curr_gene = None
    with GTFIO(gtf_path) as FH_in:
        for record in FH_in:
            if record.annot["feature"] == "exon":
                record = _castedRegionTree(record, Exon)
                # Transcript
                transcript_id = record.annot["transcript_id"]
                if transcript_id not in transcripts:
                    transcript_name = record.annot["transcript_name"] if "transcript_name" in record.annot else None
                    transcripts[transcript_id] = Transcript(
                        None, None, record.strand, record.reference, transcript_name, {"feature": "transcript", "id": transcript_id}
                    )
                    if curr_gene is None:
                        gene_id = record.annot["gene_id"]
                        if gene_id not in genes:
                            gene_name = record.annot["gene_name"] if "gene_name" in record.annot else None
                            gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                            genes[gene_id] = gene
                        curr_gene = genes[gene_id]
                    curr_gene.addChild(transcripts[transcript_id])
                transcript = transcripts[transcript_id]
                # Exon
                record.name = record.annot["transcript_id"] + "_e" + str(len(transcript.children))
                record.annot["id"] = record.annot["exon_id"] if "exon_id" in record.annot else record.name
                transcript.addChild(record)
                exons.append(record)
            elif record.annot["feature"] == "CDS":
                record = _castedRegionTree(record, CDS)
                # Transcript
                transcript_id = record.annot["transcript_id"]
                if transcript_id not in transcripts:
                    transcript_name = record.annot["transcript_name"] if "transcript_name" in record.annot else None
                    transcripts[transcript_id] = Transcript(
                        None, None, record.strand, record.reference, transcript_name, {"feature": "transcript", "id": transcript_id}
                    )
                    if curr_gene is None:
                        gene_id = record.annot["gene_id"]
                        if gene_id not in genes:
                            gene_name = record.annot["gene_name"] if "gene_name" in record.annot else None
                            gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                            genes[gene_id] = gene
                        curr_gene = genes[gene_id]
                    curr_gene.addChild(transcripts[transcript_id])
                transcript = transcripts[transcript_id]
                # Protein
                protein_id = record.annot["protein_id"] if "protein_id" in record.annot else "prot:None_tr:" + record.annot["transcript_id"]
                protein = None
                if protein_id not in proteins:
                    protein = Protein(None, None, record.strand, record.reference, protein_id, {"feature": "protein", "id": protein_id}, None, None, transcript)
                    proteins[protein_id] = protein
                else:
                    protein = proteins[protein_id]
                # CDS
                protein.addChild(record)
                cds.append(record)
            elif record.annot["feature"] == "gene":
                gene_id = record.annot["gene_id"]
                if gene_id not in genes:
                    gene_name = record.annot["gene_name"] if "gene_name" in record.annot else None
                    gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                    genes[gene_id] = gene
                curr_gene = genes[gene_id]
    # Sort transcripts because reverse strand are added without exons information and coordinate are initially None
    for gene in genes.values():
        gene.sortChildren()
    # Select the right handler
    return_obj = None
    if rtype == "genes":
        return_obj = RegionList(genes.values())
    elif rtype == "exons":
        return_obj = exons
    elif rtype == "transcripts":
        return_obj = RegionList(transcripts.values())
    elif rtype == "proteins":
        return_obj = RegionList(proteins.values())
    elif rtype == "cds":
        return_obj = cds
    return return_obj
