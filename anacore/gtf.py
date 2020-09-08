# -*- coding: utf-8 -*-
"""Classes and functions for reading/writing/processing GTF files."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.4'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
from anacore.abstractFile import AbstractFile
from anacore.region import RegionTree, RegionList
from anacore.genomicRegion import Gene, Transcript, Protein, Exon, CDS


class GTFIO(AbstractFile):
    """Class to manage read and write in GTF file."""

    def _parseLine(self):
        """
        Return a structured record from the GTF current line.

        :return: The record.
        :rtype: region.Region
        """
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
        prev = None
        for current_attr in attributes.split(";"):
            if current_attr.strip() != "":
                if prev is not None:  # Attribute value is split because it contains semicolon
                    current_attr = prev + current_attr
                matches = re.fullmatch("(.+)\s+\"([^\"]+)\"", current_attr.strip())
                if matches is None:
                    prev = current_attr + ";"
                else:
                    key, val = matches.groups(1)
                    if key != "tag":
                        record_attr[key] = val
                    else:  # "tag" contains a set of values so this key can can occurate several times (cf. www.genecodegenes/pages/data_format.html)
                        if key in record_attr:
                            record_attr[key].add(val)
                        else:
                            record_attr[key] = {val}
                    prev = None
        record.annot = record_attr
        return record

    def isRecordLine(self, line):
        """
        Return True if the line corresponds to a record (it is not a comment or an header line).

        :param line: The evaluated line.
        :type line: str
        :return: True if the line corresponds to a record.
        :rtype: bool
        """
        is_record = True
        if line.startswith("#"):
            is_record = False
        return is_record

    def recordToLine(self, record):
        """
        Return the record in GTF format.

        :param record: The region to process.
        :type record: region.Region
        :return: The GTF line.
        :rtype: str
        """
        attributes = []
        std_fields = {"source", "feature", "score", "frame"}
        for key, val in sorted(record.annot.items()):
            if key not in std_fields:
                if key != "tag" or not issubclass(val.__class__, (set, list, tuple)):
                    attributes.append('{} "{}"'.format(key, val))
                else:  # "tag" contains a set of values so this key can can occurate several times (cf. www.genecodegenes/pages/data_format.html)
                    for sub_val in sorted(val):
                        attributes.append('{} "{}"'.format(key, sub_val))
        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            ("" if record.reference is None else record.reference.name),
            ("." if "source" not in record.annot else record.annot["source"]),
            ("." if "feature" not in record.annot else record.annot["feature"]),
            record.start,
            record.end,
            ("." if "score" not in record.annot else record.annot["score"]),
            ("." if record.strand is None else record.strand),
            ("." if "frame" not in record.annot else record.annot["frame"]),
            "; ".join(attributes)
        )
        return line

    def write(self, record):
        """
        Write one line on GTF file.

        :param record: The object to write.
        :type record: region.Region
        """
        self.file_handle.write(self.recordToLine(record) + "\n")
        self.current_line_nb += 1


def _castedRegionTree(region, new_class):
    """
    Return the region casted in the selected class of genomicRegion.

    :param region: The region to cast.
    :type region: region.RegionTree
    :param new_class: The new class for the object.
    :type new_class: genomicRegion.*
    :return: The list of regions corresponding to the selected feature.
    :rtype: region.RegionList
    """
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


def loadModel(gtf_path, feature_handle="exons", restrict_to=None):
    """
    Return genomic model from a GTF. A genomic model is a tree where nodes are genes, transcripts, protein, exons and CDS.

    :param gtf_path: Path to the GTF containing the annotations.
    :type gtf_path: str
    :param feature_handle: The type of objects directly accessible in returned list. Authorized values: genes, exons, transcripts, proteins, cds.
    :type feature_handle: str
    :param restrict_to: If this parameter is used, the model is retricted to this seqname. Default: All the seqname of the GTF are used.
    :type restrict_to: str
    :return: The list of regions corresponding to the selected feature.
    :rtype: region.RegionList
    """
    if feature_handle not in ["genes", "exons", "transcripts", "proteins", "cds"]:
        raise Exception("The retrun type {} for loadModel is invalid.".format(feature_handle))
    genes = {}
    transcripts = {}
    proteins = {}
    exons = RegionList()
    cds = RegionList()
    # Build model from GTF
    with GTFIO(gtf_path) as FH_in:
        for record in FH_in:
            if restrict_to is None or restrict_to == record.reference.name:
                if record.annot["feature"] == "exon":
                    record = _castedRegionTree(record, Exon)
                    # Transcript
                    transcript_id = record.annot["transcript_id"]
                    transcript_uid = "{}:{}".format(record.reference.name, transcript_id)  # Add chromosome to prevent collision on sex chromosomes
                    if transcript_uid not in transcripts:
                        transcript_name = record.annot["transcript_name"] if "transcript_name" in record.annot else None
                        transcripts[transcript_uid] = Transcript(
                            None, None, record.strand, record.reference, transcript_name, {"feature": "transcript", "id": transcript_id}
                        )
                        # Parent gene
                        gene_id = record.annot["gene_id"]
                        gene_uid = "{}:{}".format(record.reference.name, gene_id)  # Add chromosome to prevent collision on sex chromosomes
                        if gene_uid not in genes:
                            gene_name = None
                            if "gene_name" in record.annot:
                                gene_name = record.annot["gene_name"]
                            elif "gene" in record.annot:
                                gene_name = record.annot["gene"]
                            gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                            genes[gene_uid] = gene
                        genes[gene_uid].addChild(transcripts[transcript_uid])
                    transcript = transcripts[transcript_uid]
                    # Exon
                    record.name = record.annot["transcript_id"] + "_e" + str(len(transcript.children))
                    record.annot["id"] = record.annot["exon_id"] if "exon_id" in record.annot else record.name
                    transcript.addChild(record)
                    exons.append(record)
                elif record.annot["feature"] == "CDS":
                    record = _castedRegionTree(record, CDS)
                    # Transcript
                    transcript_id = record.annot["transcript_id"]
                    transcript_uid = "{}:{}".format(record.reference.name, transcript_id)  # Add chromosome to prevent collision on sex chromosomes
                    if transcript_uid not in transcripts:
                        transcript_name = record.annot["transcript_name"] if "transcript_name" in record.annot else None
                        transcripts[transcript_uid] = Transcript(
                            None, None, record.strand, record.reference, transcript_name, {"feature": "transcript", "id": transcript_id}
                        )
                        # Parent gene
                        gene_id = record.annot["gene_id"]
                        gene_uid = "{}:{}".format(record.reference.name, gene_id)  # Add chromosome to prevent collision on sex chromosomes
                        if gene_uid not in genes:
                            gene_name = None
                            if "gene_name" in record.annot:
                                gene_name = record.annot["gene_name"]
                            elif "gene" in record.annot:
                                gene_name = record.annot["gene"]
                            gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                            genes[gene_uid] = gene
                        genes[gene_uid].addChild(transcripts[transcript_uid])
                    transcript = transcripts[transcript_uid]
                    # Protein
                    protein_id = record.annot["protein_id"] if "protein_id" in record.annot else "prot:None_tr:{}_gene:{}".format(record.annot["transcript_id"], record.annot["gene_id"])
                    protein_uid = "{}:{}".format(record.reference.name, protein_id)  # Add chromosome to prevent collision on sex chromosomes
                    protein = None
                    if protein_uid in proteins:
                        protein = proteins[protein_uid]
                    else:
                        protein = Protein(None, None, record.strand, record.reference, protein_id, {"feature": "protein", "id": protein_id}, None, None, transcript)
                        proteins[protein_uid] = protein
                    # CDS
                    protein.addChild(record)
                    cds.append(record)
                elif record.annot["feature"] == "gene":
                    gene_id = record.annot["gene_id"]
                    gene_uid = "{}:{}".format(record.reference.name, gene_id)  # Add chromosome to prevent collision on sex chromosomes
                    if gene_uid not in genes:
                        gene_name = None
                        if "gene_name" in record.annot:
                            gene_name = record.annot["gene_name"]
                        elif "gene" in record.annot:
                            gene_name = record.annot["gene"]
                        gene = Gene(None, None, record.strand, record.reference, gene_name, {"feature": "gene", "id": gene_id})
                        genes[gene_uid] = gene
    # Sort transcripts because reverse strand are added without exons information and coordinate are initially None
    for gene in genes.values():
        gene.sortChildren()
    # Select the right handler
    return_obj = None
    if feature_handle == "genes":
        return_obj = RegionList(genes.values())
    elif feature_handle == "exons":
        return_obj = exons
    elif feature_handle == "transcripts":
        return_obj = RegionList(transcripts.values())
    elif feature_handle == "proteins":
        return_obj = RegionList(proteins.values())
    elif feature_handle == "cds":
        return_obj = cds
    return return_obj
