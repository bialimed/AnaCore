#!/usr/bin/env python3
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

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.bed import BEDIO
from anacore.region import Region, RegionList
from anacore.gff3 import GFF3IO


########################################################################
#
# FUNCTIONS
#
########################################################################
#~ name_by_acc = dict()
#~ with open(sys.argv[1]) as FH_names:
    #~ for line in FH_names:
        #~ if not line.startswith("#"):
            #~ fields = [field.strip() for field in line.split("\t")]
            #~ name_by_acc[fields[1]] = fields[0]

#~ with open(sys.argv[2]) as FH_names:
    #~ for line in FH_names:
        #~ if line.startswith("#"):
            #~ print(line.strip())
        #~ elif line.startswith("NC"):
            #~ fields = [field.strip() for field in line.split("\t")]
            #~ if fields[8].endswith(";"):
                #~ fields[8] += "chrAcc=" + fields[0] + ";"
            #~ else:
                #~ fields[8] += ";chrAcc=" + fields[0]
            #~ fields[0] = name_by_acc[fields[0]]
            #~ print("\t".join(fields))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="**********************************.")
    parser.add_argument('-c', '--is-thick-based', action='store_true', help='******************************.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-t', '--input-reference-tr', required=True, help='******************************.')
    group_input.add_argument('-a', '--input-annotation', required=True, help='*****************************.')
    group_input.add_argument('-r', '--input-regions', required=True, help='*****************************.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-regions', default="amplicons.bed", help='The amplicons description (format: BED). [Default: %(default)s]')
    args = parser.parse_args()

    # Get gene by transcript ID
    gene_by_tr = dict()
    with open(args.input_reference_tr) as FH_tr:
        for line in FH_tr:
            if not line.startswith("#"):
                gene, tr_id = [field.strip() for field in line.split("\t")]
                if tr_id != "":
                    tr_id = tr_id.split(".")[0]
                    gene_by_tr[tr_id] = gene

    # Get transcript region by transcript ID (exons are in annot attribute)
    tr_by_id = dict()
    with GFF3IO(args.input_annotation) as FH_annot:
        for record in FH_annot:
            if record.type == "mRNA" and "transcript_id" in record.attributes:
                tr_id = record.attributes["transcript_id"]
                tr_id = tr_id.split(".")[0]
                if tr_id in gene_by_tr:
                    if tr_id not in tr_by_id:
                        tr_by_id[tr_id] = Region(
                            record.start,
                            record.end,
                            record.strand,
                            record.seq_id,
                            tr_id,
                            {"exons": RegionList(), "gene": gene_by_tr[tr_id]}
                        )
            if record.type == "exon" and "transcript_id" in record.attributes:
                tr_id = record.attributes["transcript_id"]
                tr_id = tr_id.split(".")[0]
                if tr_id in gene_by_tr:
                    tr_by_id[tr_id].annot["exons"].append(
                        Region(
                            record.start,
                            record.end,
                            record.strand,
                            record.seq_id,
                            None,
                            {"exon_idx": None, "type":"exon"}
                        )
                    )
    if len(gene_by_tr) != len(tr_by_id):
        raise Exception("The following transcripts are missing in {}: {}".format(
            args.input_annotation,
            set(gene_by_tr.keys()).difference(set(tr_by_id.keys()))
        ))

    # Region list for chr
    tr_by_chr = dict()
    for tr_id in tr_by_id:
        curr_tr = tr_by_id[tr_id]
        # Sort exons by transcript order
        sorted_exons = list()
        if curr_tr.strand == "-":
            sorted_exons = sorted(curr_tr.annot["exons"], key=lambda exon: (exon.end, exon.start), reverse=True)
        else:
            sorted_exons = sorted(curr_tr.annot["exons"], key=lambda exon: (exon.start, exon.end))
        for idx_exon, curr_exon in enumerate(sorted_exons):
            curr_exon.annot["exon_idx"] = idx_exon
        # Store regions on chromosome
        chrom = curr_tr.reference.name
        if chrom not in tr_by_chr:
            tr_by_chr[chrom] = RegionList()
        tr_by_chr[chrom].append(curr_tr)

    # Write renamed regions
    with BEDIO(args.input_regions) as FH_regions:
        with BEDIO(args.output_regions, "w", (6 if args.is_thick_based else 8)) as FH_out:
            for record_idx, record in enumerate(FH_regions):
                amplicon = Region(record.start, record.end, record.strand, record.chrom)
                overlapped_tr = list()
                if record.chrom in tr_by_chr:
                    overlapped_tr = tr_by_chr[record.chrom].getOverlapped( amplicon )
                if len(overlapped_tr) >= 1:
                    overlapped_exons = overlapped_tr[0].annot["exons"].getOverlapped( amplicon )
                    features = list()
                    for curr_feature in overlapped_exons:
                        features.append( "ex" + str(curr_feature.annot["exon_idx"] + 1) )
                    record.name = "_".join([
                        overlapped_tr[0].annot["gene"],
                        "-".join(features),
                        str(amplicon.start),
                        amplicon.strand,
                        "ampl" + str(record_idx)
                    ])
                else:
                    record.name = "_".join([
                        "noFeature",
                        str(amplicon.start),
                        amplicon.strand,
                        "ampl" + str(record_idx)
                    ])
                if args.is_thick_based and record.thickStart is not None and record.thickEnd is not None:
                    record.start = record.thickStart
                    record.end = record.thickEnd
                FH_out.write( record )
