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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse
from copy import deepcopy

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.bed import BEDIO
from anacore.sequenceIO import FastaIO
from anacore.region import Region, RegionList
from anacore.vcf import VCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def canBeMovedToInterest(overlapped_region, ref_seq, variant):
    """
    @summary: Returns True if the variant can be moved out of the primer toward the interest region. This movement keeps the same alternative sequence.
    @param overlapped_region: [Region] The primer. It must contains a location annotation: "upstream" if it is the first primer of the amplicon in reference strand '+' orientation and "downstream" otherwise.
    @param ref_seq: [str] The reference sequence.
    @param variant: [VCFRecord] The variant.
    @return: [bool] True if the variant can be moved out of the primer toward the interest region.
    """
    can_be_moved = False
    if overlapped_region.annot["location"] == "upstream":
        moved_variant = variant.getMostDownstream(ref_seq)
        if moved_variant.pos > overlapped_region.end:
            can_be_moved = True
    else:
        moved_variant = variant.getMostUpstream(ref_seq)
        if variant.isInsertion():
            if moved_variant.ref == ".":  # The standardization is ok
                if moved_variant.pos <= overlapped_region.start:  # Positions can be equals in case where the insertion end the interest area
                    can_be_moved = True
            else:  # The standardization has failed
                if moved_variant.pos + len(moved_variant.ref) - 1 < overlapped_region.start:
                    can_be_moved = True
        elif variant.isDeletion():
            if len(moved_variant.alt[0].replace(".", "")):  # The standardization is ok
                if moved_variant.pos + len(variant.ref) - 1 < overlapped_region.start:
                    can_be_moved = True
    return can_be_moved

def checkAmpliconsOverlap(in_regions):
    """
    @summary: Raises an exception if at least one of the region overlap an other.
    @param in_regions: [str] Path to the design (format: BED).
    """
    # Retrieve amplicons
    ampl_by_chr = dict()
    with BEDIO(in_regions) as FH_in:
        for record in FH_in:
            if record.chrom not in ampl_by_chr:
                ampl_by_chr[record.chrom] = list()
            ampl_by_chr[record.chrom].append({"start": record.start, "end": record.end})
    # Check overlaps
    for chrom in ampl_by_chr:
        prev_end = -1
        for ampl in sorted(ampl_by_chr[chrom], key=lambda x: (x["start"], x["end"])):
            if ampl["start"] <= prev_end:
                raise Exception('The amplicons defined in "' + in_regions + '" have an overlap between them.')
            prev_end = ampl["end"]

def getPrimersByChr(in_regions):
    """
    @summary: Returns the list of primers by chromosome.
    @param in_regions: [str] Path to the amplicons design with their primers (format: BED). The zone of interest is defined by thickStart and thickEnd.
    @return: [dict] By chromosome an instance of RegionList containing the primers. Each primer has an location annotation: upstream or downstream (this information is strand + based).
    """
    primers_by_chr = dict()
    with BEDIO(in_regions) as FH_in:
        for record in FH_in:
            if record.chrom not in primers_by_chr:
                primers_by_chr[record.chrom] = RegionList()
            if record.thickStart is None or record.thickEnd is None:
                raise Exception('The BED file "' + in_regions + '" does not contains thickStart and thickEnd for all the amplicons.')
            upstream_primer = Region(record.start, record.thickStart -1, record.strand, record.reference, None, {"location": "upstream"})
            primers_by_chr[record.chrom].append(upstream_primer)
            downstream_primer = Region(record.thickEnd + 1, record.end, record.strand, record.reference, None, {"location": "downstream"})
            primers_by_chr[record.chrom].append(downstream_primer)
    return primers_by_chr

def getSeqRecord(in_seq, selected_id):
    """
    @summary: Returns the selected sequence object from the sequences file.
    @param in_seq: [str] Path to the sequences file (format: fasta).
    @param selected_id: [str] The ID of the selected sequence.
    @return: [Sequence] The selected sequence object.
    """
    selected_record = None
    with FastaIO(in_seq) as FH_in:
        for record in FH_in:
            if record.id == selected_id:
                selected_record = record
    return selected_record

def getVariantRegion(variant):
    """
    @summary: Returns region object corresponding to the variant.
    @param variant: [VCFRecord] The variant.
    @return: [Region] The region object corresponding to the variant.
    @warnings: This function can only be used on variant with only one alternative allele.
    """
    std_variant = deepcopy(variant)
    std_variant.standardizeSingleAllele()
    return Region(
        std_variant.pos,
        std_variant.pos + len(std_variant.ref) - 1,  # Works also with standardized insertion
        None,
        std_variant.chrom
   )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Removes variants located on amplicons primers. The variants defined on primers but with possibility to be move on zone of interest with the same consequences are kept.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF). This file should be sorted by coordinates otherwise the execution time will be dramatically increased.')
    group_input.add_argument('-r', '--input-regions', required=True, help='Path to the amplicons design with their primers (format: BED). The zone of interest is defined by thickStart and thickEnd. The amplicons must not have any overlap between them.')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to the reference sequences file (format: fasta). The reference used to discover variants.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', default="filtered.vcf", help='Path to the outputted variants file (format: VCF). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    checkAmpliconsOverlap(args.input_regions)
    primers_by_chr = getPrimersByChr(args.input_regions)
    chr_seq = None
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.filter["PRIM"] = 'The variant is located on an amplicon primer (amplicon desgin: ' + args.input_regions + ').'
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                for alt_idx, alt in enumerate(record.alt):
                    is_kept = True
                    alt_record = getAlleleRecord(FH_in, record, alt_idx)
                    alt_region = getVariantRegion(alt_record)
                    overlapped_primers = primers_by_chr[alt_record.chrom].getOverlapped(alt_region)
                    if len(overlapped_primers) > 0:  # The variant overlaps a primer or variant is an insertion just before the downstream primer
                        is_kept = False
                        if len(overlapped_primers) == 1:  # Variants over 2 primers are removed
                            if alt_record.isIndel():
                                if chr_seq is None or chr_seq.id != alt_record.chrom:
                                    chr_seq = getSeqRecord(args.input_sequences, alt_record.chrom)
                                if canBeMovedToInterest(overlapped_primers[0], chr_seq.string, alt_record):
                                    is_kept = True
                    if is_kept:
                        FH_out.write(alt_record)
