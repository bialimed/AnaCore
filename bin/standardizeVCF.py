#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, VCFRecord, getAlleleRecord, HeaderInfoAttr
from anacore.annotVcf import AnnotVCFIO
from anacore.sequenceIO import IdxFastaIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def stdizeVCF(FH_ref, FH_in, FH_out, trace_unstandard=False, log=None):
    """
    Split alternatives alleles in multi-lines, removes unecessary reference and alternative nucleotids, move indel to most upstream position and update alt allele in annotations.

    :param FH_ref: File handle to the reference file (format: fasta with faidx).
    :type FH_ref: anacore.sequenceIO.IdxFastaIO
    :param FH_in: File handle to the variants file (format: VCF).
    :type FH_in: anacore.vcf.VCFIO
    :param FH_out: File handle to the standardized variants file (format: VCF).
    :type FH_out: anacore.vcf.VCFIO
    :param trace_unstandard: True if you want to keep the trace of the variant before standardization in INFO.
    :type trace_unstandard: bool
    :param log: Logger used.
    :type log: logging.Logger
    """
    nb_annot = {"exact": 0, "collocated": 0}
    is_annotated = issubclass(FH_out.__class__, AnnotVCFIO)
    # Header
    FH_out.copyHeader(FH_in)
    if trace_unstandard:
        FH_out.info["UNSTD"] = HeaderInfoAttr("UNSTD", type="String", number="1", description="The variant id (chromosome:position=reference/alternative) before standardization.")
    FH_out.writeHeader()
    # Records
    for record in FH_in:
        collocated_records = []
        for alt_idx, alt in enumerate(record.alt):
            alt_record = getAlleleRecord(FH_in, record, alt_idx)
            if trace_unstandard:
                alt_record.info["UNSTD"] = alt_record.getName()
            # Previous
            unstd = {"chrom": alt_record.chrom, "pos": alt_record.pos, "ref": alt_record.ref, "alt": alt_record.alt[0]}
            # Standardize pos, ref and alt
            alt_record.fastStandardize(FH_ref, 1000)
            # Update annotations
            if is_annotated and FH_in.annot_field in alt_record.info:
                cleaned_annot = []
                for idx_ann, annot in enumerate(alt_record.info[FH_in.annot_field]):
                    if unstd["alt"] == annot["Allele"]:
                        nb_annot["exact"] += 1
                        annot["Allele"] = alt_record.alt[0]
                        cleaned_annot.append(annot)
                    else:
                        nb_annot["collocated"] += 1
                alt_record.info[FH_in.annot_field] = cleaned_annot
            collocated_records.append(alt_record)
        if len(collocated_records) == 1:
            FH_out.write(collocated_records[0])
        else:
            for alt_record in sorted(collocated_records, key=lambda elt: (elt.refStart(), elt.refEnd())):  # Sorted splitted alleleles
                FH_out.write(alt_record)
    if log is not None and nb_annot["collocated"] != 0:
        log.warning("{}/{} annotations have been deleted because they concern collocated variant.".format(nb_annot["collocated"], nb_annot["exact"] + nb_annot["collocated"]))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Splits alternatives alleles in multi-lines and removes unecessary reference and alternative nucleotids, move indel to most upstream position and update alt allele in annotations.')
    parser.add_argument('-a', '--annotations-field', help='Field used to store annotations if the variants file is annotated.')
    parser.add_argument('-t', '--trace-unstandard', action='store_true', help='Use this option to add "UNSTD" tag in record INFO. This tag contains the trace of the variant before standardization: chromosome:position=reference/alternative.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_input.add_argument('-r', '--input-reference', required=True, help='The path to the reference used in variant calling to produced the input VCF (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the output file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with IdxFastaIO(args.input_reference, use_cache=True) as FH_ref:
        if args.annotations_field is not None:
            with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotations_field) as FH_out:
                with AnnotVCFIO(args.input_variants, annot_field=args.annotations_field) as FH_in:
                    stdizeVCF(FH_ref, FH_in, FH_out, args.trace_unstandard, log)
        else:
            with VCFIO(args.output_variants, "w") as FH_out:
                with VCFIO(args.input_variants) as FH_in:
                    stdizeVCF(FH_ref, FH_in, FH_out, args.trace_unstandard, log)

    log.info("End of job")
