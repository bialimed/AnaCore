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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.region import Region
from anacore.annotVcf import AnnotVCFIO, VCFIO, HeaderFilterAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
# def strandsAreOpposite(alt):
#     return True if alt.startswith("[") or alt.endswith("]") else False  # Shard is in front of the open of bracket


def shardStrand(record):
    """
    Return the strand of the shard used in the fusion.

    :param record: A breakend record.
    :type record: anacore.vcf.VCFRecord
    :return: Strand of the shard ('+' or '-').
    :rtype: str
    """
    strand = None
    alt_starts = set([alt[0] for alt in record.alt])
    if len(alt_starts) > 1:
        raise Exception("The alternatives for the record {} have opposites configurations: {}.".format(record.id, record.alt))
    if "RNA_FIRST" in record.info:  # First shard
        if record.alt[0].startswith("[") or record.alt[0].startswith("]"):
            strand = "-"
        else:
            strand = "+"
    else:  # Second shard
        if record.alt[0].startswith("[") or record.alt[0].startswith("]"):
            strand = "+"
        else:
            strand = "-"
    return strand


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter fusions do not merge the upstream of a gene with the downstream of another gene with a strand corresponding of the reading strand of each partner.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the file containing annotated variants. (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the filtered file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Filter breakends
    log.info("Filter breakends.")
    kept_fusions = 0
    nb_fusions = 0
    with AnnotVCFIO(args.output_variants, "w") as FH_out:
        with AnnotVCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.filter["not_expr_model"] = HeaderFilterAttr("not_expr_model", "Filter fusions do not merge the upstream of a gene with the downstream of another gene with a strand corresponding of the reading strand of each partner.")
            FH_out.writeHeader()
            # Records
            mate_by_id = {}
            for bnd in FH_in:
                if bnd.info["SVTYPE"] == "BND":
                    if len(bnd.alt) > 1:
                        log.error("Breakend {} contains several alternatives and this configuration is not implemented.".format(bnd.id))
                    for mate_bnd_id in bnd.info["MATEID"]:
                        if mate_bnd_id not in mate_by_id:
                            mate_by_id[bnd.id] = bnd
                        else:
                            nb_fusions += 1
                            mate_bnd = mate_by_id[mate_bnd_id]
                            if args.annotation_field in bnd.info and args.annotation_field in mate_bnd.info:  # The two shards contain a gene
                                kept_ann_first_shard = set()
                                kept_ann_second_shard = set()
                                # Shards in transcript order
                                first_shard = bnd
                                second_shard = mate_bnd
                                if "RNA_FIRST" in mate_bnd.info:
                                    first_shard = mate_bnd
                                    second_shard = bnd
                                # Find couple of annotations compatible with the expression model
                                for idx_ann_first, ann_first in enumerate(first_shard.info[args.annotation_field]):
                                    if ann_first["GENE_SHARD"] is not None and ann_first["GENE_SHARD"] == "up":  # The first shard on transcript contains the start of a gene
                                        ann_first_strand = "+" if ann_first["STRAND"] == "1" else "-"
                                        if shardStrand(first_shard) == ann_first_strand:  # The shard is transcripted in gene strand
                                            for idx_ann_second, ann_second in enumerate(second_shard.info[args.annotation_field]):
                                                if ann_first["Gene"] != ann_second["Gene"]:  # The two genes implicated in fusion must be differents
                                                    if ann_second["GENE_SHARD"] is not None and ann_second["GENE_SHARD"] == "down":  # The second shard on transcript contains the downtream section of a gene
                                                        ann_second_strand = "+" if ann_second["STRAND"] == "1" else "-"
                                                        if shardStrand(second_shard) == ann_second_strand:  # The shard is transcripted in gene strand
                                                            kept_ann_first_shard.add(idx_ann_first)
                                                            kept_ann_second_shard.add(idx_ann_second)
                                # Delete invalid annot and write breakends
                                if len(kept_ann_first_shard) != 0 and len(kept_ann_second_shard) != 0:
                                    first_shard.info[args.annotation_field] = [first_shard.info[args.annotation_field][idx] for idx in kept_ann_first_shard]
                                    second_shard.info[args.annotation_field] = [second_shard.info[args.annotation_field][idx] for idx in kept_ann_second_shard]
                                    FH_out.write(bnd)
                                    FH_out.write(mate_bnd)
                                    kept_fusions += 1
    log.info("Kept {}/{} ({:.1f}%) fusions.".format(kept_fusions, nb_fusions, 100 * kept_fusions / nb_fusions))
    log.info("End of job")
