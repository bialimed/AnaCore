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
import time
import random
import logging
import argparse
import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
LIB_DIR = os.path.abspath(os.path.join(APP_DIR, "lib"))
sys.path.append(LIB_DIR)

from anacore.sequenceIO import FastaIO, Sequence
from anacore.sv import HashedSVIO
from anacore.bed import getAreasByChr
from anacore.region import Region


########################################################################
#
# FUNCTIONS
#
########################################################################
def revcom(seq):
    """
    Return the reverse complement the nucleic sequence.

    :param seq: The sequence.
    :type seq: str
    :return: The reverse complement of the sequence.
    :rtype: str
    """
    complement_rules = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N',
                        'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a', 'n': 'n'}

    return("".join([complement_rules[base] for base in seq[::-1]]))


def getVariantsProfile(profile_path, min_allele_freq=None):
    profiles = list()
    with HashedSVIO(profile_path) as FH_profile:
        for record in FH_profile:  # Type   Occurence Freq_min    Freq_max    Lg_min  Lg_max
            curr_profile = {
                "type": record["Type"],
                "occurence": float(record["Occurence"]),
                "AF": {"min": float(record["Freq_min"]), "max": float(record["Freq_max"])},
                "length": {"min": int(record["Lg_min"]), "max": int(record["Lg_max"])}
            }
            profiles.append(curr_profile)
            if min_allele_freq is not None:
                if min_allele_freq > curr_profile["AF"]["min"]:
                    log.error("The minimum allele frequency in {} must be >= {}.".format(profile_path, min_allele_freq))
                if curr_profile["AF"]["min"] != round(curr_profile["AF"]["min"], int(1 / min_allele_freq)):
                    log.error("The allele frequency precision must be >= {}.".format(min_allele_freq))
                if curr_profile["AF"]["max"] != round(curr_profile["AF"]["max"], int(1 / min_allele_freq)):
                    log.error("The allele frequency precision must be >= {}.".format(min_allele_freq))
    return profiles


def getCoveredList(areas_by_chr):
    """
    must be non overlapping
    """
    covered_pos = list()
    for chrom_id in sorted(areas_by_chr):
        for curr_area in sorted(areas_by_chr[chrom_id], key=lambda x: x.start):
            for pos in range(curr_area.start, curr_area.end + 1):
                covered_pos.append("{}:{}".format(chrom_id, pos))
    return(covered_pos)


def setAlt(variant, region_sequence, start_idx):
    """
    """
    variant["ref"] = ""
    variant["alt"] = ""
    if variant["type"] == "substitution":  # Substitution
        for idx in range(variant["length"]):
            # Ref
            variant["ref"] += region_sequence[start_idx + idx]
            # Alt
            choice_opt = ["A", "T", "G", "C"]
            choice_opt.remove(variant["ref"][-1].upper())
            variant["alt"] += random.choice(choice_opt)
    elif variant["type"] == "deletion":  # Deletion
        for idx in range(variant["length"]):
            variant["ref"] += region_sequence[start_idx + idx]
    else:  # Insertion
        variant["ref"] = region_sequence[start_idx]
        variant["alt"] = variant["ref"]
        for idx in range(variant["length"]):
            variant["alt"] += random.choice(["A", "T", "G", "C"])


def writeVariants(out_path, variants_by_chr):
    """
    """
    with open(out_path, "w") as FH_out_vcf:
        FH_out_vcf.write('##fileformat=VCFv4.2\n')
        FH_out_vcf.write('##source=simuReads\n')
        FH_out_vcf.write('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency in simulation">\n')
        FH_out_vcf.write('##INFO=<ID=EAF,Number=.,Type=Float,Description="Allele Frequency expected by profile (before simulation optimisation)">\n')
        FH_out_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
        for chrom_id, variants_by_start in variants_by_chr.items():
            for start in sorted(variants_by_start):
                variant = variants_by_start[start]
                FH_out_vcf.write(
                    "{}\t{}\t.\t{}\t{}\t.\t.\tAF={};EAF={}\t.\n".format(
                        chrom_id,
                        start,
                        variant["ref"],
                        variant["alt"],
                        variant["simu_AF"],
                        variant["AF"]
                    )
                )


def selectMatchedVariants(region, variant_by_pos):
    """
    """
    matched_variants = []
    for start, variant in variant_by_pos.items():
        if start >= region.start and start <= region.end:  # Variants starts on region
            matched_variants.append(variant)
        elif variant["type"] == "substitution" or variant["type"] == "deletion":
            variant_end = variant["start"] + variant["length"] - 1
            if variant_end >= region.start and variant_end <= region.end:  # Deletion not starts on the region but ends on it (overlapping two targets separated by non targeted)
                log.debug("The variant {} overlap 2 targets ; it ends on {}.".format(variant["id"], region))
                # matched_variants.append(variant)  ################################################################################## pb no management by the update of getTargetReads
    return sorted(matched_variants, key=lambda x: x["start"])


def applyVariants(region, ref_sequence, alt_records, variant_by_pos):
    """
    """
    nb_alt_sequences = len(alt_records)
    for pos in range(region.start, region.end + 1):
        if pos in variant_by_pos:
            curr_var = variant_by_pos[pos]
            # Set alternative variant
            setAlt(curr_var, ref_sequence, pos - region.start)
            alt_allele = curr_var["alt"]
            ref_allele = curr_var["ref"]
            # Random selection of alternative sequence
            nb_alterated_seq = int(curr_var["AF"] * nb_alt_sequences)
            # Apply variant in selected alternatives
            for curr_alt in random.sample(alt_records, nb_alterated_seq):
                curr_alt_seq = curr_alt["seq"]
                curr_alt["mut"][pos] = curr_var  # Possible because it does not exist any overlap between variants
                # Apply variant
                if "-" in alt_allele or len(alt_allele) < len(ref_allele):  # Deletion
                    del_length = len(ref_allele)
                    for idx in range(del_length):
                        curr_alt_seq[pos - region.start + idx] = ""
                elif "-" in ref_allele or len(alt_allele) > len(ref_allele):  # Insertion
                    curr_alt_seq[pos - region.start] = alt_allele
                else:  # Substitution
                    subst_length = len(alt_allele)
                    for idx in range(subst_length):
                        curr_alt_seq[pos - region.start + idx] = alt_allele[idx]


def mergeOverlapped(regions, padding=0, trace=False):
    """
    """
    sorted_regions = sorted(regions, key=lambda x: (x.start, x.end))
    deleted_idx = []
    prev_region = Region(-1, -1)
    # Extend regions
    for idx, curr_region in enumerate(sorted_regions):
        curr_start = max(1, curr_region.start - padding)
        prev_end = curr_region.end + padding
        if curr_start <= prev_end:  # Overlap between regions
            if trace:
                if "merge_traceback" not in prev_region.annot:
                    prev_region.annot["merge_traceback"] = [Region(
                        prev_region.start,
                        prev_region.end,
                        prev_region.strand,
                        prev_region.reference,
                        prev_region.name
                    )]
                prev_region.annot["merge_traceback"].append(Region(
                    curr_region.start,
                    curr_region.end,
                    curr_region.strand,
                    curr_region.reference,
                    curr_region.name
                ))
            prev_region.end = max(curr_region.end, prev_region.end)  # Max to manage included regions
            deleted_idx.append(idx)
        else:
            prev_region = curr_region
    # Delete useless regions
    for idx in sorted(deleted_idx, reverse=True):
        del(sorted_regions[idx])


def mergeOverlappedWithPadding(regions_by_chr, padding=0):
    """
    """
    for chrom_id, regions in regions_by_chr.items():
        mergeOverlapped(regions, padding, True)


def getChrSeq(in_ref, chrom_id):
    """
    """
    chrom_seq = None
    with FastaIO(in_ref) as FH_seq:
        for record in FH_seq:
            if record.id == chrom_id:
                chrom_seq = record.string
    return chrom_seq


def setReadsPair(fragment_record, reads_length, end_R1_adpter, end_R2_adapter):
    fragment_length = len(fragment_record.annot["seq"])
    reads_id = "frag={}_coord={}:{}-{}_strand={}".format(
        fragment_record.name,
        fragment_record.reference.name,
        fragment_record.start,
        fragment_record.end,
        fragment_record.strand
    )
    # Reads from matrix
    R1 = Sequence(reads_id, fragment_record.annot["seq"][:reads_length])
    R2 = Sequence(reads_id, fragment_record.annot["seq"][-reads_length:])
    R2.string = revcom(R2.string)
    if fragment_record.strand == "-":
        R1_old = R1
        R1 = R2
        R2 = R1_old
    fragment_record.annot["reads"] = (R1, R2)
    fragment_record.annot["read_up"] = {
        "start": fragment_record.start,
        "end": fragment_record.start + reads_length
    }
    fragment_record.annot["read_down"] = {
        "start": fragment_record.end - reads_length,
        "end": fragment_record.end
    }
    # Add adapter for fragment < reads length
    if reads_length > fragment_length:
        fragment_record.annot["read_up"]["end"] = fragment_record.start + len(R1.string)
        fragment_record.annot["read_down"]["start"] = fragment_record.end - len(R1.string)
        R1.string += end_R1_adpter[:reads_length - fragment_length]
        R2.string += end_R2_adapter[:reads_length - fragment_length]


def getPartialFragment(region_seq, start_idx_in_reg, fragment_len):
    """
    seq, end_idx_in_reg, missing_len
    """
    idx_in_reg = start_idx_in_reg
    curr_len_fragment = 0
    len_region = len(region_seq)
    while idx_in_reg < len_region and curr_len_fragment < fragment_len:
        nt_on_pos = region_seq[idx_in_reg]  ################################################## pb perf
        curr_len_fragment += len(nt_on_pos)
        idx_in_reg += 1
    end_idx_in_reg = idx_in_reg - 1
    missing = fragment_len - curr_len_fragment
    fragment_seq = "".join(region_seq[start_idx_in_reg:end_idx_in_reg + 1])
    if missing == 0:
        fragment_seq = fragment_seq[:fragment_len]  # If the fragment end with an insertion
    return fragment_seq, end_idx_in_reg, missing


def getFragmentRegion(chrom_seq, target, target_seq, start_pos, fragment_len):
    fragment_seq = ""
    end_pos = None
    if start_pos > target.end:  # Fragment starts after target
        end_pos = start_pos + fragment_len - 1
        fragment_seq = chrom_seq[start_pos - 1:end_pos]  # Position is 1-based indexes are 0-based
    elif start_pos + fragment_len - 1 < target.start:  # Fragment ends before target
        end_pos = start_pos + fragment_len - 1
        fragment_seq = chrom_seq[start_pos - 1:end_pos]  # Position is 1-based indexes are 0-based
    else:  # Fragment overlap target
        start_idx_on_target = start_pos - target.start
        # Before target
        if start_pos < target.start:  # Fragment starts before target
            start_idx_on_target = 0
            add_start_pos = start_pos
            add_end_pos = target.start - 1
            fragment_seq = chrom_seq[add_start_pos - 1:add_end_pos]  # Position is 1-based indexes are 0-based
        # On target
        fragment_seq_on_target, end_idx, missing_len = getPartialFragment(
            target_seq,
            start_idx_on_target,
            fragment_len - len(fragment_seq)
        )
        fragment_seq += fragment_seq_on_target
        end_pos = target.start + end_idx + missing_len
        # After target
        if missing_len > 0:  # Fragment ends after target
            add_start_pos = target.end + 1
            add_end_pos = add_start_pos + missing_len - 1
            fragment_seq += chrom_seq[add_start_pos - 1:add_end_pos]  # Position is 1-based indexes are 0-based
    return Region(
        start_pos,
        end_pos,
        None,
        target.reference,
        None,
        {"seq": fragment_seq}
    )


def updateOverlapping(overlapping, fragment):
    for curr_read in ["up", "down"]:
        for curr_pos in range(fragment.annot["read_up"]["start"], fragment.annot["read_up"]["end"] + 1):
            if curr_pos in fragment.annot["mut"]:
                overlapping[curr_pos].append((
                    curr_read,
                    fragment.annot["mut"][curr_pos]["id"]
                ))
            else:
                overlapping[curr_pos].append((curr_read, None))


def writeTargetReads(out_R1_path, out_R2_path, reads_pairs):
    """
    """
    with FastaIO(out_R1_path, "a") as FH_out_R1:
        with FastaIO(out_R2_path, "a") as FH_out_R2:
            for R1, R2 in reads_pairs:
                FH_out_R1.write(R1)
                FH_out_R2.write(R2)


def simulationPenaltyScore(res_simulated, res_expected):
    penalty_score = 0
    if len(res_expected):
        nb_lost = res_simulated.count(0)
        ratio_lost = nb_lost / len(res_expected)
        nb_out_of_threshold = 0
        error_rates = []
        for expected, simulated in zip(res_expected, res_simulated):
            rate = min(expected, simulated) / max(expected, simulated)
            error_rates.append(1 - rate)
            if rate < 0.95:
                nb_out_of_threshold += 1
        ratio_out_of_05 = nb_out_of_threshold / len(res_expected)
        error_rates_median = np.median(error_rates)
        penalty_score = 0.99 * ratio_lost + 0.008 * ratio_out_of_05 + 0.002 * error_rates_median
        log.debug('Simulation nb_lost: {}, ratio_out_of_05: {}, error_rates_median: {}, penalty_score: {}'.format(nb_lost, ratio_out_of_05, error_rates_median, penalty_score))
    return penalty_score


def getTargetReads(chrom_seq, chrom_len, target, alt_records, fragments_len, args):
    padded_start = max(1, target.start - args.targets_padding)
    padded_end = min(chrom_len, target.end + args.targets_padding)
    large_padded_start = max(1, padded_start - args.reads_length)
    large_padded_end = min(chrom_len, padded_end + args.reads_length)
    mut_by_pos = {}
    for alt in alt_records:
        for pos, variant in alt["mut"].items():
            mut_by_pos[pos] = variant
    nb_random_simulations = 1 if len(mut_by_pos) == 0 else args.nb_random_simulations
    log.debug('Start {} random simulations on target {} cointaining {} variants'.format(nb_random_simulations, target, len(mut_by_pos)))
    best_simu_reads = []
    best_simu_AF_by_var = {}
    lowest_simu_error = np.inf
    for heuristic_idx in range(nb_random_simulations):
        log.debug('Simulation {}/{}'.format(heuristic_idx + 1, nb_random_simulations))
        cp_fragments_length = [elt for elt in fragments_len]
        res_expected = []
        res_simulated = []
        target_reads = []
        res_simu_by_var = {}
        start_dp = random.randrange(args.min_depth, args.max_depth, 1)
        # Simulate
        fragment_idx = 1
        reads_overlapping = {curr_pos: [] for curr_pos in range(large_padded_start, large_padded_end + 1)}  # by pos permet d'ajouter le mate quand on fait le premier et que autre apres un insert
        for curr_pos in range(padded_start, padded_end + 1):
            # Proba to add reads
            adding_proba = None
            if curr_pos < target.start:  # On upstream padding
                distance_to_target = target.start - curr_pos
                depth_ratio = distance_to_target / args.targets_padding
                adding_proba = random.randrange(0, int(start_dp * depth_ratio) + 1, 1)
            elif curr_pos > target.end:  # On downstream padding
                distance_to_target = curr_pos - target.end
                depth_ratio = (args.targets_padding - distance_to_target) / args.targets_padding
                adding_proba = random.randrange(0, int(args.min_depth * depth_ratio) + 1, 1)
            else:  # On target
                adding_proba = random.randrange(args.min_depth, args.max_depth + 1, 1)
            # Add reads
            reads_overlapping_pos = reads_overlapping[curr_pos]
            if adding_proba >= len(reads_overlapping_pos):
                nb_added = random.randrange(args.min_add_depth, args.max_add_depth, 1)
                added_alt_record = random.sample(alt_records, nb_added)
                for curr_alt_record in added_alt_record:
                    # Create fragment
                    curr_fragment_len = cp_fragments_length.pop(0)
                    fragment_record = getFragmentRegion(chrom_seq, target, curr_alt_record["seq"], curr_pos, curr_fragment_len)
                    fragment_record.name = fragment_idx
                    fragment_record.strand = random.choice(["+", "-"])
                    # Create reads pair
                    setReadsPair(fragment_record, args.reads_length, args.R1_end_adapter, args.R2_end_adapter)
                    target_reads.append(fragment_record.annot["reads"])
                    # Add in reads_overlapping (del are occupied by reads)
                    fragment_record.annot["mut"] = curr_alt_record["mut"]
                    updateOverlapping(reads_overlapping, fragment_record)
                    # Next fragment
                    fragment_idx += 1
            # Score simulation
            if curr_pos >= target.start and curr_pos <= target.end:
                if curr_pos in mut_by_pos:
                    variant = mut_by_pos[curr_pos]
                    AD = 0
                    for read, var_in_read in reads_overlapping_pos:
                        if var_in_read is not None:
                            AD += 1
                    AF = round(AD / len(reads_overlapping_pos), 6)
                    res_expected.append(variant["AF"])
                    res_simulated.append(AF)
                    res_simu_by_var[variant["id"]] = AF
            # Next pos
            reads_overlapping[curr_pos] = None
        # Compare simulation
        curr_simu_error = simulationPenaltyScore(res_simulated, res_expected)
        if curr_simu_error <= lowest_simu_error:
            lowest_simu_error = curr_simu_error
            best_simu_reads = target_reads
            best_simu_AF_by_var = res_simu_by_var
    return best_simu_reads, best_simu_AF_by_var


def updateVariantsAF(variants_by_pos, simu_AF_by_id):
    """
    """
    # Get variants objects by ID
    variant_by_id = {}
    for pos, variant in variants_by_pos.items():
        variant_by_id[variant["id"]] = variant
    # Add simulated AF in variants objects
    for variant_id, simu_AF in simu_AF_by_id.items():
        variant_by_id[variant_id]["simu_AF"] = simu_AF
    return variants_by_pos


def getFragmentsLengths(covered_len, args):
    """
    """
    sequenced_nt = int(covered_len * np.mean([args.min_depth, args.max_depth]))  # without padding
    approx_nb_reads = sequenced_nt / (2 * args.reads_length)
    approx_nb_reads = int(approx_nb_reads * 1.2) ################################################
    log.debug('Simulation targeted_nt:{}, sequenced_nt_on_target:{}, approx_nb_reads: {}'.format(covered_len, sequenced_nt, approx_nb_reads))
    fragments_len = np.random.normal(args.fragments_length, args.fragments_length_sd, approx_nb_reads)
    fragments_len = [int(round(elt, 0)) for elt in fragments_len]  # To int
    return fragments_len


class LoggerAction(argparse.Action):
    """Manages logger level parameters (The value "INFO" becomes logging.info and so on)."""
    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Creates simulated reads sequences corresponding to the reference and a mutation profile. Sequencing error, PCR duplications are not applied in this step.')
    parser.add_argument('--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-s', '--random-seed', type=int, default=int(time.time()), help="The seed used for the random generator. If you want reproduce results of one execution: use the same parameters AND the same random-seed. [Default: auto]")
    parser.add_argument('-n', '--nb-random-simulations', type=int, default=30, help='Number of simulations used to find the best reads set to optimize variants frequencies simulated compared to expected. [Default: %(default)s]')
    parser.add_argument('-d', '--min-distance', type=int, default=3, help="The minimum distance between two variants. [Default: %(default)s]")
    parser.add_argument('-g', '--targets-padding', type=int, default=150, help="Padding introduced by panel around the primary targets. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_reads = parser.add_argument_group('Reads')  # Reads
    group_reads.add_argument('-1', '--R1-end-adapter', default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGCGGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA", help='The sequence of the Illumina p7 adapter + flowcell anchor (stretch of 10 A). It is found at the end of the R1 when the read length is superior to the fragment length. [Default: %(default)s]')
    group_reads.add_argument('-2', '--R2-end-adapter', default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTGACAAGCGCTTGTCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAA", help='The reverse complemented sequence of the Illumina p5 adapter + flowcell anchor (stretch of 10 A). It is found at the end of the R2 when the read length is superior to the fragment length. [Default: %(default)s]')
    group_reads.add_argument('-l', '--reads-length', type=int, default=150, help='The reads length. [Default: %(default)s]')
    group_reads.add_argument('-f', '--fragments-length', type=int, default=180, help='The fragments lengths. [Default: %(default)s]')
    group_reads.add_argument('-e', '--fragments-length-sd', type=int, default=30, help='The standard deviation for fragments lengths. [Default: %(default)s]')
    group_depths = parser.add_argument_group('Depths')  # Depths
    group_depths.add_argument('--min-depth', type=int, default=200, help='The minimum depth on targets. [Default: %(default)s]')
    group_depths.add_argument('--max-depth', type=int, default=15000, help='The maximum depth on targets. [Default: %(default)s]')
    group_depths.add_argument('--min-add-depth', type=int, default=5, help='When the depth increase between two consecutive positions it increase to this minimum number of reads. [Default: %(default)s]')
    group_depths.add_argument('--max-add-depth', type=int, default=15, help="When the depth increase between two consecutive positions it can't increase more than this number of reads. [Default: %(default)s]")
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-t', '--input-targets', required=True, help='Path to the targets (format: BED).')
    group_input.add_argument('-r', '--input-reference', required=True, help='Path to the reference sequences file (format: fasta).')
    group_input.add_argument('-p', '--input-profile', required=True, help='Path to the variants profile file (format: TSV).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('--output-variants', required=True, help='Path to the file describing variants introduced in reads and their AF (format: VCF).')
    group_output.add_argument('--output-R1', required=True, help='Path to the R1 output file.')
    group_output.add_argument('--output-R2', required=True, help='Path to the R2 output file.')
    args = parser.parse_args()

    # Logger initialisation
    logging.basicConfig(level=args.logging_level, format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info(" ".join(sys.argv))

    # Set random seed
    log.info("Random seed used: {}".format(args.random_seed))
    random.seed(args.random_seed)
    np.random.seed(args.random_seed)

    # Load variant profile
    log.info("Get variants profile")
    min_allele_freq = 0.01
    models = getVariantsProfile(args.input_profile, min_allele_freq)

    # Get reference sequences
    log.info("Get covered sequences")
    targets_by_chr = getAreasByChr(args.input_targets)
    for chrom_id, targets in targets_by_chr.items():
        mergeOverlapped(targets)

    # Find variants positions
    log.info("Get variants positions")
    mutable_pos = getCoveredList(targets_by_chr)
    covered_len = len(mutable_pos)
    variant_by_pos = dict()
    var_idx = 0
    for model in models:  # For each type of variant
        nb_variants = int(model["occurence"] * covered_len)  # Number of variants on targets
        for idx in range(nb_variants):
            var_idx += 1
            mutable_pos_len = len(mutable_pos)
            variant_len = random.randint(model["length"]["min"], model["length"]["max"])
            # Find mutated positions
            ########################################### fct
            pos_start_idx = None
            pos_end_idx = None
            variant_start_chr = None
            variant_end_chr = None
            valid_position = False
            while not valid_position:
                pos_start_idx = random.randint(0, mutable_pos_len - 1)
                pos_end_idx = pos_start_idx + variant_len - 1 if model["type"] != "insertion" else pos_start_idx
                variant_start_chr, variant_start_pos = mutable_pos[pos_start_idx].split(":")
                variant_start_pos = int(variant_start_pos)
                if pos_end_idx < mutable_pos_len:  # not out of index
                    variant_end_chr, variant_end_pos = mutable_pos[pos_end_idx].split(":")
                    variant_end_pos = int(variant_end_pos)
                    if variant_len == 1 or model["type"] == "insertion":
                        valid_position = True
                    elif variant_end_chr == variant_start_chr:  # start and end are on the same chromosome
                        if variant_end_pos == (variant_start_pos + variant_len - 1):  # The positions are continuous on chromosome
                            valid_position = True
            variant_chr = variant_start_chr
            ###########################################
            # Remove mutated pos and margin from mutable positions
            masked_indexes = [pos_idx for pos_idx in range(pos_start_idx, pos_end_idx + 1)]  # List of indexes to remove from the mutable positions
            for margin_offset in range(1, (args.min_distance + 1)):
                if pos_end_idx + margin_offset < mutable_pos_len:
                    end_chr = mutable_pos[pos_end_idx + margin_offset].split(":")[0]
                    if variant_chr == end_chr:
                        masked_indexes.append(pos_end_idx + margin_offset)
                if pos_start_idx - margin_offset >= 0:
                    start_chr = mutable_pos[pos_start_idx - margin_offset].split(":")[0]
                    if variant_chr == start_chr:
                        masked_indexes.insert(0, pos_start_idx - margin_offset)
            masked_indexes = sorted(masked_indexes, reverse=True)
            clean_variant_len = variant_len if model["type"] != "insertion" else 1
            if(len(masked_indexes) < (clean_variant_len + 2 * args.min_distance)):
                log.debug(
                    "The mask has been limited to {} on position {}:{}-{}.".format(
                        ", ".join([mutable_pos[idx] for idx in masked_indexes]),
                        variant_chr,
                        variant_start_pos,
                        variant_end_pos
                    )
                )
            for idx in masked_indexes:
                del mutable_pos[idx]
            # Store variant information
            variant_freq = random.randint(
                int((1 / min_allele_freq) * model["AF"]["min"]),
                int((1 / min_allele_freq) * model["AF"]["max"])
            ) / (1 / min_allele_freq)
            if variant_chr not in variant_by_pos:
                variant_by_pos[variant_chr] = dict()
            variant_by_pos[variant_chr][variant_start_pos] = {
                "id": var_idx,
                "chrom": variant_chr,
                "start": variant_start_pos,
                "alt": None,
                "ref": None,
                "type": model["type"],
                "length": variant_len,
                "AF": variant_freq
            }

    # Clean previous files
    for filepath in [args.output_R1, args.output_R2]:
        if os.path.exists(filepath):
            os.remove(filepath)

    # Create reads
    log.info("Create reads")
    fragments_len = getFragmentsLengths(covered_len, args)
    targets_by_chr = getAreasByChr(args.input_targets)
    mergeOverlappedWithPadding(targets_by_chr, args.targets_padding)
    for chrom_id, targets in sorted(targets_by_chr.items()):
        log.info("Get sequence from region {}.".format(chrom_id))
        chrom_seq = getChrSeq(args.input_reference, chrom_id)
        chrom_len = len(chrom_seq)
        for curr_target in targets:
            log.info("Create reads from targeted region {}".format(curr_target))
            target_ref_seq = chrom_seq[curr_target.start - 1:curr_target.end]  # Position is 1-based indexes are 0-based
            # Apply variants
            alt_records = [{"id": idx + 1, "mut": {}, "seq": list(target_ref_seq)} for idx in range(int(1 / min_allele_freq))]
            variants_on_target = []
            if chrom_id in variant_by_pos:
                variants_on_target = selectMatchedVariants(curr_target, variant_by_pos[chrom_id])
            applyVariants(curr_target, list(target_ref_seq), alt_records, {mut["start"]: mut for mut in variants_on_target})
            # Create reads
            simulated_reads_pairs, simu_AF_by_var = getTargetReads(chrom_seq, chrom_len, curr_target, alt_records, fragments_len, args)
            writeTargetReads(args.output_R1, args.output_R2, simulated_reads_pairs)
            if len(simu_AF_by_var) > 0:
                updateVariantsAF(variant_by_pos[chrom_id], simu_AF_by_var)

    # Write variants
    log.info("Write variants trace")
    writeVariants(args.output_variants, variant_by_pos)

    log.info("Simulation success")
