#!/usr/bin/env python3
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
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def getBiggerCluster(freq_by_spl, max_freq_diff):
    """
    Returns the bigger group of samples with similar frequencies.

    :param freq_by_spl: The frequency by sample.
    :type freq_by_spl: dict
    :param max_freq_diff: The maximum difference of frequencies in a cluster.
    :type max_freq_diff: float
    :return: The list of samples in bigger cluster. The samples are sorted by frequency.
    :rtype: list
    """
    bigger_cluster_size = 0
    bigger_clster_samples = list()
    asc_freq_spl = sorted(freq_by_spl, key=freq_by_spl.get)
    spl_len = len(asc_freq_spl)
    for start_idx, start_spl in enumerate(asc_freq_spl):
        clstr_size = 0
        clstr_samples = list()  # The start sample itself will be added in comparison below
        start_AF = freq_by_spl[start_spl]
        if len(asc_freq_spl) - start_idx > bigger_cluster_size:
            curr_idx = start_idx
            out_of_range = False
            while curr_idx < spl_len - 1 and not out_of_range:
                curr_spl = asc_freq_spl[curr_idx]
                curr_spl_AF = freq_by_spl[curr_spl]
                if curr_spl_AF - start_AF <= max_freq_diff:
                    clstr_size += 1
                    clstr_samples.append(curr_spl)
                else:
                    out_of_range = True
                curr_idx += 1
        if clstr_size > bigger_cluster_size:
            bigger_cluster_size = clstr_size
            bigger_clster_samples = clstr_samples
    return bigger_clster_samples


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Find constitutive variants in population of samples coming from one VCF. An constitutive variant is a variant coming from polymerase error rate at this position, or workflow artifact, ... . It is common in almost all samples with an allele frequency very similar in each samples. If this variant is present in a superiror frequency in a particular sample that can signify that this variant is not an artefact in this sample only.')
    parser.add_argument('-n', '--noise-offset', type=float, default=0.01, help='The value added to the maximum constitutive frequency of the variant.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_thresholds = parser.add_argument_group('Thresolds')  # Thresholds
    group_thresholds.add_argument('-r', '--min-presence-ratio', type=float, default=0.75, help='Minimum ratio between samples with variant at the same AF and the samples with valid depth (see "--min-DP") for declare variant as constitutive.')
    group_thresholds.add_argument('-c', '--min-presence-count', type=int, default=4, help='Minimum number of samples with variant at the same AF for declare variant as constitutive.')
    group_thresholds.add_argument('-d', '--min-DP', type=int, default=150, help='Minimum number of reads at the variant position for use the sample in constitutive variant evaluation.')
    group_thresholds.add_argument('-a', '--max-AF-var', type=float, default=0.015, help='Maximum variation between samples AF in constitutive variant. Constitutive variant have an AF similar between all samples, for determine if a variant is constitutive only the samples where the variant has an AF similar in comparison to this value are counted as support.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-TSV', required=True, help='The path to the outputted file containing the constitutive variants (format: TSV).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.input_variants) as FH_in:
        with open(args.output_variants, "w") as FH_out:
            # Header
            FH_out.write("## PARAMETERS: {}\n".format(
                " ".join(sys.argv)
            ))
            FH_out.write("## VERSION: {}\n".format(__version__))
            FH_out.write("\t".join([
                "#Chromosome",
                "Position",
                "Reference_allele",
                "Alternative_allele",
                "Noise_rate",
                "Nb_input_spl",
                "Nb_usable_spl",
                "Nb_support_spl",
                "Nb_constit_spl",
                "Constit_spl",
                "Constit_AF"
            ]))
            # Records
            for record in FH_in:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(FH_in, record, idx)
                    nb_spl = len(FH_in.samples)
                    nb_usable_spl = 0  # Nb samples with valid DP at the variant position
                    nb_support_spl = 0  # Nb samples containing the variant
                    nb_support_constit_spl = 0
                    AF_by_spl = dict()
                    for curr_spl in FH_in.samples:
                        curr_DP = curr_allele.getDP(curr_spl)
                        if curr_DP is None:
                            raise Exception('The DP for the variant "{}" in sample "{}" must be set when you merge samples.'.format(curr_allele.getName(), curr_spl))
                        elif curr_DP > args.min_DP:  # Skip samples with limited DP at the variant position
                            nb_usable_spl += 1
                            curr_AF = curr_allele.getAF(curr_spl)[0]
                            if curr_AF is not None and curr_AF > 0:  # The sample contain the variant
                                nb_support_spl += 1
                                AF_by_spl[curr_spl] = curr_AF
                    if nb_support_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_spl >= args.min_presence_count:  # Prefilter by support number for reduce processing time
                        clstr_samples = getBiggerCluster(AF_by_spl, args.max_AF_var)
                        nb_support_constit_spl = len(clstr_samples)
                        if nb_support_constit_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_constit_spl >= args.min_presence_count:  # The variant is constitutive
                            clstr_AF = [AF_by_spl[clstr_samples] for spl in sorted(clstr_samples)]
                            FH_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                curr_allele.chrom,
                                curr_allele.pos,
                                curr_allele.ref,
                                curr_allele.alt[0],
                                max(clstr_AF) + args.noise_offset,
                                nb_spl,
                                nb_usable_spl,
                                nb_support_spl,
                                nb_support_constit_spl,
                                ";".join(sorted(clstr_samples)),
                                ";".join(clstr_AF)
                            ))
