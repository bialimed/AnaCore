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
__version__ = '1.0.0'
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
        start_AF = AF_by_spl[start_spl]
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
    parser = argparse.ArgumentParser(description='Tags or removes the constitutive variants. An constitutive variant is a variant coming from sequencing polymerase error rate at this position, or workflow artifact, ... . A constitutive variant is common in almost all samples with an allele frequency very similar in each samples. If this variant is present in a superiror frequency in a particular sample that can signify that this variant is not an artefact in this sample only.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant is constitutive and any samples really contain this variant the tag "popConst" is added in FILTER. If at least one sample really contains this variant The AD and DP of the variant in samples where the variant is artifact are reduce to 0 and the tag "popConst" is added in sample format tag "filter". In mode "remove" if the variant is constitutive and any samples really contain this variant it is removed from the output. If at least one sample really contains this variant The AD and DP of the variant in samples where the variant is artifact are reduce to 0 and the tag "popConst" is added in sample format tag "filter". [Default: %(default)s]')
    group_thresholds = parser.add_argument_group('Thresolds')  # Thresholds
    group_thresholds.add_argument('-r', '--min-presence-ratio', type=float, default=0.8, help='Minimum ratio between samples with variant at the same AF and the samples with valid depth (see "--min-DP") for declare variant as constitutive.')
    group_thresholds.add_argument('-c', '--min-presence-count', type=int, default=4, help='Minimum number of samples with variant at the same AF for declare variant as constitutive.')
    group_thresholds.add_argument('-d', '--min-DP', type=int, default=150, help='Minimum number of reads at the variant position for use the sample in constitutive variant evaluation.')
    group_thresholds.add_argument('-a', '--max-AF-var', type=float, default=0.015, help='Maximum variation between samples AF in constitutive variant. Constitutive variant have an AF similar between all samples, for determine if a variant is constitutive only the samples where the variant has an AF similar in comparison to this value are counted as support. Samples with a constitutive variant with an AF superior than this reference AF are not tagged.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.info["CONSCT"] = {
                "type": int,
                "type_tag": "Integer",
                "number": 1,
                "number_tag": "1",
                "description": "The number of analysed samples where the variant correspond to a constitutive detection: sample AF is <= CONSAF."
            }
            FH_out.info["CONSAF"] = {
                "type": float,
                "type_tag": "Float",
                "number": 2,
                "number_tag": "2",
                "description": "The min and max AF for the constitutive detection. The variant is seen at this AF in majority of samples (see the filter popConst)."
            }
            FH_out.filter["popConst"] = {
                "description": "The variant correspond to a constitutive detection (sequencing polymerase error rate at this position, workflow artifact, ...). Criteria: the variant is present in at least {0.2%}% of samples with an AF variation <= {0.2%}%.".format(args.min_presence_ratio, args.max_AF_var)
            }
            FH_out.format["filter"] = {
                "type": int,
                "type_tag": "Integer",
                "number": 1,
                "number_tag": "1",
                "description": "Specifics filters of the sample."
            }
            FH_out._writeHeader()
            # Records
            print("#Variant\tRemove_type\tMin_constit_AF\tMax_constit_AF\tNb_spl_valid_depth\tNb_spl_valid_depth\tNb_support_spl\tNb_support_spl_with_constAF\tSupport_spl_with_constAF")
            for record in FH_in:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(FH_in, record, idx)
                    curr_allele.info["CONSCT"] = 0
                    nb_valid_spl = 0
                    nb_support_spl = 0
                    AF_by_spl = dict()
                    for curr_spl in FH_in.samples:
                        curr_DP = curr_allele.getDP(curr_spl)
                        if curr_DP is None:
                            raise Exception('The DP for the variant "{}" in sample "{}" must be set when you merge samples.'.format(curr_allele.getName(), curr_spl))
                        elif curr_DP > args.min_DP:
                            nb_valid_spl += 1
                            curr_AF = curr_allele.getAF(curr_spl)[0]
                            if curr_AF is not None and curr_AF > 0:
                                nb_support_spl += 1
                                AF_by_spl[curr_spl] = curr_AF
                    if nb_support_spl / nb_valid_spl >= args.min_presence_ratio and nb_support_spl >= args.min_presence_count:
                        clstr_samples = getBiggerCluster(AF_by_spl, args.max_AF_var)
                        curr_allele.info["CONSCT"] = len(clstr_samples)
                        if curr_allele.info["CONS"] / nb_valid_spl >= args.min_presence_ratio and curr_allele.info["CONSCT"] >= args.min_presence_count:
                            curr_allele.info["CONSAF"] = [
                                AF_by_spl[clstr_samples[0]],
                                AF_by_spl[clstr_samples[-1]]
                            ]
                            remove_type = None
                            if curr_allele.info["CONSCT"] > nb_support_spl:  # The variant is constitutive but certains samples have this variant at an other AF
                                remove_type = "samples"
                                # Remove AD and AF from samples with AF < constit AF
                                if "AF" in curr_allele.samples[curr_spl]:
                                    curr_allele.samples[curr_spl]["AF"] = 0
                                if "AD" in curr_allele.samples[curr_spl]:
                                    curr_allele.samples[curr_spl]["AD"] = 0
                                if "popConst" not in curr_allele.samples[curr_spl]["filter"]:
                                    if len(curr_allele.samples[curr_spl]["filter"]) == 0:
                                        curr_allele.samples[curr_spl]["filter"] = ["popConst"]
                                    else:
                                        curr_allele.samples[curr_spl]["filter"].append("popConst")
                            else:  # The variant is only view with an equal proportion in majority of samples and is removed
                                remove_type = "variant"
                                if curr_allele.filter is None:
                                    curr_allele.filter = []
                                if "popConst" not in curr_allele.filter:
                                    if len(curr_allele.filter) == 0 or curr_allele.filter[0] == "PASS":
                                        curr_allele.filter = ["popConst"]
                                    else:
                                        curr_allele.filter.append("popConst")
                            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                record.getName(), remove_type, curr_allele.info["CONSAF"][0],
                                curr_allele.info["CONSAF"][1], nb_valid_spl, nb_support_spl,
                                curr_allele.info["CONSCT"], ";".join(sorted(clstr_samples))
                            ))
                    if curr_allele.filter is None:
                        curr_allele.filter = ["PASS"]
                    if args.mode == "tag":
                        FH_out.write(curr_allele)
                    elif "popConst" not in curr_allele.filter:
                        FH_out.write(curr_allele)
