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
from anacore.sv import HashedSVIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getNoise(input_noise):
    """
    Return by variant id ("chrom:pos=ref/alt") the noise rate.

    :param input_noise: The path to the file containing artifactual variants with their maximum frequency (format: TSV). The header line of the file must be "#Chromosome<tab>Possition<tab>Reference_allele<tab>Alternative_allele<tab>Noise_rate".
    :type input_noise: str
    :return: By variant id ("chrom:pos=ref/alt") the noise rate.
    :rtype: dict
    """
    expected_titles = ["Chromosome", "Position", "Reference_allele", "Alternative_allele", "Noise_rate"]
    noise_by_var = dict()
    with HashedSVIO(input_noise) as FH_noise:
        if FH_noise.titles != expected_titles:
            raise Exception(
                'The header line in "{}" does not correpond to "#{}".'.format(
                    input_noise, "\t".join(expected_titles)
                )
            )
        for record in FH_noise:
            variant_id = "{}:{}={}/{}".format(record["Chromosome"], record["Position"], record["Reference_allele"], record["Alternative_allele"])
            noise_by_var[variant_id] = float(record["Noise_rate"])
    return noise_by_var


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Tags or removes the artifactual variants. If a variant is known as noise and it has an AF lower or equal than the noise rate it is tagged/removed. In contrast, if it is known as noise but his AF is upper than the noise rate it is kept (signal is upper than the noise). This script is very useful to tag/remove variants coming from sequencing polymerase error rate for example on homopolymer.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant is noise in all samples present in the VCF the a tag is added in FILTER (cf. "--tag-name"). In mode "remove" if the variant is noise all samples present in the VCF this variant is removed. [Default: %(default)s]')
    group_filter.add_argument('-tn', '--tag-name', default="popConst", help='The name of the tag added on variant if it correspond to noise. [Default: %(default)s]')
    group_filter.add_argument('-td', '--tag-description', default="The variant correspond to a constitutive detection (sequencing polymerase error rate at this position, workflow artifact, ...).", help='The description of the applied filter. This description is stored in header of the VCF. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_input.add_argument('-n', '--input-noises', required=True, help='The path to the file containing artifactual variants with their maximum frequency (format: TSV). The header line of the file must be "#Chromosome<tab>Possition<tab>Reference_allele<tab>Alternative_allele<tab>Noise_rate".')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Process
    noise_by_variant = getNoise(args.input_noises)
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.filter[args.tag_name] = args.tag_description
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(FH_in, record, idx)
                    # Compare signal to noise
                    if curr_allele.getName() in noise_by_variant:
                        nb_spl_over_noise = 0
                        for curr_spl in curr_allele.samples:
                            if curr_allele.getAltAF(curr_spl)[0] > noise_by_variant[curr_allele.getName()]:
                                nb_spl_over_noise += 1
                        if nb_spl_over_noise == 0:
                            if curr_allele.filter is None or len(curr_allele.filter) == 0:
                                curr_allele.filter = [args.tag_name]
                            else:
                                if "PASS" in curr_allele.filter:
                                    curr_allele.filter.remove("PASS")
                                curr_allele.filter.append(args.tag_name)
                    # Update empty filter
                    if curr_allele.filter is None or len(curr_allele.filter) == 0:
                        curr_allele.filter = ["PASS"]
                    # Write record
                    if args.mode == "tag":
                        FH_out.write(curr_allele)
                    elif args.tag_name not in curr_allele.filter:
                        FH_out.write(curr_allele)
