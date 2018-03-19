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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.    If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import warnings
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
def getWorkflowVariants(wf_folder):
    wf_variants = dict()
    for file_name in os.listdir(wf_folder):
        filepath = os.path.join(wf_folder, file_name)
        if os.path.isdir(filepath):
            wf_variants[file_name] = dict()
            variants_folder = os.path.join(filepath, "variants")
            if not os.path.exists(variants_folder):
                warnings.warn('The folder "' + variants_folder + '" does not exist.')
            else:
                for variants_file_name in os.listdir(variants_folder):
                    variants_filepath = os.path.join(variants_folder, variants_file_name)
                    if variants_filepath.endswith("_filtered.vcf"):
                        spl_name = variants_file_name.split("_filtered.vcf")[0]
                        spl_variants = dict()
                        with VCFIO(variants_filepath) as FH_variants:
                            for record in FH_variants:
                                for idx in range(len(record.alt)):
                                    curr_allele = getAlleleRecord(FH_variants, record, idx)
                                    variant_id = curr_allele.chrom + ":" + str(curr_allele.pos) + "=" + curr_allele.ref + "/" + "/".join(curr_allele.alt)
                                    spl_variants[variant_id] = curr_allele
                        wf_variants[file_name][spl_name] = spl_variants
    return wf_variants



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Compares DSVF and ADIVaR filtered variants.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--adivar-folder', required=True, help='The path to the folder containing one folder by run outputs for ADIVaR.')
    group_input.add_argument('-d', '--dsvf-folder', required=True, help='The path to the folder containing one folder by run outputs for DSVF.')
    args = parser.parse_args()

    # Process
    dsvf_variants = getWorkflowVariants(args.dsvf_folder)
    adivar_variants = getWorkflowVariants(args.adivar_folder)

    print(
      "#Run",
      "Sample",
      "Position",
      "Ref/Alt",
      "ADIVaR_AF",
      "DSVF_AF",
      "Diff_AF",
      "Ratio_diff_AF",
      sep="\t"
    )
    for curr_run in adivar_variants:
        if curr_run in dsvf_variants:
            for curr_spl in adivar_variants[curr_run]:
                if curr_spl in dsvf_variants[curr_run]:
                    dsvf_spl_var = dsvf_variants[curr_run][curr_spl]
                    adivar_spl_var = adivar_variants[curr_run][curr_spl]
                    processed = dict()
                    for variant_id in adivar_spl_var:
                        processed[variant_id] = 1
                        adivar_AF = round(adivar_spl_var[variant_id].getPopAF()[0], 2)
                        dsvf_AF = round((dsvf_spl_var[variant_id].getPopAF()[0] if variant_id in dsvf_spl_var else 0), 2)
                        print(
                            curr_run,
                            curr_spl,
                            variant_id.split("=")[0],
                            variant_id.split("=")[1],
                            adivar_AF,
                            dsvf_AF,
                            max(adivar_AF, dsvf_AF) - min(adivar_AF, dsvf_AF),
                            1 - (min(adivar_AF, dsvf_AF) / max(adivar_AF, dsvf_AF)),
                            sep="\t"
                        )
                    for variant_id in dsvf_spl_var:
                        if variant_id not in processed:
                            adivar_AF = 0.00
                            dsvf_AF = round(dsvf_spl_var[variant_id].getPopAF()[0], 2)
                            print(
                                curr_run,
                                curr_spl,
                                variant_id.split("=")[0],
                                variant_id.split("=")[1],
                                adivar_AF,
                                dsvf_AF,
                                max(adivar_AF, dsvf_AF) - min(adivar_AF, dsvf_AF),
                                1 - (min(adivar_AF, dsvf_AF) / max(adivar_AF, dsvf_AF)),
                                sep="\t"
                            )
