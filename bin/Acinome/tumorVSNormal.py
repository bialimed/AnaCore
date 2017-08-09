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
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import csv
import json
import warnings
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import *



########################################################################
#
# FUNCTIONS
#
########################################################################
def addVCFVariants( variants, vcf_path, vcf_idx, spl_name=None ):
    """
    @summary: Add variant from VCF in dict.
    @param variants: [dict] By uniq ID the variants. The content of this variable is set by the call of this function.
                     Content example:
                     {
                       "chr1:10=A/T":{
                         "chrom":"chr1",
                         "pos":10,
                         "ref":"A",
                         "alt":"T",
                         "freq":[0.2, 0.5] },
                       "chr1:10=A/G":{
                         "chrom":"chr1",
                         "pos":10,
                         "ref":"A",
                         "alt":"G",
                         "freq":[0.01, 0] },
                       "chr3:20=T/G":{
                         "chrom":"chr3",
                         "pos":20,
                         "ref":"G",
                         "alt":"T",
                         "freq":[0, 0.4] }
                     }
                     The list of frequencies is appended by each call of the function with a vcf_idx different.
    @param vcf_path: [str] Path to the VCF file to add.
    @param vcf_idx: [int] Index used to store the frequency of each vrariants of the VCF in frequencies list (start from 0).
    @param spl_name: [str] The frequency of the variants came from this sample. This parameters is optional when the VCF file contain 0 to 1 sample.
    """
    FH_vcf = VCFIO( vcf_path )
    try:
        if spl_name is None:
            spl_name = FH_vcf.samples[0]
        for record in FH_vcf:
            # Get allele frequency
            allele_freq = list()
            #~ if len(record.samples) <= 1 and "AF" in record.info:
                #~ allele_freq = record.info["AF"]
            #~ elif "AF" in record.samples[spl_name]:
            if "AF" in record.samples[spl_name]:
                allele_freq = record.samples[spl_name]["AF"]
            elif "AD" in record.samples[spl_name]:
                for allele_depth in record.samples[spl_name]["AD"][1:]: # Skip the reference allele depth
                    allele_freq.append( allele_depth/float(record.samples[spl_name]["DP"]) )
            ####################################
            elif len(record.samples) <= 1 and "AF" in record.info:
                warnings.warn( 'The allele frequency for variant "' + record.chrom + ":" + str(record.pos) + '" cannot be retrieved from the sample "' + spl_name + '". The AF will be extract from INFO.' )
                allele_freq = record.info["AF"]
            ####################################
            else:
                raise Exception( 'The allele frequency cannot be retrieved in variant "' + record.chrom + ":" + str(record.pos) + '".' )
            if not isinstance(allele_freq, (list, tuple)):
                allele_freq = [allele_freq]
            # Add alleles
            for idx_alt, alt in enumerate(record.alt):
                variant_id = record.chrom + ":" + str(record.pos) + "=" + record.ref + "/" + alt
                if variant_id not in variants:
                    variants[variant_id] = {
                        "chrom": record.chrom,
                        "pos": record.pos,
                        "ref": record.ref,
                        "alt": alt,
                        "id": record.id,
                        "filters": list(),
                        "freq": list(),
                        "AD": list()
                    }
                # Complete variants missing in previous VCF
                while len(variants[variant_id]["freq"]) <= vcf_idx:
                    variants[variant_id]["freq"].append(0)
                    variants[variant_id]["AD"].append(0)
                    variants[variant_id]["filters"].append(None)
                # Add allele frequency
                variants[variant_id]["freq"][vcf_idx] = allele_freq[idx_alt]
                if "DP" in record.samples[spl_name]:
                    variants[variant_id]["AD"][vcf_idx] = record.samples[spl_name]["DP"] * allele_freq[idx_alt]
                else:
                    variants[variant_id]["AD"][vcf_idx] = record.info["DP"] * allele_freq[idx_alt]
                variants[variant_id]["filters"][vcf_idx] = record.filter
    finally:
        FH_vcf.close()
    # Complete variants missing in current VCF
    for variant_id in variants:
        while len(variants[variant_id]["freq"]) <= vcf_idx:
            variants[variant_id]["freq"].append(0)
            variants[variant_id]["AD"].append(0)
            variants[variant_id]["filters"].append(None)

def filterVarNotInTumor( variants, min_ratio, min_support ):
    """
    @summary: Removes all variants not present in the specified VCF.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param min_ratio: [float] Minimum allele ratio between tumor and normal VCF to keep an allele. To keep the variant the AF_tumor must be >= min_ratio * AF_normal.
    @param min_support: [int] Minimum number of reads supporting the allele in tumor VCF.
    """
    removed = list()
    for variant_id in variants:
        if not(variants[variant_id]["freq"][0] >= (min_ratio * variants[variant_id]["freq"][1]) and variants[variant_id]["AD"][0] >= min_support):
            removed.append( variant_id )
    for variant_id in removed:
        del( variants[variant_id] )

def writeTSVResults(variants, out_path, separator="\t"):
    """
    @summary: Writes expected and detected fequency for each expected variant in TSV file.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param out_path: [str] Path to the output file.
    @param error_threshold: [float] The minimum similarity ratio between expected and detected.
    @param separator: [str] The  column separator in output file.
    """
    # Write
    with open(out_path, "w") as FH_out:
        FH_csv = csv.writer(FH_out, delimiter=separator)
        FH_csv.writerow( ["#ID", "Chr:Pos", "Ref/Alt", "Tumor_freq", "Tumor_AD", "Tumor_filters", "Normal_freq", "Normal_AD", "Normal_filters", "dbSNP"] )
        for variant_id in sorted(variants):
            current_variant = variants[variant_id]
            FH_csv.writerow([
                variant_id,
                current_variant["chrom"] + ":" + str(current_variant["pos"]),
                current_variant["ref"] + "/" + current_variant["alt"],
                '{:.4f}'.format(current_variant["freq"][0]),
                current_variant["AD"][0],
                current_variant["filters"][0],
                '{:.4f}'.format(current_variant["freq"][1]),
                current_variant["AD"][1],
                current_variant["filters"][1],
                current_variant["id"] ])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Compare variant calling result to expected variants. This comparison is only processed on expected variants.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_filter = parser.add_argument_group( 'Filters' ) # Filters
    group_filter.add_argument( '-r', '--min-ratio', default=1.5, type=float, help='Minimum allele ratio between tumor and normal VCF to keep an allele. To keep the variant the AF_tumor must be >= min_ratio * AF_normal. [Default: %(default)s]' )
    group_filter.add_argument( '-s', '--min-support', default=20, type=int, help='Minimum number of reads supporting the allele in tumor VCF. [Default: %(default)s]' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-t', '--tumor-file', required=True, help='The path to the file containing the expected variants (format: VCF).' )
    group_input.add_argument( '-n', '--normal-file', required=True, help='The path to the file containing the detected variants (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', default="results.tsv", help='The path to the file containing the results of the comparison (format: TSV or JSON depends on extension). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    all_variants = dict()
    addVCFVariants( all_variants, args.tumor_file, 0 )
    addVCFVariants( all_variants, args.normal_file, 1 )
    filterVarNotInTumor( all_variants, args.min_ratio, args.min_support )
    writeTSVResults( all_variants, args.output_file )
