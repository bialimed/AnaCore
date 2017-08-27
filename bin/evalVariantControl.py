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
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from vcf import VCFIO, getAlleleRecord



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
                       "chr1:10=T":{
                         "chrom":"chr1",
                         "pos":10,
                         "ref":"A",
                         "alt":"T",
                         "freq":[0.2, 0.5] },
                       "chr1:10=G":{
                         "chrom":"chr1",
                         "pos":10,
                         "ref":"A",
                         "alt":"G",
                         "freq":[0.01, 0] },
                       "chr3:20=T":{
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
    with VCFIO( vcf_path ) as FH_vcf:
        if spl_name is None:
            spl_name = FH_vcf.samples[0]
        for record in FH_vcf:
            allele_freq = record.getAF( spl_name )
            # For each alternative allele
            for idx_alt, alt in enumerate(record.alt):
                allele_record = getAlleleRecord( FH_vcf, record, idx_alt )
                allele_record.standardizeSingleAllele()
                variant_id = allele_record.chrom + ":" + str(allele_record.pos) + "=" + allele_record.alt[0]
                if variant_id not in variants:
                    variants[variant_id] = {
                        "chrom": allele_record.chrom,
                        "pos": allele_record.pos,
                        "ref": allele_record.ref,
                        "alt": allele_record.alt[0],
                        "freq": list()
                    }
                # Complete variants missing in previous VCF
                while len(variants[variant_id]["freq"]) <= vcf_idx:
                    variants[variant_id]["freq"].append(0)
                # Add allele frequency
                variants[variant_id]["freq"][vcf_idx] = allele_freq[idx_alt]
    # Complete variants missing in current VCF
    for variant_id in variants:
        while len(variants[variant_id]["freq"]) <= vcf_idx:
            variants[variant_id]["freq"].append(0)

def filterVarNotIn( variants, vcf_idx ):
    """
    @summary: Removes all variants not present in the specified VCF.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param vcf_idx: [int] Index used to store the frequency of each vrariants of the VCF in frequencies list (start from 0).
    """
    removed = list()
    for variant_id in variants:
        if variants[variant_id]["freq"][vcf_idx] == 0:
            removed.append( variant_id )
    for variant_id in removed:
        del( variants[variant_id] )

def writeTSVResults(variants, out_path, error_threshold=0.2, separator="\t"):
    """
    @summary: Writes expected and detected fequency for each expected variant in TSV file.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param out_path: [str] Path to the output file.
    @param error_threshold: [float] The maximum percentage difference between expected and detected.
    @param separator: [str] The  column separator in output file.
    """
    idx_expec = 0
    idx_detec = 1

    # Process summary
    nb_checked = len(variants)
    nb_out_threshold = 0
    error_ratio_sum = 0
    error_values_sum = 0
    for variant_id in variants:
        current_variant = variants[variant_id]
        error = current_variant["freq"][idx_expec] - current_variant["freq"][idx_detec]
        error_ratio = abs(1 - (current_variant["freq"][idx_detec]/current_variant["freq"][idx_expec]) )
        # Dissimilarity
        error_values_sum += abs(error)
        error_ratio_sum += error_ratio
        # Error out of threshold
        if error_ratio > error_threshold:
            nb_out_threshold += 1

    # Write
    with open(out_path, "w") as FH_out:
        FH_csv = csv.writer(FH_out, delimiter=separator)
        # Summary section
        FH_csv.writerow( ["[Summary]"] )
        FH_csv.writerow( ["#Nb_checked", "Errors_sum", "Errors_ratio_sum", "Error_out_of_threshold_(" + str(error_threshold) + ")"] )
        FH_csv.writerow([
            nb_checked,
            error_values_sum,
            error_ratio_sum,
            str(nb_out_threshold) + "/" + str(nb_checked) ])
        # Section separator
        FH_csv.writerow([])
        # Details section
        FH_csv.writerow( ["[Details]"] )
        FH_csv.writerow( ["#Chr:Pos", "Ref/Alt", "Expected", "Detected", "Error", "Error_ratio", "Out_of_threshold"] )
        for variant_id in sorted(variants):
            current_variant = variants[variant_id]
            error = current_variant["freq"][idx_expec] - current_variant["freq"][idx_detec]
            error_ratio = abs(1 - (current_variant["freq"][idx_detec]/current_variant["freq"][idx_expec]) )
            FH_csv.writerow([
                current_variant["chrom"] + ":" + str(current_variant["pos"]),
                current_variant["ref"] + "/" + current_variant["alt"],
                current_variant["freq"][idx_expec],
                current_variant["freq"][idx_detec],
                error,
                error_ratio,
                (error_ratio > error_threshold) ])

def writeJSONResults(variants, out_path):
    """
    @summary: Writes expected and detected fequency for each expected variant in JSON file.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param out_path: [str] Path to the output file.
    """
    data = list()
    for variant_id in variants:
        current_variant = variants[variant_id]
        data.append({
            "chrom": current_variant["chrom"],
            "pos": current_variant["pos"],
            "ref": current_variant["ref"],
            "alt": current_variant["alt"],
            "expected": current_variant["freq"][0],
            "detected": current_variant["freq"][1]
        })
    with open(out_path, "w") as FH_out:
        FH_out.write( json.dumps(data, default=lambda o: o.__dict__, sort_keys=True ) )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Compare variant calling result to expected variants. This comparison is only processed on expected variants.' )
    parser.add_argument( '-t', '--error-threshold', type=float, default=0.2, help='In TSV output, variants with a percentage difference compared to expected frequency superior than this value are tagged "out of threshold". Difference percentage calculation: abs(1 - detected_freq/expected_freq). Examples: diff=50prct for expected_freq=0.1 and detected_freq=0.05 ; diff=300prct for expected_freq=0.1 and detected_freq=0.4. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-e', '--expected-file', required=True, help='The path to the file containing the expected variants (format: VCF).' )
    group_input.add_argument( '-d', '--detected-file', required=True, help='The path to the file containing the detected variants (format: VCF).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', default="results.tsv", help='The path to the file containing the results of the comparison (format: TSV or JSON depends on extension). [Default: %(default)s]')
    args = parser.parse_args()

    # Process
    all_variants = dict()
    addVCFVariants( all_variants, args.expected_file, 0 )
    addVCFVariants( all_variants, args.detected_file, 1 )
    filterVarNotIn( all_variants, 0 )
    if args.output_file.endswith("json"):
        writeJSONResults( all_variants, args.output_file )
    else:
        writeTSVResults( all_variants, args.output_file, args.error_threshold )
