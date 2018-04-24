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
__status__ = 'dev'


import os
import sys
import csv
import json
import argparse
import statistics

CURRENT_DIR = os.path.dirname(__file__)
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
        if spl_name is None and len(FH_vcf.samples) != 0:
            spl_name = FH_vcf.samples[0]
        for record in FH_vcf:
            if record.chrom.startswith("chr"):
                record.chrom = record.chrom[3:]
            # Get allele frequency
            allele_freq = list()
            if len(record.samples) <= 1 and "AF" in record.info:
                allele_freq = record.info["AF"]
            elif "AF" in record.samples[spl_name]:
                allele_freq = record.samples[spl_name]["AF"]
            elif "AD" in record.samples[spl_name]:
                for allele_depth in record.samples[spl_name]["AD"][1:]: # Skip the reference allele depth
                    allele_freq.append( allele_depth/float(record.info["DP"]) )
            else:
                raise Exception( 'The allele frequency cannot be retrieved in variant "' + record.chrom + ":" + str(record.pos) + '".' )
            if not isinstance(allele_freq, (list, tuple)):
                allele_freq = [allele_freq]
            # Add alleles
            for idx_alt, alt in enumerate(record.alt):
                alt_record = getAlleleRecord( FH_vcf, record, idx_alt )
                alt_record.ref = record.ref.upper()
                alt_record.alt[0] = alt_record.alt[0].upper()
                alt_record.standardizeSingleAllele()
                variant_id = alt_record.chrom + ":" + str(alt_record.pos) + "=" + alt_record.ref + "/" + alt_record.alt[0]
                if variant_id not in variants:
                    variants[variant_id] = {
                        "chrom": record.chrom,
                        "pos": alt_record.pos,
                        "ref": alt_record.ref,
                        "type": alt_record.type(),
                        "alt": alt_record.alt[0],
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

def writeTSVResults(variants, out_path, error_threshold=0.8, separator="\t"):
    """
    @summary: Writes expected and detected fequency for each expected variant in TSV file.
    @param variants: [dict] By uniq ID the variants (see addVCFVariants).
    @param out_path: [str] Path to the output file.
    @param error_threshold: [float] The minimum similarity ratio between expected and detected.
    @param separator: [str] The  column separator in output file.
    """
    idx_expec = 0
    idx_detec = 1

    # Process summary
    nb_TP = 0
    nb_FP = 0
    nb_FN = 0
    freq_FN = list()
    freq_FP = list()
    nb_checked = len(variants)
    nb_out_threshold = 0
    expect_error_count = list()
    #~ expect_error_count_sum = 0
    expect_error_ratio = list()
    #~ expect_error_ratio_sum = 0
    for variant_id in variants:
        current_variant = variants[variant_id]
        current_variant["error"] = current_variant["freq"][idx_expec] - current_variant["freq"][idx_detec]
        current_variant["error_ratio"] = 1
        if current_variant["freq"][idx_detec] != 0 and current_variant["freq"][idx_expec] != 0:
            current_variant["error_ratio"] = abs(1 - (current_variant["freq"][idx_detec]/current_variant["freq"][idx_expec]) )
        # TP, FP and FN
        current_variant["status"] = None
        if current_variant["freq"][idx_expec] > 0 and current_variant["freq"][idx_detec] > 0:
            current_variant["status"] = "TP"
            nb_TP += 1
        elif current_variant["freq"][idx_expec] > 0:
            current_variant["status"] = "FN"
            freq_FN.append( current_variant["freq"][idx_expec] )
            nb_FN += 1
        elif current_variant["freq"][idx_detec] > 0:
            current_variant["status"] = "FP"
            freq_FP.append( current_variant["freq"][idx_detec] )
            nb_FP += 1
        # Dissimilarity
        if current_variant["freq"][idx_expec] != 0:
            #~ expect_error_count_sum += abs(current_variant["error"])
            expect_error_count.append( abs(current_variant["error"]) )
            #~ expect_error_ratio_sum += current_variant["error_ratio"]
            expect_error_ratio.append( current_variant["error_ratio"] )
        # Error out of threshold
        if current_variant["error_ratio"] > error_threshold:
            nb_out_threshold += 1

    # Write
    with open(out_path, "w") as FH_out:
        FH_csv = csv.writer(FH_out, delimiter=separator)
        # Summary section
        FH_csv.writerow( ["[Summary]"] )
        FH_csv.writerow( ["#Nb_checked", "True_positive", "False_positive", "False_negative", "Errors_sum", "Errors_median", "Errors_ratio_sum", "Errors_ratio_median", "Error_out_of_threshold_(" + str(error_threshold) + ")", "Median_freq_FP", "Median_freq_FN"] )
        FH_csv.writerow([
            nb_checked,
            nb_TP,
            nb_FP,
            nb_FN,
            sum(expect_error_count),
            (0 if len(expect_error_count) == 0 else statistics.median(expect_error_count)),
            sum(expect_error_ratio),
            (0 if len(expect_error_ratio) == 0 else statistics.median(expect_error_ratio)),
            str(nb_out_threshold) + "/" + str(nb_checked),
            (0 if len(freq_FP) == 0 else statistics.median(freq_FP)),
            (0 if len(freq_FN) == 0 else statistics.median(freq_FN)) ])
        # Section separator
        FH_csv.writerow([])
        # Details section
        FH_csv.writerow( ["[Details]"] )
        FH_csv.writerow( ["#Chr:Pos", "Ref/Alt", "Type", "Expected", "Detected", "Status", "Error", "Error_ratio", "Out_of_threshold"] )
        for variant_id in sorted(variants):
            current_variant = variants[variant_id]
            FH_csv.writerow([
                current_variant["chrom"] + ":" + str(current_variant["pos"]),
                current_variant["ref"] + "/" + current_variant["alt"],
                current_variant["type"],
                current_variant["freq"][idx_expec],
                current_variant["freq"][idx_detec],
                current_variant["status"],
                current_variant["error"],
                current_variant["error_ratio"],
                (current_variant["error_ratio"] > error_threshold) ])

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
    parser.add_argument( '-t', '--error-threshold', type=float, default=0.8, help='The minimum similarity ratio between expected and detected. This value is used on TSV output to count the number of the expected variants with an error out of threshold. Simitarity ratio calculation: |1 - detected_freq/expected_freq| (e.g. ratio=0.5 for expected_freq=0.1 and detected_freq=0.05 ; ratio=0.25 for expected_freq=0.1 and detected_freq=0.4). [Default: %(default)s]' )
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
    if args.output_file.endswith("json"):
        writeJSONResults( all_variants, args.output_file )
    else:
        writeTSVResults( all_variants, args.output_file, args.error_threshold )
