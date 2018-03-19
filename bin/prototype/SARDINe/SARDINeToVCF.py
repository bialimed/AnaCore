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
__version__ = '0.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'dev'

import argparse



########################################################################
#
# FUNCTIONS
#
########################################################################
def getVCFHeader( spl_name ):
    return "\n".join([
        '##fileformat=VCFv4.2',
        '##INFO=<ID=VT,Number=1,Type=String,Description="Variant Type : SNP(simple nucleotide polymorphism), INS(insertion), DEL(deletion)">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Mean Allele Frequency">',
        '##INFO=<ID=AM,Number=1,Type=String,Description="Amplicon">',
        '##INFO=<ID=NM,Number=1,Type=String,Description="Genomique Sequence Version">',
        '##INFO=<ID=LI,Number=1,Type=String,Description="Library : A, B, AB">',
        '##FILTER=<ID=LI,Description="Variants detected on a single library ">',
        '##FORMAT=<ID=DPFa1,Number=1,Type=Integer,Description="Depth of forward reads supporting bases on the first library">',
        '##FORMAT=<ID=AFFa1,Number=1,Type=Float,Description="Allele frequency of variants on the first library for the forward reads">',
        '##FORMAT=<ID=DPRa1,Number=1,Type=Integer,Description="Depth of reverse reads supporting bases on the first library">',
        '##FORMAT=<ID=AFRa1,Number=1,Type=Float,Description="Allele frequency of variants on the first library for the reverse reads">',
        '##FORMAT=<ID=DPFa2,Number=1,Type=Integer,Description="Depth of forward reads supporting bases on the second library">',
        '##FORMAT=<ID=AFFa2,Number=1,Type=Float,Description="Allele frequency of variants on the second library for the forward reads">',
        '##FORMAT=<ID=DPRa2,Number=1,Type=Integer,Description="Depth of reverse reads supporting bases on the second library">',
        '##FORMAT=<ID=AFRa2,Number=1,Type=Float,Description="Allele frequency of variants on the second library for the reverse reads">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + spl_name
    ])



########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Transforms SARDINe CSV annotation file in VCF.' )
    parser.add_argument( '-s', '--sample-name', required=True, help='The sample name.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-annot', required=True, help="The path to the SARDINe's annotation file (format: CSV)." )
    args = parser.parse_args()

    # Header
    print( getVCFHeader(args.sample_name) )

    # Variants
    with open(args.input_annot) as FH_in:
        for line in FH_in:
            if not line.startswith('#'):
                fields = [elt.strip() for elt in line.split(";")]
                gene = fields[0]
                chrom, exon, start, stop, ref, alt, NM, nt_change, aa_change, library, amplicon_gp, amplicon_name, DP_tot, FA_mean, DP_ar1, FA_ar1, DP_ar2, FA_ar2, DP_br1, FA_br1, DP_br2, FA_br2 = fields[-22:]
                function = " ".join(fields[1:-22])
                filter = ("PASS" if library == "AB" else "OneLibrary")
                VT = "SNP"
                if ref == "-":
                    VT = "INS"
                elif alt == "-":
                    VT = "DEL"
                VCF_fiels = [ 
                    chrom,
                    str(start),
                    '.',
                    ref,
                    alt,
                    '.',
                    filter,
                    "VT=" + VT + ";DP=" + DP_tot + ";AF=" + str(FA_mean) + ";AM=" + amplicon_name + ";NM=" + NM + ";LI=" + library,
                    "AFFa1:DPFa1:AFRa1:DPRa1:AFFa2:DPFa2:AFRa2:DPRa2",
                    FA_ar1 + ":" + DP_ar1 + ":" + FA_ar2 + ":" + DP_ar2 + ":" + FA_br1 + ":" + DP_br1 + ":" + FA_br2 + ":" + DP_br2    ]
                print( "\t".join(VCF_fiels) )
