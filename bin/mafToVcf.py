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
__version__ = '1.0.3'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import sys
import logging
import argparse
from anacore.maf import MAFIO, getName
from anacore.annotVcf import AnnotVCFIO, VCFRecord


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Convert MAF to VCF.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: MAF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))
    log.info("Version: " + str(__version__))

    # Get samples names and number of records by variant (several samples with the same variant)
    log.info("Get number of records by variant.")
    samples = set()
    occur_by_id = {}
    nb_records = 0
    with MAFIO(args.input_variants) as FH_in:
        for record in FH_in:
            nb_records += 1
            samples.add(record["Tumor_Sample_Barcode"])
            variants_id = getName(record)
            if variants_id not in occur_by_id:
                occur_by_id[variants_id] = {"nb_expec": 0, "data": {}}
            occur_by_id[variants_id]["nb_expec"] += 1
    samples = sorted(samples)
    log.info("{} samples, {} variants and {} records.".format(len(samples), len(occur_by_id), nb_records))

    # Convert
    log.info("Convert to VCF.")
    with MAFIO(args.input_variants) as FH_in:
        with AnnotVCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.samples = samples
            FH_out.ANN_titles = ["Allele", "Consequence", "SYMBOL", "Feature_type", "Feature", "HGVSc", "HGVSp", "RefSeq"]
            FH_out.info = {
                "SC": {"type": int, "type_tag": "Integer", "number": 1, "number_tag": "1", "description": "Mutated sample count"},
                "AD": {"type": int, "type_tag": "Integer", "number": 1, "number_tag": "1", "description": "Allele depth in tumor"},
                "AF": {"type": float, "type_tag": "Float", "number": 1, "number_tag": "1", "description": "Allele frequency in tumor"},
                "ANN": {"type": str, "type_tag": "String", "number": None, "number_tag": ".", "description": "Consequence annotations. Format: #TOREPLACE#"},
                "DP": {"type": int, "type_tag": "Integer", "number": 1, "description": "Depth in tumor"},
                "VarType": {"type": str, "type_tag": "String", "number": 1, "description": "Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)"}
            }
            FH_out.format = {
                "AD": {"type": int, "type_tag": "Integer", "number": 1, "number_tag": "1", "description": "Allele depth in tumor"},
                "AF": {"type": float, "type_tag": "Float", "number": 1, "number_tag": "1", "description": "Allele frequency in tumor"},
                "DP": {"type": int, "type_tag": "Integer", "number": 1, "description": "Depth in tumor"}
            }
            FH_out._writeHeader()
            # Records
            for record in FH_in:
                variants_id = getName(record)
                variant_by_spl = occur_by_id[variants_id]["data"]
                if int(record["t_depth"]) == 0:
                    log.warning(
                        "{} has a depth equal to 0 on line {} in file {}.".format(
                            variants_id,
                            FH_in.current_line_nb,
                            FH_in.filepath
                        )
                    )
                    variant_by_spl[record["Tumor_Sample_Barcode"]] = {
                        "AD": 0, "AF": 0, "DP": 0
                    }
                else:
                    variant_by_spl[record["Tumor_Sample_Barcode"]] = {
                        "AD": int(record["t_alt_count"]),
                        "AF": int(record["t_alt_count"]) / int(record["t_depth"]),
                        "DP": int(record["t_depth"]),
                    }
                if len(variant_by_spl) == occur_by_id[variants_id]["nb_expec"]:
                    for spl in FH_out.samples:
                        if spl not in variant_by_spl:
                            variant_by_spl[spl] = {"AD": None, "AF": None, "DP": None}
                    FH_out.write(
                        VCFRecord(
                            region=record["Chromosome"],
                            position=int(record["Start_Position"]),
                            knownSNPId=(record["dbSNP_RS"] if record["dbSNP_RS"] != "novel" and record["dbSNP_RS"] != "" else None),  # Empty is used for variants do not present in dbSNP but present in other db
                            refAllele=record["Reference_Allele"],
                            altAlleles=[record["Tumor_Seq_Allele2"]],
                            qual=(float(record["Score"]) if record["Score"] != "" else None),
                            pFilter=record["FILTER"].split(","),
                            info={
                                "SC": sum([1 for spl, var in variant_by_spl.items() if var["AD"] is not None and var["AD"] > 0]),
                                "AD": sum([var["AD"] for spl, var in variant_by_spl.items() if var["AD"] is not None]),
                                "AF": sum([var["AF"] for spl, var in variant_by_spl.items() if var["AF"] is not None]),
                                "DP": sum([var["DP"] for spl, var in variant_by_spl.items() if var["DP"] is not None]),
                                "VarType": record["Variant_Type"],
                                "ANN": [{
                                    "Allele": record["Allele"],
                                    "Consequence": record["Consequence"].replace(";", "&"),
                                    "SYMBOL": record["Hugo_Symbol"],
                                    "Feature_type": record["Feature_type"],
                                    "Feature": record["Feature"],
                                    "HGVSc": record["HGVSc"],
                                    "HGVSp": record["HGVSp"],
                                    "RefSeq": record["RefSeq"].replace(";", "&")
                                }]
                            },
                            pFormat=["AD", "AF", "DP"],
                            samples=variant_by_spl
                        )
                    )
                    occur_by_id[variants_id]["data"] = None

    log.info("End of job")
