import os
import sys
import json
import argparse

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from VEPvcf import VEPVCFIO, getAlleleRecord


if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='**********************************.' )
    group_filter = parser.add_argument_group( 'Filters' ) # Filters
    group_filter.add_argument( '-s', '--shared-threshold', default=2, help='Minimum number of samples with this variant. [Default: %(default)s]' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '--variants-files', nargs='+', required=True, help='The path to the file file containing variants and annotated with VEP v88+ (format: VCF).' )
    group_input.add_argument( '--ids-files', nargs='+', required=True, help='T*******************************.' )
    args = parser.parse_args()

    # Process
    kept_variants = dict()
    for current_ids in args.ids_files:
        with open(current_ids) as FH_ids:
            for line in FH_ids:
                ident = line.strip()
                if ident not in kept_variants:
                    kept_variants[ident] = 0
                kept_variants[ident] += 1
    removed_ids = list()
    for ident in kept_variants:
        if kept_variants[ident] < args.shared_threshold:
            removed_ids.append(ident)
    for ident in removed_ids:
        del(kept_variants[ident])


    for idx_vcf, current_vcf in enumerate(args.variants_files):
        with VEPVCFIO(current_vcf) as FH_vcf:
            if idx_vcf == 0:
                print(
                    "#ID",
                    "Chromosome",
                    "Position",
                    "Reference_allele",
                    "Alternative_allele",
                    "Variant_quality",
                    "Nb_share",
                    "\t".join( FH_vcf.CSQ_titles ),
                    sep="\t"    
                )
            for record in FH_vcf:
                for idx_alt, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord( FH_vcf, record, idx_alt )
                    alt_record.standardizeSingleAllele( "-" )
                    variant_id = alt_record.chrom + ":" + str(alt_record.pos) + "=" + alt_record.ref + "/" + alt
                    if variant_id in kept_variants:
                        for csq in alt_record.info["CSQ"]:
                            csq_values = list()
                            for title in FH_vcf.CSQ_titles:
                                if csq[title] is None:
                                    csq_values.append("")
                                else:
                                    csq_values.append( str(csq[title]) )
                            print(
                                variant_id,
                                alt_record.chrom,
                                alt_record.pos,
                                alt_record.ref,
                                alt_record.alt[0],
                                alt_record.qual,
                                kept_variants[variant_id],
                                "\t".join( csq_values ),
                                sep="\t"
                            )
                        del(kept_variants[variant_id])
    if len(kept_variants) != 0:
        raise Exception( "The following variant cannot be retrieved from VCF: " + ", ".join(kept_variants) )
